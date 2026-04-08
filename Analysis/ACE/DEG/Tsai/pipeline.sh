#!/bin/bash
# pipeline.sh -- End-to-end ACE DEG analysis pipeline with AI-generated report
#
# Chains all stages via SLURM dependencies:
#   prep -> DEG -> summary -> figures + GO enrichment -> report (Claude agent)
#
# Usage:
#   bash pipeline.sh --integration derived_batch --author "Mahmoud Abdelmoneum"
#   bash pipeline.sh --integration projid --stage summary   # run single stage
#   bash pipeline.sh --integration derived_batch --dry-run   # show plan only
#
# Stages: prep, deg, summary, figures, go, report, all (default)

set -euo pipefail

# ── Resolve paths ─────────────────────────────────────────────────────────────
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  SCRIPT_DIR="${SLURM_SUBMIT_DIR}"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"

# Source centralized config
source "${REPO_ROOT}/config/paths.sh"
source "${SCRIPT_DIR}/pipeline.conf"

# ── Parse arguments ──────────────────────────────────────────────────────────
INTEGRATION=""
STAGE="all"
DRY_RUN=false

while [[ $# -gt 0 ]]; do
  case "$1" in
    --integration) INTEGRATION="$2"; shift 2 ;;
    --stage)       STAGE="$2"; shift 2 ;;
    --author)      AUTHOR="$2"; shift 2 ;;
    --dry-run)     DRY_RUN=true; shift ;;
    *) echo "Unknown argument: $1"; exit 1 ;;
  esac
done

if [[ -z "$INTEGRATION" ]]; then
  echo "ERROR: --integration required (derived_batch or projid)"
  exit 1
fi

# Normalize integration name
case "$INTEGRATION" in
  batch|derived_batch) INTEGRATION="derived_batch" ;;
  projid) ;;
  *) echo "ERROR: Unknown integration '$INTEGRATION'"; exit 1 ;;
esac

# ── Output paths ──────────────────────────────────────────────────────────────
OUTPUT_ROOT="${ANALYSIS_OUTPUT_ROOT:-${SCRIPT_DIR}}/ACE/DEG/Tsai"
# Fall back to local results dirs if ANALYSIS_OUTPUT_ROOT is not set
if [[ ! -d "$(dirname "$OUTPUT_ROOT")" ]]; then
  OUTPUT_ROOT="${SCRIPT_DIR}"
fi

RESULTS_DIR="${SCRIPT_DIR}/results_${INTEGRATION}"
SPLIT_DIR="${SCRIPT_DIR}/celltype_splits_${INTEGRATION}"
FIGURES_DIR="${RESULTS_DIR}/figures"
GO_DIR="${RESULTS_DIR}/go_enrichment"
REPORT_DIR="${SCRIPT_DIR}/report"
LOG_DIR="${SCRIPT_DIR}/logs"
MANIFEST="${SCRIPT_DIR}/pipeline_manifest_${INTEGRATION}.json"

mkdir -p "$LOG_DIR" "$REPORT_DIR"

CONDA_SETUP="source ${CONDA_INIT_SCRIPT} && conda activate ${NEBULA_ENV} && export HDF5_USE_FILE_LOCKING=FALSE"
PYTHON="${NEBULA_ENV}/bin/python3"

echo "=== ACE DEG Pipeline ==="
echo "Integration:  $INTEGRATION"
echo "Stage:        $STAGE"
echo "Author:       $AUTHOR"
echo "Results dir:  $RESULTS_DIR"
echo "Partitions:   $SLURM_PARTITION_SMALL"
echo "Dry run:      $DRY_RUN"
echo ""

# ── Helper: submit or print ──────────────────────────────────────────────────
submit_job() {
  if $DRY_RUN; then
    echo "[DRY RUN] sbatch $*"
    echo "DRY_RUN_JOB_$(date +%s%N)"
  else
    sbatch --parsable "$@"
  fi
}

should_run() {
  local stage="$1"
  [[ "$STAGE" == "all" || "$STAGE" == "$stage" ]]
}

# ── Track job IDs for dependencies ────────────────────────────────────────────
declare -a DEG_JOBS=()
PREP_JOB=""
SUMMARY_JOB=""
FIGURES_JOB=""
GO_JOB=""

# ══════════════════════════════════════════════════════════════════════════════
# STAGE 1: PREP (cell type splitting)
# ══════════════════════════════════════════════════════════════════════════════
if should_run "prep"; then
  if [[ -d "$SPLIT_DIR" ]] && ls "$SPLIT_DIR"/*.h5ad &>/dev/null; then
    echo "[PREP] Cell type splits already exist in $SPLIT_DIR, skipping."
    echo "       Delete the directory to force re-run."
  else
    echo "[PREP] Submitting cell type splitting..."
    PREP_JOB=$(submit_job \
      -p "$SLURM_PARTITION_LARGE" \
      -n 16 --mem=200G -t 24:00:00 \
      -o "${LOG_DIR}/%j_prep_${INTEGRATION}.out" \
      -e "${LOG_DIR}/%j_prep_${INTEGRATION}.err" \
      --job-name="ace_prep_${INTEGRATION}" \
      --wrap="${CONDA_SETUP} && ${PYTHON} ${SCRIPT_DIR}/prep_celltype_splits.py --integration ${INTEGRATION} --output-dir ${SPLIT_DIR}")
    echo "  Job: $PREP_JOB"
  fi
fi

# ══════════════════════════════════════════════════════════════════════════════
# STAGE 2: DEG (per phenotype)
# ══════════════════════════════════════════════════════════════════════════════
if should_run "deg"; then
  echo "[DEG] Submitting DEG analysis jobs..."
  DEP_FLAG=""
  [[ -n "$PREP_JOB" ]] && DEP_FLAG="--dependency=afterok:${PREP_JOB}"

  for PHENO in $PHENOTYPES; do
    JOB=$(submit_job \
      -p "$SLURM_PARTITION_MEDIUM" \
      -t 16:00:00 --mem=100G \
      $DEP_FLAG \
      -o "${LOG_DIR}/%j_deg_${PHENO}_${INTEGRATION}.out" \
      -e "${LOG_DIR}/%j_deg_${PHENO}_${INTEGRATION}.err" \
      --job-name="ace_deg_${PHENO}_${INTEGRATION}" \
      "${SCRIPT_DIR}/run_deg.sh" "${INTEGRATION}" "${PHENO}")
    DEG_JOBS+=("$JOB")
    echo "  ${PHENO}: $JOB"
  done
fi

# Build dependency string for post-DEG stages
DEG_DEP=""
if [[ ${#DEG_JOBS[@]} -gt 0 ]]; then
  DEG_DEP="--dependency=afterok:$(IFS=:; echo "${DEG_JOBS[*]}")"
fi

# ══════════════════════════════════════════════════════════════════════════════
# STAGE 3: SUMMARIZE
# ══════════════════════════════════════════════════════════════════════════════
if should_run "summary"; then
  echo "[SUMMARY] Submitting result summarization..."
  SUMMARY_JOB=$(submit_job \
    -p "$SLURM_PARTITION_SMALL" \
    -t 1:00:00 --mem=16G -n 2 \
    $DEG_DEP \
    -o "${LOG_DIR}/%j_summary_${INTEGRATION}.out" \
    -e "${LOG_DIR}/%j_summary_${INTEGRATION}.err" \
    --job-name="ace_summary_${INTEGRATION}" \
    --wrap="${CONDA_SETUP} && Rscript ${SCRIPT_DIR}/summarize_results.R --results-dir ${RESULTS_DIR} --output-dir ${RESULTS_DIR}")
  echo "  Job: $SUMMARY_JOB"
fi

SUMMARY_DEP=""
[[ -n "${SUMMARY_JOB:-}" ]] && SUMMARY_DEP="--dependency=afterok:${SUMMARY_JOB}"

# ══════════════════════════════════════════════════════════════════════════════
# STAGE 4: FIGURES (parallel with GO)
# ══════════════════════════════════════════════════════════════════════════════
if should_run "figures"; then
  echo "[FIGURES] Submitting figure generation..."
  FIGURES_JOB=$(submit_job \
    -p "$SLURM_PARTITION_SMALL" \
    -t 2:00:00 --mem=32G -n 4 \
    $SUMMARY_DEP \
    -o "${LOG_DIR}/%j_figures_${INTEGRATION}.out" \
    -e "${LOG_DIR}/%j_figures_${INTEGRATION}.err" \
    --job-name="ace_figures_${INTEGRATION}" \
    --wrap="${CONDA_SETUP} && Rscript ${SCRIPT_DIR}/generate_figures.R --results-dir ${RESULTS_DIR} --summary-csv ${RESULTS_DIR}/deg_summary.csv --output-dir ${FIGURES_DIR} --integration-label ${INTEGRATION}")
  echo "  Job: $FIGURES_JOB"
fi

# ══════════════════════════════════════════════════════════════════════════════
# STAGE 5: GO ENRICHMENT (parallel with figures)
# ══════════════════════════════════════════════════════════════════════════════
if should_run "go"; then
  echo "[GO] Submitting GO enrichment..."
  GO_JOB=$(submit_job \
    -p "$SLURM_PARTITION_SMALL" \
    -t 4:00:00 --mem=32G -n 4 \
    $SUMMARY_DEP \
    -o "${LOG_DIR}/%j_go_${INTEGRATION}.out" \
    -e "${LOG_DIR}/%j_go_${INTEGRATION}.err" \
    --job-name="ace_go_${INTEGRATION}" \
    --wrap="${CONDA_SETUP} && Rscript ${SCRIPT_DIR}/run_go_enrichment.R --results-dir ${RESULTS_DIR} --summary-json ${RESULTS_DIR}/deg_summary.json --output-dir ${GO_DIR} --min-degs ${GO_MIN_DEGS} --top-n ${GO_TOP_N}")
  echo "  Job: $GO_JOB"
fi

# ══════════════════════════════════════════════════════════════════════════════
# STAGE 6: REPORT (Claude agent, runs on login node after all SLURM jobs)
# ══════════════════════════════════════════════════════════════════════════════
if should_run "report"; then
  # Build dependency on both figures and GO
  REPORT_DEPS=""
  REPORT_DEP_JOBS=()
  [[ -n "${FIGURES_JOB:-}" ]] && REPORT_DEP_JOBS+=("$FIGURES_JOB")
  [[ -n "${GO_JOB:-}" ]] && REPORT_DEP_JOBS+=("$GO_JOB")

  if [[ ${#REPORT_DEP_JOBS[@]} -gt 0 ]]; then
    REPORT_DEPS="--dependency=afterok:$(IFS=:; echo "${REPORT_DEP_JOBS[*]}")"
    echo "[REPORT] Will run after figures + GO complete."
    echo "  Dependencies: ${REPORT_DEP_JOBS[*]}"

    if ! $DRY_RUN; then
      echo "  Submitting report generation as SLURM job (login-node compatible)..."
      # Use srun to wait for dependencies, then run on login node
      # Or submit as a lightweight SLURM job
      submit_job \
        -p "$SLURM_PARTITION_SMALL" \
        -t 1:00:00 --mem=8G -n 1 \
        $REPORT_DEPS \
        -o "${LOG_DIR}/%j_report_${INTEGRATION}.out" \
        -e "${LOG_DIR}/%j_report_${INTEGRATION}.err" \
        --job-name="ace_report_${INTEGRATION}" \
        --wrap="bash ${SCRIPT_DIR}/generate_report.sh --results-dir ${RESULTS_DIR} --figures-dir ${FIGURES_DIR} --go-dir ${GO_DIR} --integration '${INTEGRATION}' --author '${AUTHOR}' --affiliation '${AFFILIATION}' --output-dir ${REPORT_DIR}"
    fi
  else
    echo "[REPORT] Running report generation now (no SLURM dependencies)..."
    if ! $DRY_RUN; then
      bash "${SCRIPT_DIR}/generate_report.sh" \
        --results-dir "$RESULTS_DIR" \
        --figures-dir "$FIGURES_DIR" \
        --go-dir "$GO_DIR" \
        --integration "$INTEGRATION" \
        --author "$AUTHOR" \
        --affiliation "$AFFILIATION" \
        --output-dir "$REPORT_DIR"
    else
      echo "[DRY RUN] Would run generate_report.sh"
    fi
  fi
fi

# ── Save manifest ─────────────────────────────────────────────────────────────
if ! $DRY_RUN; then
  cat > "$MANIFEST" <<EOF
{
  "pipeline": "ACE DEG Analysis",
  "integration": "$INTEGRATION",
  "timestamp": "$(date -Iseconds)",
  "author": "$AUTHOR",
  "stage": "$STAGE",
  "jobs": {
    "prep": "${PREP_JOB:-null}",
    "deg": [$(IFS=,; echo "\"${DEG_JOBS[*]:-}\"" | sed 's/,/","/g')],
    "summary": "${SUMMARY_JOB:-null}",
    "figures": "${FIGURES_JOB:-null}",
    "go": "${GO_JOB:-null}"
  },
  "paths": {
    "results": "$RESULTS_DIR",
    "figures": "$FIGURES_DIR",
    "go": "$GO_DIR",
    "report": "$REPORT_DIR"
  }
}
EOF
  echo ""
  echo "Manifest saved: $MANIFEST"
fi

echo ""
echo "=== Pipeline submitted ==="
echo "Monitor with: squeue -u \$USER"
echo "Results will be in: $RESULTS_DIR"
echo "Report will be in:  $REPORT_DIR"
