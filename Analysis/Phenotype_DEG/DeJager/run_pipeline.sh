#!/bin/bash
# run_pipeline.sh -- Submit the full Phenotype DEG comparison pipeline.
#
# Phase D1: 5 prep jobs (one per integration), each 2-4h on pi_lhtsai
# Phase D2: 60 DESeq2 jobs (5 integrations × 2 phenotypes × 6 cell types),
#           each ~30 min, can run in parallel after D1 for that integration
# Phase D3: 1 aggregator job, depends on all D2 jobs

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

OUTPUT_ROOT="${ANALYSIS_OUTPUT_ROOT:-${WORKSPACE_ROOT}/Analysis_Outputs}/Phenotype_DEG/DeJager"
LOG_DIR="${OUTPUT_ROOT}/logs"
mkdir -p "${LOG_DIR}"

INTEGRATIONS=(library_id patient_id pool_batch derived_batch sequencing_date)
PHENOTYPES=(msex cogdx_binary)
CELL_TYPES=(Exc Inh Ast Mic Oli OPC)

echo "=== Phenotype DEG Pipeline ==="
echo "Integrations: ${INTEGRATIONS[*]}"
echo "Phenotypes:   ${PHENOTYPES[*]}"
echo "Cell types:   ${CELL_TYPES[*]}"
echo "Output root:  ${OUTPUT_ROOT}"

ALL_DEG_JOBS=()

for INTEGRATION in "${INTEGRATIONS[@]}"; do
  echo ""
  echo "--- Integration: ${INTEGRATION} ---"

  # ----- Phase D1: prep -----
  PREP_OUT="${OUTPUT_ROOT}/celltype_splits_${INTEGRATION}"
  if [[ -d "${PREP_OUT}" && $(ls -1 "${PREP_OUT}"/*.h5ad 2>/dev/null | wc -l) -ge 6 ]]; then
    echo "  Prep already complete (${PREP_OUT} has h5ad files). Skipping D1."
    PREP_DEP=""
  else
    PREP_JOB=$(sbatch --parsable \
      -o "${LOG_DIR}/%j_prep_${INTEGRATION}.out" \
      -e "${LOG_DIR}/%j_prep_${INTEGRATION}.err" \
      --job-name="phDEG_prep_${INTEGRATION}" \
      "${SCRIPT_DIR}/prep_splits.sh" "${INTEGRATION}")
    echo "  Phase D1 prep: job ${PREP_JOB}"
    PREP_DEP="--dependency=afterok:${PREP_JOB}"
  fi

  # ----- Phase D2: per (phenotype, celltype) DEG -----
  for PHENO in "${PHENOTYPES[@]}"; do
    for CT in "${CELL_TYPES[@]}"; do
      DEG_JOB=$(sbatch --parsable ${PREP_DEP} \
        -o "${LOG_DIR}/%j_${PHENO}_${CT}_${INTEGRATION}.out" \
        -e "${LOG_DIR}/%j_${PHENO}_${CT}_${INTEGRATION}.err" \
        --job-name="phDEG_${PHENO}_${CT}_${INTEGRATION}" \
        "${SCRIPT_DIR}/phenoDeg.sh" "${INTEGRATION}" "${PHENO}" "${CT}")
      echo "    D2 ${PHENO}/${CT}/${INTEGRATION}: job ${DEG_JOB}"
      ALL_DEG_JOBS+=("${DEG_JOB}")
    done
  done
done

# ----- Phase D3: aggregator depends on all D2 jobs -----
DEP_LIST=$(IFS=:; echo "${ALL_DEG_JOBS[*]}")
echo ""
echo "--- Phase D3: aggregator ---"
AGG_JOB=$(sbatch --parsable \
  -p ou_bcs_normal,mit_normal \
  -n 4 \
  --mem=32G \
  -t 02:00:00 \
  --dependency=afterany:${DEP_LIST} \
  -o "${LOG_DIR}/%j_aggregate.out" \
  -e "${LOG_DIR}/%j_aggregate.err" \
  --job-name="phDEG_aggregate" \
  --wrap="set -euo pipefail; source ${REPO_ROOT}/config/paths.sh; set +u; source \${CONDA_INIT_SCRIPT}; conda activate \${NEBULA_ENV}; set -u; export PYTHONUNBUFFERED=1; python -u ${SCRIPT_DIR}/aggregate_effects.py --results-root ${OUTPUT_ROOT} --output-dir ${OUTPUT_ROOT}/comparison")

echo "Phase D3 aggregator: job ${AGG_JOB}"
echo ""
echo "All jobs submitted. ${#ALL_DEG_JOBS[@]} DEG jobs + 1 aggregator."
echo "Monitor with: squeue -u \$USER"
