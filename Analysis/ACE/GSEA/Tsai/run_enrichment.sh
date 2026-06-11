#!/bin/bash
#SBATCH -p pi_lhtsai,pi_manoli,ou_bcs_low,mit_normal
#SBATCH -n 8
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH -o gsea_%j.out
#SBATCH -e gsea_%j.err

set -euo pipefail

# Under SLURM, BASH_SOURCE points to a temp copy in /var/spool/slurmd/;
# use SLURM_SUBMIT_DIR (the directory where sbatch was run) instead.
# Prefer LAUNCHER_SCRIPT_DIR (exported from launcher) over SLURM_SUBMIT_DIR
if [[ -n "${LAUNCHER_SCRIPT_DIR:-}" ]]; then
  SCRIPT_DIR="${LAUNCHER_SCRIPT_DIR}"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  SCRIPT_DIR="${SLURM_SUBMIT_DIR}"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

# ── Validate arguments ──────────────────────────────────────────────────────

INTEGRATION="${1:?ERROR: arg 1 (integration) is required (e.g. derived_batch)}"
PHENOTYPE="${2:?ERROR: arg 2 (phenotype) is required (e.g. tot_adverse_exp)}"
SEX="${3:?ERROR: arg 3 (sex) is required (Female or Male)}"
DEG_ROOT="${4:?ERROR: arg 4 (DEG results root) is required}"
OUTPUT_DIR="${5:?ERROR: arg 5 (output directory) is required}"
# Optional arg 6: male AD-model arm (e.g. MaleContAD). When set, DEG_ROOT is the
# per-arm DEG dir (results_derived_batch_<ARM>_AllCellTypes) and the analysis reads
# the arm-suffixed .rda files instead of the sex-stratified ones.
ARM="${6:-}"

# ── Create output directory ──────────────────────────────────────────────────

if [[ -n "${ARM}" ]]; then
  RESULTS_DIR="${OUTPUT_DIR}/results_${INTEGRATION}_${ARM}/${PHENOTYPE}"
else
  RESULTS_DIR="${OUTPUT_DIR}/results_${INTEGRATION}/${PHENOTYPE}/${SEX}"
fi
mkdir -p "${RESULTS_DIR}"

echo "=== ACE GSEA Job ==="
echo "Integration: ${INTEGRATION}"
echo "Phenotype:   ${PHENOTYPE}"
echo "Sex:         ${SEX}"
echo "Arm:         ${ARM:-(baseline)}"
echo "DEG root:    ${DEG_ROOT}"
echo "Results dir: ${RESULTS_DIR}"
echo ""

# ── Activate environment ─────────────────────────────────────────────────────

set +u
activate_env "${GSEA_ANALYSIS_ENV}"
set -u

# ── Run GSEA analysis ────────────────────────────────────────────────────────

GSEA_ARGS=(
  --deg-results-dir "${DEG_ROOT}/${PHENOTYPE}"
  --phenotype "${PHENOTYPE}"
  --sex "${SEX}"
  --output-dir "${RESULTS_DIR}"
)
if [[ -n "${ARM}" ]]; then
  GSEA_ARGS+=(--arm "${ARM}")
fi
if [[ -n "${GSEA_CELLTYPES:-}" ]]; then
  GSEA_ARGS+=(--celltypes "${GSEA_CELLTYPES}")
fi
if [[ "${SMOKE_FLAG:-}" == "--smoke" ]]; then
  GSEA_ARGS+=(--smoke)
fi

"${GSEA_ANALYSIS_ENV}/bin/Rscript" "${SCRIPT_DIR}/gsea_analysis.R" "${GSEA_ARGS[@]}"

echo ""
echo "GSEA job complete for ${PHENOTYPE} / ${SEX} / ${ARM:-baseline}"
