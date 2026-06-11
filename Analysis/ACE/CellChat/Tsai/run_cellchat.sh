#!/bin/bash
#SBATCH -p pi_lhtsai,pi_manoli
#SBATCH -n 16
#SBATCH --mem=200G
#SBATCH -t 12:00:00
#SBATCH -o %j_cellchat.out
#SBATCH -e %j_cellchat.err

# =============================================================================
# SLURM Job: ACE CellChat Differential Interaction Analysis
#
# Runs CellChat on ACE-high vs ACE-low groups for a given sex.
# =============================================================================

set -euo pipefail

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

INTEGRATION="${1:?ERROR: integration argument required}"
PHENOTYPE="${2:?ERROR: phenotype argument required}"
SEX="${3:?ERROR: sex argument required (Male or Female)}"
INPUT_H5AD="${4:?ERROR: input h5ad path required}"
OUTPUT_ROOT="${5:?ERROR: output root directory required}"

RESULTS_DIR="${OUTPUT_ROOT}/results_${INTEGRATION}/${PHENOTYPE}/${SEX}"
mkdir -p "${RESULTS_DIR}"

echo "=== CellChat Analysis ==="
echo "Integration: ${INTEGRATION}"
echo "Phenotype:   ${PHENOTYPE}"
echo "Sex:         ${SEX}"
echo "Input:       ${INPUT_H5AD}"
echo "Output:      ${RESULTS_DIR}"
echo ""

export HDF5_USE_FILE_LOCKING=FALSE

# Use env's Rscript directly to bypass any login-shell PATH pollution
"${CELLCHAT_ENV}/bin/Rscript" "${SCRIPT_DIR}/cellchat_analysis.R" \
    --sex "${SEX}" \
    --phenotype "${PHENOTYPE}" \
    --input-h5ad "${INPUT_H5AD}" \
    --pheno-csv "${ACE_SCORES_CSV}" \
    --output-dir "${RESULTS_DIR}"
