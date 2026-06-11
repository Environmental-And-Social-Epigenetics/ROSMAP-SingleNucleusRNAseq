#!/bin/bash
#SBATCH -p pi_lhtsai,pi_manoli
#SBATCH -n 16
#SBATCH --mem=300G
#SBATCH -t 24:00:00
#SBATCH -o %j_wgcna.out
#SBATCH -e %j_wgcna.err

# =============================================================================
# SLURM Job: ACE hdWGCNA Co-Expression Network Analysis
#
# Runs hdWGCNA on one cell type × sex combination.
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

INTEGRATION="${1:?ERROR: integration required}"
PHENOTYPE="${2:?ERROR: phenotype required}"
SEX="${3:?ERROR: sex required (Male or Female)}"
CELL_TYPE="${4:?ERROR: cell type required}"
INPUT_H5AD="${5:?ERROR: input h5ad required}"
OUTPUT_ROOT="${6:?ERROR: output root required}"

RESULTS_DIR="${OUTPUT_ROOT}/results_${INTEGRATION}/${PHENOTYPE}/${SEX}_${CELL_TYPE}"
mkdir -p "${RESULTS_DIR}"

echo "=== hdWGCNA Analysis ==="
echo "Integration: ${INTEGRATION}"
echo "Phenotype:   ${PHENOTYPE}"
echo "Sex:         ${SEX}"
echo "Cell type:   ${CELL_TYPE}"
echo "Input:       ${INPUT_H5AD}"
echo "Output:      ${RESULTS_DIR}"
echo ""

export HDF5_USE_FILE_LOCKING=FALSE

# Use env's Rscript directly to bypass any login-shell PATH pollution
"${WGCNA_ENV}/bin/Rscript" "${SCRIPT_DIR}/wgcna_analysis.R" \
    --cell-type "${CELL_TYPE}" \
    --sex "${SEX}" \
    --phenotype "${PHENOTYPE}" \
    --input-h5ad "${INPUT_H5AD}" \
    --pheno-csv "${ACE_SCORES_CSV}" \
    --output-dir "${RESULTS_DIR}" \
    --deg-results-dir "${ACE_OUTPUT_ROOT}/DEG/Tsai/results_${INTEGRATION}/${PHENOTYPE}"
