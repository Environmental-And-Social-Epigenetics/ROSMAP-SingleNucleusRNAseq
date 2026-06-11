#!/bin/bash
#SBATCH -p pi_lhtsai,pi_manoli,ou_bcs_low,mit_normal
#SBATCH -n 8
#SBATCH --mem=128G
#SBATCH -t 6:00:00
#SBATCH -o %j_micstate.out
#SBATCH -e %j_micstate.err

# =============================================================================
# SLURM Job: ACE Microglial State Analysis
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
INPUT_H5AD="${4:?ERROR: input h5ad required}"
OUTPUT_ROOT="${5:?ERROR: output root required}"

RESULTS_DIR="${OUTPUT_ROOT}/results_${INTEGRATION}/${PHENOTYPE}/${SEX}"
mkdir -p "${RESULTS_DIR}"

echo "=== Microglial State Analysis ==="
echo "Integration: ${INTEGRATION}"
echo "Phenotype:   ${PHENOTYPE}"
echo "Sex:         ${SEX}"
echo "Input:       ${INPUT_H5AD}"
echo "Output:      ${RESULTS_DIR}"
echo ""

set +u
activate_env "${DECOUPLER_ENV}"
set -u
export HDF5_USE_FILE_LOCKING=FALSE

"${DECOUPLER_ENV}/bin/python" "${SCRIPT_DIR}/mic_state_analysis.py" \
    --sex "${SEX}" \
    --phenotype "${PHENOTYPE}" \
    --input-h5ad "${INPUT_H5AD}" \
    --pheno-csv "${ACE_SCORES_CSV}" \
    --output-dir "${RESULTS_DIR}"
