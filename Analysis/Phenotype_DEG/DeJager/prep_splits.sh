#!/bin/bash
# prep_splits.sh -- Phase D1 wrapper: prep per-cell-type raw-count h5ad files
# for one integration. Uses ACE/DEG/DeJager/prep_celltype_splits.py.
#
# Usage (from sbatch):
#   sbatch -p pi_lhtsai -n 16 --mem=350G -t 12:00:00 prep_splits.sh <integration>

#SBATCH -p ou_bcs_normal
#SBATCH -n 16
#SBATCH --mem=350G
#SBATCH -t 12:00:00
#SBATCH -o %j.out
#SBATCH -e %j.err

set -euo pipefail

if [[ -z "${REPO_ROOT:-}" ]]; then
  if git -C "${PWD}" rev-parse --show-toplevel >/dev/null 2>&1; then
    REPO_ROOT="$(git -C "${PWD}" rev-parse --show-toplevel)"
  else
    echo "ERROR: REPO_ROOT is not set and could not be inferred from the current directory."
    exit 1
  fi
fi
SCRIPT_DIR="${REPO_ROOT}/Analysis/Phenotype_DEG/DeJager"
source "${REPO_ROOT}/config/paths.sh"

INTEGRATION="${1:?ERROR: integration argument required (library_id|patient_id|pool_batch|derived_batch|sequencing_date)}"

OUTPUT_ROOT="${ANALYSIS_OUTPUT_ROOT:-${WORKSPACE_ROOT}/Analysis_Outputs}/Phenotype_DEG/DeJager"
SPLIT_DIR="${OUTPUT_ROOT}/celltype_splits_${INTEGRATION}"
mkdir -p "${SPLIT_DIR}"

# Reuse ACE prep script (already supports all 5 integrations now)
ACE_PREP="${REPO_ROOT}/Analysis/ACE/DEG/DeJager/prep_celltype_splits.py"

set +u
source "${CONDA_INIT_SCRIPT}"
conda activate "${NEBULA_ENV}"
set -u
export HDF5_USE_FILE_LOCKING=FALSE
export PYTHONUNBUFFERED=1

echo "=== Phenotype DEG prep: ${INTEGRATION} ==="
echo "Output dir: ${SPLIT_DIR}"

python -u "${ACE_PREP}" \
  --integration "${INTEGRATION}" \
  --output-dir "${SPLIT_DIR}"

echo "Done."
