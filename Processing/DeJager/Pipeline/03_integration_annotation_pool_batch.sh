#!/bin/bash
#SBATCH -J dej_integrate_pool
#SBATCH -t 48:00:00
#SBATCH -n 32
#SBATCH --mem=350G
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

# When submitted via submit_pipeline.sh, REPO_ROOT is already exported.
# Fall back to deriving from BASH_SOURCE for standalone use.
if [[ -z "${REPO_ROOT:-}" ]]; then
    _SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    REPO_ROOT="$(cd "${_SCRIPT_DIR}/../../.." && pwd)"
fi
SCRIPT_DIR="${REPO_ROOT}/Processing/DeJager/Pipeline"

source "${REPO_ROOT}/config/paths.sh"

mkdir -p "${DEJAGER_PROCESSING_LOGS}"

export HDF5_USE_FILE_LOCKING=FALSE
export PYTHONUNBUFFERED=1
export PYTHONPATH="${REPO_ROOT}/src:${PYTHONPATH:-}"

# Temporarily relax nounset for conda activation (some activate.d scripts
# reference unset variables like ADDR2LINE).
set +u
source "${CONDA_INIT_SCRIPT}"
conda activate "${BATCHCORR_ENV}"
set -u
if [[ -z "${CONDA_PREFIX:-}" ]]; then
    echo "ERROR: Failed to activate conda environment: ${BATCHCORR_ENV}"
    exit 1
fi
export PATH="${CONDA_PREFIX}/bin:${PATH}"

python -m rosmap_tx.processing \
    --dataset dejager \
    --stage 3 \
    --variant pool_batch
