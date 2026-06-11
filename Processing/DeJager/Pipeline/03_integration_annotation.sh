#!/bin/bash
#SBATCH -J dej_integrate
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

# NOTE: SLURM log redirection is handled by submit_pipeline.sh via
# command-line --output/--error flags.
# Stage 3 is a single integration job (not a SLURM array) because it
# must load and concatenate all samples simultaneously.
mkdir -p "${DEJAGER_PROCESSING_LOGS}"

export HDF5_USE_FILE_LOCKING=FALSE
# Disable numba JIT to avoid CPUDispatcher crash on Python 3.14
export NUMBA_DISABLE_JIT=1
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
# Ensure the activated env's bin is first on PATH (overrides any stale
# conda env directories that may be baked into the user's shell profile).
export PATH="${CONDA_PREFIX}/bin:${PATH}"

python -m rosmap_tx.processing \
    --dataset dejager \
    --stage 3 \
    --variant library_id
