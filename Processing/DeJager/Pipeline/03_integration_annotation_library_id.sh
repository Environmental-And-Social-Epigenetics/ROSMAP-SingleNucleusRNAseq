#!/bin/bash
#SBATCH -J dej_integrate_libid
#SBATCH -t 48:00:00
#SBATCH -n 32
#SBATCH --mem=500G
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

INPUT_DIR="${DEJAGER_DOUBLET_REMOVED}"
OUTPUT_DIR="${DEJAGER_INTEGRATED_LIBRARY_ID}"
MARKERS_RDS="${DEJAGER_MARKERS_RDS}"

python "${SCRIPT_DIR}/03_integration_annotation.py" \
    --input-dir "${INPUT_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --markers-rds "${MARKERS_RDS}" \
    --harmony-batch-key library_id
