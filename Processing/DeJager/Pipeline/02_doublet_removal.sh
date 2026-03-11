#!/bin/bash
#SBATCH -J dej_doublets
#SBATCH -t 12:00:00
#SBATCH -n 4
#SBATCH --mem=32G
#SBATCH --array=1-200%32
#SBATCH --mail-type=FAIL

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
mkdir -p "${DEJAGER_PROCESSING_LOGS}"

export HDF5_USE_FILE_LOCKING=FALSE

# Temporarily relax nounset for conda activation (some activate.d scripts
# reference unset variables like ADDR2LINE).
set +u
source "${CONDA_INIT_SCRIPT}"
conda activate "${SINGLECELL_ENV}"
set -u
if [[ -z "${CONDA_PREFIX:-}" ]]; then
    echo "ERROR: Failed to activate conda environment: ${SINGLECELL_ENV}"
    exit 1
fi

INPUT_DIR="${DEJAGER_QC_FILTERED}"
OUTPUT_DIR="${DEJAGER_DOUBLET_REMOVED}"

mapfile -t SAMPLE_IDS < <(
    Rscript "${SCRIPT_DIR}/02_doublet_removal.Rscript" \
        --input-dir "${INPUT_DIR}" \
        --list-samples
)

if [ "${#SAMPLE_IDS[@]}" -eq 0 ]; then
    echo "No QC-filtered inputs were found."
    exit 1
fi

INDEX=$((SLURM_ARRAY_TASK_ID - 1))
if [ "${INDEX}" -lt 0 ] || [ "${INDEX}" -ge "${#SAMPLE_IDS[@]}" ]; then
    echo "Array index ${SLURM_ARRAY_TASK_ID} is out of range for ${#SAMPLE_IDS[@]} samples."
    exit 0
fi

SAMPLE_ID="${SAMPLE_IDS[${INDEX}]}"
echo "Running Stage 2 doublet removal for library ${SAMPLE_ID}"

Rscript "${SCRIPT_DIR}/02_doublet_removal.Rscript" \
    --input-dir "${INPUT_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --sample-ids "${SAMPLE_ID}" \
    --threads 4
