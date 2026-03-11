#!/bin/bash
#SBATCH -J dej_qc
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
# command-line --output/--error flags.  Setting them here as env vars
# does not work (SLURM reads log paths at submission time, not runtime).
mkdir -p "${DEJAGER_PROCESSING_LOGS}"

export HDF5_USE_FILE_LOCKING=FALSE

# Temporarily relax nounset for conda activation (some activate.d scripts
# reference unset variables like ADDR2LINE).
set +u
source "${CONDA_INIT_SCRIPT}"
conda activate "${QC_ENV}"
set -u
if [[ -z "${CONDA_PREFIX:-}" ]]; then
    echo "ERROR: Failed to activate conda environment: ${QC_ENV}"
    exit 1
fi

INPUT_DIR="${DEJAGER_PREPROCESSED}"
OUTPUT_DIR="${DEJAGER_QC_FILTERED}"

mapfile -t SAMPLE_IDS < <(
    python "${SCRIPT_DIR}/01_qc_filter.py" \
        --input-dir "${INPUT_DIR}" \
        --list-samples
)

if [ "${#SAMPLE_IDS[@]}" -eq 0 ]; then
    echo "No CellBender outputs were found."
    exit 1
fi

INDEX=$((SLURM_ARRAY_TASK_ID - 1))
if [ "${INDEX}" -lt 0 ] || [ "${INDEX}" -ge "${#SAMPLE_IDS[@]}" ]; then
    echo "Array index ${SLURM_ARRAY_TASK_ID} is out of range for ${#SAMPLE_IDS[@]} samples."
    exit 0
fi

SAMPLE_ID="${SAMPLE_IDS[${INDEX}]}"
echo "Running Stage 1 QC for library ${SAMPLE_ID}"

python "${SCRIPT_DIR}/01_qc_filter.py" \
    --input-dir "${INPUT_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --patient-map-csv "${DEJAGER_PATIENT_MAP}" \
    --patient-id-overrides-json "${DEJAGER_PATIENT_ID_OVERRIDES}" \
    --sample-ids "${SAMPLE_ID}"
