#!/bin/bash
#SBATCH -J tsai_qc
#SBATCH -t 12:00:00
#SBATCH -n 4
#SBATCH --mem=32G
#SBATCH --array=1-447%32
#SBATCH --mail-type=FAIL

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

source "${REPO_ROOT}/config/paths.sh"

# Redirect SLURM logs (create directory if needed)
mkdir -p "${TSAI_PROCESSING_LOGS}"
export SLURM_OUTPUT="${TSAI_PROCESSING_LOGS}/tsai_qc_%A_%a.out"
export SLURM_ERROR="${TSAI_PROCESSING_LOGS}/tsai_qc_%A_%a.err"

source "${CONDA_INIT_SCRIPT}"
conda activate "${QC_ENV}"

INPUT_DIR="${WORKSPACE_ROOT}/Tsai_Data/Cellbender_Outputs"
OUTPUT_DIR="${TSAI_QC_FILTERED}"
METADATA_CSV="${TSAI_METADATA_CSV}"

mapfile -t SAMPLE_IDS < <(
    python "${SCRIPT_DIR}/01_qc_filter.py" \
        --input-dir "${INPUT_DIR}" \
        --metadata-csv "${METADATA_CSV}" \
        --list-samples
)

if [ "${#SAMPLE_IDS[@]}" -eq 0 ]; then
    echo "No complete CellBender outputs were found."
    exit 1
fi

INDEX=$((SLURM_ARRAY_TASK_ID - 1))
if [ "${INDEX}" -lt 0 ] || [ "${INDEX}" -ge "${#SAMPLE_IDS[@]}" ]; then
    echo "Array index ${SLURM_ARRAY_TASK_ID} is out of range for ${#SAMPLE_IDS[@]} samples."
    exit 1
fi

SAMPLE_ID="${SAMPLE_IDS[${INDEX}]}"
echo "Running Stage 1 QC for sample ${SAMPLE_ID}"

python "${SCRIPT_DIR}/01_qc_filter.py" \
    --input-dir "${INPUT_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --metadata-csv "${METADATA_CSV}" \
    --sample-ids "${SAMPLE_ID}"
