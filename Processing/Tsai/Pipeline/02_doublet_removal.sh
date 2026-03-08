#!/bin/bash
#SBATCH -J tsai_doublets
#SBATCH -t 12:00:00
#SBATCH -n 4
#SBATCH --mem=32G
#SBATCH --array=1-447%32
#SBATCH --mail-type=FAIL

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

source "${REPO_ROOT}/config/paths.sh"

# NOTE: SLURM log redirection is handled by submit_pipeline.sh via
# command-line --output/--error flags.
mkdir -p "${TSAI_PROCESSING_LOGS}"

source "${CONDA_INIT_SCRIPT}"
conda activate "${SINGLECELL_ENV}"
if [[ -z "${CONDA_PREFIX:-}" ]]; then
    echo "ERROR: Failed to activate conda environment: ${SINGLECELL_ENV}"
    exit 1
fi

INPUT_DIR="${TSAI_QC_FILTERED}"
OUTPUT_DIR="${TSAI_DOUBLET_REMOVED}"
METADATA_CSV="${TSAI_METADATA_CSV}"

mapfile -t SAMPLE_IDS < <(
    Rscript "${SCRIPT_DIR}/02_doublet_removal.Rscript" \
        --input-dir "${INPUT_DIR}" \
        --metadata-csv "${METADATA_CSV}" \
        --list-samples
)

if [ "${#SAMPLE_IDS[@]}" -eq 0 ]; then
    echo "No QC-filtered inputs were found."
    exit 1
fi

INDEX=$((SLURM_ARRAY_TASK_ID - 1))
if [ "${INDEX}" -lt 0 ] || [ "${INDEX}" -ge "${#SAMPLE_IDS[@]}" ]; then
    echo "Array index ${SLURM_ARRAY_TASK_ID} is out of range for ${#SAMPLE_IDS[@]} samples."
    exit 1
fi

SAMPLE_ID="${SAMPLE_IDS[${INDEX}]}"
echo "Running Stage 2 doublet removal for sample ${SAMPLE_ID}"

Rscript "${SCRIPT_DIR}/02_doublet_removal.Rscript" \
    --input-dir "${INPUT_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --metadata-csv "${METADATA_CSV}" \
    --sample-ids "${SAMPLE_ID}" \
    --threads 4
