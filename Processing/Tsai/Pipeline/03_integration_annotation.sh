#!/bin/bash
#SBATCH -J tsai_integrate
#SBATCH -t 48:00:00
#SBATCH -n 32
#SBATCH --mem=256G
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

source "${REPO_ROOT}/config/paths.sh"

mkdir -p "${TSAI_PROCESSING_LOGS}"
export SLURM_OUTPUT="${TSAI_PROCESSING_LOGS}/tsai_integrate_%j.out"
export SLURM_ERROR="${TSAI_PROCESSING_LOGS}/tsai_integrate_%j.err"

source "${CONDA_INIT_SCRIPT}"
conda activate "${BATCHCORR_ENV}"

INPUT_DIR="${TSAI_DOUBLET_REMOVED}"
OUTPUT_DIR="${TSAI_INTEGRATED}"
METADATA_CSV="${TSAI_METADATA_CSV}"
MARKERS_RDS="${TSAI_MARKERS_RDS}"

python "${SCRIPT_DIR}/03_integration_annotation.py" \
    --input-dir "${INPUT_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --metadata-csv "${METADATA_CSV}" \
    --markers-rds "${MARKERS_RDS}"
