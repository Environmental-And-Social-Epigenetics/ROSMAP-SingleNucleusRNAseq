#!/bin/bash
#
# CellBender array job - max 10 concurrent
#
# Source config/paths.sh before submitting, or ensure the environment variables
# are set.  SBATCH directives below use values from the sourced config.
#

# Source config/paths.sh relative to this script's location (depth 4 from repo root)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

PIPELINE_DIR="${REPO_ROOT}/Preprocessing/Tsai/03_Cellbender"

#SBATCH --job-name=cellbender_array
#SBATCH --array=1-36%10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=47:00:00
#SBATCH --gres=gpu:1
#SBATCH --mem=128G
#SBATCH --output=${PIPELINE_DIR}/Logs/Outs/array_%A_%a.out
#SBATCH --error=${PIPELINE_DIR}/Logs/Errs/array_%A_%a.err
#SBATCH --mail-type=FAIL

# CellBender array job - max 10 concurrent

# Disable HDF5 file locking to avoid shared filesystem locking issues
export HDF5_USE_FILE_LOCKING=FALSE

TRACKING_DIR="${PIPELINE_DIR}/Tracking"
NEEDS_RERUN="${TRACKING_DIR}/needs_rerun.txt"

# Get projid for this array task
PROJID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$NEEDS_RERUN")

if [[ -z "$PROJID" ]]; then
    echo "ERROR: No projid found for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "Processing projid: $PROJID (task ${SLURM_ARRAY_TASK_ID})"

# Paths — use variables from config/paths.sh
INPUT_H5="${TSAI_CELLRANGER_OUTPUT}/${PROJID}/outs/raw_feature_bc_matrix.h5"
OUTPUT_DIR="${TSAI_PREPROCESSED}/${PROJID}"
OUTPUT_H5="${OUTPUT_DIR}/processed_feature_bc_matrix.h5"

# Check if already completed
if [[ -f "$OUTPUT_H5" ]]; then
    echo "Output already exists, skipping: $OUTPUT_H5"
    exit 0
fi

# Check input exists
if [[ ! -f "$INPUT_H5" ]]; then
    echo "ERROR: Input not found: $INPUT_H5"
    exit 1
fi

# Initialize conda
source "${CONDA_INIT_SCRIPT}"
conda activate "${CELLBENDER_ENV}"

# Create output directory
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

# Run CellBender
cellbender remove-background --cuda --input "$INPUT_H5" --fpr 0 --output "$OUTPUT_H5"

# Mark completion if successful
if [[ -f "$OUTPUT_H5" ]]; then
    echo "$PROJID" >> "${TRACKING_DIR}/cellbender_completed.txt"
    echo "CellBender completed for $PROJID"
else
    echo "$PROJID" >> "${TRACKING_DIR}/cellbender_failed.txt"
    echo "CellBender failed for $PROJID"
    exit 1
fi
