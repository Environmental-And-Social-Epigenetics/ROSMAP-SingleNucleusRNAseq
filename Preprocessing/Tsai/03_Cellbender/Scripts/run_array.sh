#!/bin/bash
#SBATCH --job-name=cellbender_array
#SBATCH --array=1-243%10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=47:00:00
#SBATCH --gres=gpu:1
#SBATCH --mem=128G
#SBATCH --output=/om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/03_Cellbender/Logs/Outs/array_%A_%a.out
#SBATCH --error=/om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/03_Cellbender/Logs/Errs/array_%A_%a.err
#SBATCH --mail-user=mabdel03@mit.edu
#SBATCH --mail-type=FAIL

# CellBender array job - max 10 concurrent

# Disable HDF5 file locking to avoid shared filesystem locking issues
export HDF5_USE_FILE_LOCKING=FALSE

TRACKING_DIR="/om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/03_Cellbender/Tracking"
NEEDS_RERUN="${TRACKING_DIR}/needs_rerun.txt"

# Get projid for this array task
PROJID=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$NEEDS_RERUN")

if [[ -z "$PROJID" ]]; then
    echo "ERROR: No projid found for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "Processing projid: $PROJID (task ${SLURM_ARRAY_TASK_ID})"

# Paths
INPUT_H5="/om/scratch/Mon/mabdel03/Tsai_Data/Cellranger_Outputs/${PROJID}/outs/raw_feature_bc_matrix.h5"
OUTPUT_DIR="/om/scratch/Mon/mabdel03/Tsai_Data/Cellbender_Outputs/${PROJID}"
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
source /om2/user/mabdel03/anaconda/etc/profile.d/conda.sh
conda activate /om2/user/mabdel03/conda_envs/Cellbender_env

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
