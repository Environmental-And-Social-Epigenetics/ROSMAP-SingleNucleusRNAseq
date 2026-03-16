#!/bin/bash
# =============================================================================
# SLURM array job — Upload DeJager CellRanger outputs to NAS (47 parallel workers)
# =============================================================================
#SBATCH --job-name=nas_cellranger_dejager
#SBATCH --array=0-46
#SBATCH --time=2-00:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=/home/%u/nas_cellranger_dejager_%A_%a.out
#SBATCH --error=/home/%u/nas_cellranger_dejager_%A_%a.err

set -euo pipefail

# Use SLURM_SUBMIT_DIR to find the original script location
SCRIPT_DIR="${SLURM_SUBMIT_DIR}/DeJager"

# CellRanger only
export DATA_TYPES="Cellranger_Output"

# 47 parallel workers (one per CellRanger folder)
export SFTP_SHARD_TOTAL=47
export SFTP_SHARD_INDEX="${SLURM_ARRAY_TASK_ID}"

echo "=== SLURM Array Task ${SLURM_ARRAY_TASK_ID} / 47 ==="
echo "Data types: ${DATA_TYPES}"
echo "Shard: ${SFTP_SHARD_INDEX} / ${SFTP_SHARD_TOTAL}"
echo "Host: $(hostname)"
echo "Start: $(date)"
echo ""

"${SCRIPT_DIR}/upload_to_nas.sh" upload

echo ""
echo "Finished: $(date)"
