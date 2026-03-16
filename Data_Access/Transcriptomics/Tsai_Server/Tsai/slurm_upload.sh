#!/bin/bash
# =============================================================================
# slurm_upload.sh — SLURM array job wrapper for Tsai NAS upload
# =============================================================================
# Submit with:
#   sbatch slurm_upload.sh
#
# Launches 20 parallel workers, each responsible for a shard of the folder
# list.  Folder-to-worker assignment is deterministic (hash-based), so the
# job is safely re-submittable: completed folders are skipped via state files.
# =============================================================================

#SBATCH --job-name=nas_upload_tsai
#SBATCH --array=0-19
#SBATCH --time=2-00:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --output=/home/%u/nas_upload_tsai_%A_%a.out
#SBATCH --error=/home/%u/nas_upload_tsai_%A_%a.err

set -euo pipefail

# Use SLURM_SUBMIT_DIR to find the original script location
SCRIPT_DIR="${SLURM_SUBMIT_DIR}/Tsai"

# ---------------------------------------------------------------------------
# Set sharding environment variables for the upload script
# ---------------------------------------------------------------------------
export SFTP_SHARD_TOTAL=20
export SFTP_SHARD_INDEX="${SLURM_ARRAY_TASK_ID}"

echo "=== SLURM Array Task ${SLURM_ARRAY_TASK_ID} of ${SLURM_ARRAY_TASK_COUNT:-20} ==="
echo "Shard: ${SFTP_SHARD_INDEX} / ${SFTP_SHARD_TOTAL}"
echo "Host: $(hostname)"
echo "Start: $(date)"
echo ""

# ---------------------------------------------------------------------------
# Run the upload script in upload mode
# ---------------------------------------------------------------------------
"${SCRIPT_DIR}/upload_to_nas.sh" upload

echo ""
echo "Finished: $(date)"
