#!/bin/bash
#SBATCH --job-name=copy_data
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --array=0-5

# Copy transcriptomics data into the canonical Data/Transcriptomics/ layout.
# Each array task handles one dataset. Rsync is used so re-running is safe
# (only copies new/changed files).
#
# NOTE: The original Openmind->repo copy (March 2026) is complete. This script
# is retained for future use if data needs to be re-synced from the primary
# storage at ROSMAP_Data/Single_Nucleus/ into the in-repo layout.
#
# Usage:
#   sbatch copy_data.sbatch
#
# Monitor:
#   squeue -u $USER -n copy_data
#   tail -f ~/copy_data_<jobid>_<taskid>.out

set -euo pipefail

# Source centralized config for all paths
_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${_SCRIPT_DIR}/../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

DEST_BASE="${REPO_ROOT}/Data/Transcriptomics"

# SLURM log location (set dynamically since we can't use variables in #SBATCH)
# Submit with: sbatch --output="${HOME}/copy_data_%A_%a.out" --error="${HOME}/copy_data_%A_%a.err" copy_data.sbatch

# Define source -> destination pairs (indexed by SLURM_ARRAY_TASK_ID)
# Sources use config variables pointing to ROSMAP_Data/Single_Nucleus/
declare -a SOURCES=(
    "${TSAI_FASTQS_DIR}/"
    "${TSAI_CELLRANGER}/"
    "${TSAI_CELLBENDER}/"
    "${DEJAGER_FASTQS}/"
    "${DEJAGER_CELLRANGER}/"
    "${DEJAGER_CELLBENDER}/"
)

declare -a DESTS=(
    "${DEST_BASE}/Tsai/FASTQs/"
    "${DEST_BASE}/Tsai/Cellranger_Output/"
    "${DEST_BASE}/Tsai/Cellbender_Output/"
    "${DEST_BASE}/DeJager/FASTQs/"
    "${DEST_BASE}/DeJager/Cellranger_Output/"
    "${DEST_BASE}/DeJager/Cellbender_Output/"
)

declare -a LABELS=(
    "Tsai FASTQs"
    "Tsai Cellranger"
    "Tsai CellBender"
    "DeJager FASTQs"
    "DeJager Cellranger"
    "DeJager CellBender"
)

IDX="${SLURM_ARRAY_TASK_ID}"
SRC="${SOURCES[$IDX]}"
DST="${DESTS[$IDX]}"
LABEL="${LABELS[$IDX]}"

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

log "=========================================="
log "Dataset: ${LABEL} (task ${IDX})"
log "Source:  ${SRC}"
log "Dest:   ${DST}"
log "=========================================="

# Verify source exists
if [[ ! -d "${SRC}" ]]; then
    log "ERROR: Source directory does not exist: ${SRC}"
    exit 1
fi

# Create destination if needed
mkdir -p "${DST}"

# Count source items
SRC_COUNT=$(ls -1 "${SRC}" | wc -l)
log "Source contains ${SRC_COUNT} items"

# Rsync with progress. Using -a (archive) preserves permissions, timestamps, etc.
# --info=progress2 gives a compact overall progress indicator.
log "Starting rsync..."
rsync -a \
    --info=progress2 \
    --human-readable \
    "${SRC}" "${DST}" \
    2>&1

RC=$?
if [[ ${RC} -eq 0 ]]; then
    DST_COUNT=$(ls -1 "${DST}" | wc -l)
    log "DONE: ${LABEL} - ${DST_COUNT} items in destination"
else
    log "ERROR: rsync exited with code ${RC}"
    exit ${RC}
fi
