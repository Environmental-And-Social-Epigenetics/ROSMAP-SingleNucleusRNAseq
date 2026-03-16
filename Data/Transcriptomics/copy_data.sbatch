#!/bin/bash
#SBATCH --job-name=copy_data
#SBATCH --partition=normal
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --array=0-5
#SBATCH --output=/home/mabdel03/copy_data_%A_%a.out
#SBATCH --error=/home/mabdel03/copy_data_%A_%a.err

# Copy transcriptomics data into the canonical Data/Transcriptomics/ layout.
# Each array task handles one dataset. Rsync is used so re-running is safe
# (only copies new/changed files).
#
# Usage:
#   sbatch copy_data.sbatch
#
# Monitor:
#   squeue -u $USER -n copy_data
#   tail -f ~/copy_data_<jobid>_<taskid>.out

set -euo pipefail

REPO_ROOT="/om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq"
DEST_BASE="${REPO_ROOT}/Data/Transcriptomics"

# Define source -> destination pairs (indexed by SLURM_ARRAY_TASK_ID)
declare -a SOURCES=(
    "/om/scratch/Mon/mabdel03/Tsai_Data/FASTQs/"
    "/om/scratch/Mon/mabdel03/Tsai_Data/Cellranger_Outputs/"
    "/om/scratch/Mon/mabdel03/Tsai_Data/Cellbender_Outputs/"
    "/om/scratch/Mon/mabdel03/FASTQs/"
    "/om/scratch/Mon/mabdel03/Counts/"
    "/om2/user/mabdel03/files/ACE_Analysis/Data/DeJager/Preprocessed_Counts/"
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
