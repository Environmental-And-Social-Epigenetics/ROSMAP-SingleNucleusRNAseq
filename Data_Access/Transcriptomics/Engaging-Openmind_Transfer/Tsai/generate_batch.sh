#!/bin/bash
#
# Generate Globus batch file for Tsai transcriptomics data
#
# Scans the local Tsai data directories (FASTQs, Cellranger_Output,
# Cellbender_Output) and generates a Globus batch file with one line per
# subfolder.  The batch file can then be consumed by send_globus.sh or
# receive_globus.sh.
#
# Usage:
#   ./generate_batch.sh [DATA_TYPE]
#
# Arguments:
#   DATA_TYPE   One of: FASTQs, Cellranger_Output, Cellbender_Output, all
#               (default: all)
#
# Output:
#   Batch_Files/tsai_<DATA_TYPE>_batch.txt
#

set -euo pipefail

# =============================================================================
# Source project paths
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

# =============================================================================
# Configuration
# =============================================================================

COHORT="Tsai"
DATA_TYPES=("FASTQs" "Cellranger_Output" "Cellbender_Output")
LOCAL_BASE="${REPO_ROOT}/Data/Transcriptomics/${COHORT}"
BATCH_DIR="${SCRIPT_DIR}/Batch_Files"

# =============================================================================
# Functions
# =============================================================================

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2
}

usage() {
    echo "Usage: $0 [DATA_TYPE]"
    echo ""
    echo "DATA_TYPE: FASTQs | Cellranger_Output | Cellbender_Output | all (default: all)"
}

generate_batch_for_type() {
    local data_type="$1"
    local batch_file="$2"
    local source_dir="${LOCAL_BASE}/${data_type}"
    local count=0

    if [[ ! -d "$source_dir" ]]; then
        error "Directory not found: ${source_dir}"
        return 1
    fi

    # Scan for subdirectories (each is a sample/library)
    for entry in "${source_dir}"/*/; do
        # Skip if no matches (glob didn't expand)
        [[ -d "$entry" ]] || continue

        local folder_name
        folder_name="$(basename "$entry")"

        # Write batch line: source_path dest_path --recursive
        # Paths are relative to the endpoint root; we use the full absolute
        # path on each cluster since both endpoints expose the filesystem root.
        local rel_path="Data/Transcriptomics/${COHORT}/${data_type}/${folder_name}"
        echo "${REPO_ROOT}/${rel_path} ${REPO_ROOT}/${rel_path} --recursive" >> "$batch_file"
        count=$((count + 1))
    done

    echo "$count"
}

# =============================================================================
# Parse arguments
# =============================================================================

FILTER="${1:-all}"

if [[ "$FILTER" != "all" ]]; then
    # Validate the filter
    valid=false
    for dt in "${DATA_TYPES[@]}"; do
        if [[ "$dt" == "$FILTER" ]]; then
            valid=true
            break
        fi
    done
    if [[ "$valid" == false ]]; then
        error "Invalid data type: $FILTER"
        usage
        exit 1
    fi
    DATA_TYPES=("$FILTER")
fi

# =============================================================================
# Main
# =============================================================================

log "Generating Globus batch file for ${COHORT} data"
log "Data types: ${DATA_TYPES[*]}"
log "Local base: ${LOCAL_BASE}"
echo

mkdir -p "$BATCH_DIR"

TIMESTAMP="$(date +%Y%m%d_%H%M%S)"
if [[ "$FILTER" == "all" ]]; then
    BATCH_FILE="${BATCH_DIR}/tsai_all_batch_${TIMESTAMP}.txt"
else
    BATCH_FILE="${BATCH_DIR}/tsai_${FILTER}_batch_${TIMESTAMP}.txt"
fi

# Start with an empty file
> "$BATCH_FILE"

TOTAL=0
for data_type in "${DATA_TYPES[@]}"; do
    log "Scanning ${data_type}..."
    n=$(generate_batch_for_type "$data_type" "$BATCH_FILE")
    log "  Found ${n} entries in ${data_type}"
    TOTAL=$((TOTAL + n))
done

echo
log "Batch file written: ${BATCH_FILE}"
log "Total entries: ${TOTAL}"

if [[ "$TOTAL" -eq 0 ]]; then
    log "WARNING: No entries generated. Check that data directories exist and contain subdirectories."
fi

# Also create/update a symlink to the latest batch file for convenience
LATEST_LINK="${BATCH_DIR}/tsai_latest_batch.txt"
ln -sf "$(basename "$BATCH_FILE")" "$LATEST_LINK"
log "Symlink updated: ${LATEST_LINK} -> $(basename "$BATCH_FILE")"
