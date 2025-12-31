#!/bin/bash
#
# Worker script: List all FASTQ files for a given path and Library_ID
# This mirrors the list_remote_files logic from the original notebooks
#
# Usage: ./list_fastqs_for_path.sh <row_index> <projid> <library_id> <tsai_path> <backup_path> <output_file>
#

set -e

# Arguments
ROW_INDEX="$1"
PROJID="$2"
LIBRARY_ID="$3"
TSAI_PATH="$4"
BACKUP_PATH="$5"
OUTPUT_FILE="$6"

# Function to list FASTQ files matching Library_ID pattern
list_fastqs() {
    local dir_path="$1"
    local lib_id="$2"
    
    if [[ -d "$dir_path" ]]; then
        # Find all fastq.gz files matching the Library_ID
        find "$dir_path" -maxdepth 1 -type f -name "${lib_id}*.fastq.gz" 2>/dev/null || true
    fi
}

# Function to get file size
get_file_size() {
    local file_path="$1"
    if [[ -f "$file_path" ]]; then
        stat -c %s "$file_path" 2>/dev/null || echo "0"
    else
        echo "0"
    fi
}

# Try primary path first
FASTQ_FILES=$(list_fastqs "$TSAI_PATH" "$LIBRARY_ID")

# If primary is empty/missing, try backup
if [[ -z "$FASTQ_FILES" ]] && [[ -n "$BACKUP_PATH" ]] && [[ "$BACKUP_PATH" != "nan" ]] && [[ "$BACKUP_PATH" != "NaN" ]]; then
    FASTQ_FILES=$(list_fastqs "$BACKUP_PATH" "$LIBRARY_ID")
    if [[ -n "$FASTQ_FILES" ]]; then
        # Update source path to backup
        TSAI_PATH="$BACKUP_PATH"
    fi
fi

# Output each file as a CSV row
if [[ -n "$FASTQ_FILES" ]]; then
    while IFS= read -r fastq_file; do
        if [[ -n "$fastq_file" ]]; then
            filename=$(basename "$fastq_file")
            filesize=$(get_file_size "$fastq_file")
            # Output: projid,Library_ID,source_dir,fastq_filename,full_path,file_size
            echo "${PROJID},${LIBRARY_ID},${TSAI_PATH},${filename},${fastq_file},${filesize}" >> "$OUTPUT_FILE"
        fi
    done <<< "$FASTQ_FILES"
else
    # Log missing files
    echo "WARNING: No FASTQ files found for ${LIBRARY_ID} in ${TSAI_PATH} or backup" >&2
fi

