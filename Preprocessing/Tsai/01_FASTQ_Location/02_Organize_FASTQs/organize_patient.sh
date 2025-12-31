#!/bin/bash
#
# Worker script: Organize FASTQ files for a single patient (projid)
# Creates symbolic links in the output directory structure
#
# Usage: ./organize_patient.sh <projid> <fastq_csv> <output_dir>
#

set -e

# Arguments
PROJID="$1"
FASTQ_CSV="$2"
OUTPUT_BASE_DIR="$3"

if [[ -z "$PROJID" ]] || [[ -z "$FASTQ_CSV" ]] || [[ -z "$OUTPUT_BASE_DIR" ]]; then
    echo "Usage: $0 <projid> <fastq_csv> <output_dir>" >&2
    exit 1
fi

# Create patient directory
PATIENT_DIR="${OUTPUT_BASE_DIR}/${PROJID}"
mkdir -p "$PATIENT_DIR"

# Track Library_ID occurrences to create _1, _2, etc. subdirectories
declare -A LIB_COUNTS

# Read the FASTQ CSV and process entries for this projid
# Expected columns: projid,Library_ID,source_dir,fastq_filename,full_path,file_size
# Skip header and filter for this projid

grep "^${PROJID}," "$FASTQ_CSV" | while IFS=',' read -r proj lib_id source_dir fastq_filename full_path file_size; do
    # Skip if empty
    if [[ -z "$lib_id" ]] || [[ -z "$full_path" ]]; then
        continue
    fi
    
    # Determine the occurrence number for this Library_ID at this source_dir
    # The key is library_id + source_dir to handle multiple sequencing runs
    KEY="${lib_id}|${source_dir}"
    
    # Read current count from a tracking file (since we're in a subshell)
    COUNT_FILE="${PATIENT_DIR}/.lib_counts"
    touch "$COUNT_FILE"
    
    # Get or initialize count for this key
    CURRENT_COUNT=$(grep "^${KEY}=" "$COUNT_FILE" 2>/dev/null | cut -d'=' -f2 || echo "0")
    if [[ -z "$CURRENT_COUNT" ]] || [[ "$CURRENT_COUNT" == "0" ]]; then
        CURRENT_COUNT=1
        echo "${KEY}=${CURRENT_COUNT}" >> "$COUNT_FILE"
    fi
    
    # Create subdirectory for this Library_ID occurrence
    LIB_SUBDIR="${PATIENT_DIR}/${lib_id}_${CURRENT_COUNT}"
    mkdir -p "$LIB_SUBDIR"
    
    # Create symbolic link if source file exists and link doesn't already exist
    LINK_PATH="${LIB_SUBDIR}/${fastq_filename}"
    if [[ -f "$full_path" ]] && [[ ! -e "$LINK_PATH" ]]; then
        ln -s "$full_path" "$LINK_PATH"
    fi
done

# Cleanup tracking file
rm -f "${PATIENT_DIR}/.lib_counts"

# Output success message
NUM_LINKS=$(find "$PATIENT_DIR" -type l 2>/dev/null | wc -l)
echo "Organized ${NUM_LINKS} FASTQ files for projid ${PROJID}"

