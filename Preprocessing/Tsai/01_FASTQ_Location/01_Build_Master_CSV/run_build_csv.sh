#!/bin/bash
#SBATCH --job-name=build_fastq_csv
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --output=/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/Tsai/Preprocessing/FASTQ_Transfer/New/Logs/build_csv_%j.out
#SBATCH --error=/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/Tsai/Preprocessing/FASTQ_Transfer/New/Logs/build_csv_%j.err
#SBATCH --mail-user=mabdel03@mit.edu
#SBATCH --mail-type=BEGIN,END,FAIL

#
# SLURM batch script to build master FASTQ CSV
# Uses GNU parallel to process all paths from Tsai_To_Openmind.csv
#

set -e

# Define absolute paths (required for SLURM which copies the script)
PIPELINE_ROOT="/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/Tsai/Preprocessing/FASTQ_Transfer"
SCRIPT_DIR="${PIPELINE_ROOT}/New/Scripts/01_Build_Master_CSV"

# Load configuration
source "${PIPELINE_ROOT}/New/Scripts/Config/config.sh"

# Activate conda environment
source "$CONDA_INIT_SCRIPT"
conda activate "$CONDA_ENV"

# Ensure output directories exist
ensure_dirs

echo "=============================================="
echo "Building Master FASTQ CSV"
echo "Started: $(date)"
echo "=============================================="

# Create temporary directory for this run
RUN_TEMP="${TEMP_DIR}/build_csv_$$"
mkdir -p "$RUN_TEMP"

# Prepare jobs file from master CSV
# Skip header, extract: row_index, projid, Library_ID, Tsai_path
JOBS_FILE="${RUN_TEMP}/jobs.txt"
echo "Preparing jobs from ${MASTER_CSV}..."

# Load backups into an associative array for quick lookup
declare -A BACKUPS_MAP
if [[ -f "$BACKUPS_CSV" ]]; then
    echo "Loading backup paths from ${BACKUPS_CSV}..."
    # Read backups CSV and build lookup (Origin_1 -> Backup_1, Origin_2 -> Backup_2)
    while IFS=',' read -r idx projid libid orig status numpaths backup1 null backup2 origin1 origin2; do
        if [[ -n "$origin1" ]] && [[ "$origin1" != "Origin_1" ]]; then
            BACKUPS_MAP["$origin1"]="$backup1"
        fi
        if [[ -n "$origin2" ]] && [[ "$origin2" != "Origin_2" ]]; then
            BACKUPS_MAP["$origin2"]="$backup2"
        fi
    done < "$BACKUPS_CSV"
fi

# Process master CSV and create jobs file
ROW_NUM=0
tail -n +2 "$MASTER_CSV" | while IFS=',' read -r idx idx2 batch projid library_id tsai_path openmind_dest; do
    ROW_NUM=$((ROW_NUM + 1))
    
    # Clean up fields (remove quotes if present)
    projid=$(echo "$projid" | tr -d '"')
    library_id=$(echo "$library_id" | tr -d '"')
    tsai_path=$(echo "$tsai_path" | tr -d '"')
    
    # Get backup path if available
    backup_path="${BACKUPS_MAP[$tsai_path]:-}"
    
    # Write job: row_num|projid|library_id|tsai_path|backup_path
    echo "${ROW_NUM}|${projid}|${library_id}|${tsai_path}|${backup_path}"
done > "$JOBS_FILE"

NUM_JOBS=$(wc -l < "$JOBS_FILE")
echo "Created ${NUM_JOBS} jobs for processing"

# Create output file with header
OUTPUT_CSV="${RUN_TEMP}/fastqs_combined.csv"
echo "projid,Library_ID,source_dir,fastq_filename,full_path,file_size" > "$OUTPUT_CSV"

# Worker output directory
WORKER_OUTPUT_DIR="${RUN_TEMP}/worker_outputs"
mkdir -p "$WORKER_OUTPUT_DIR"

# Run jobs in parallel using xargs (GNU parallel not available)
echo "Running parallel jobs with ${NUM_PARALLEL_JOBS} workers..."
WORKER_SCRIPT="${SCRIPT_DIR}/list_fastqs_for_path.sh"

# Use xargs for parallelization
cat "$JOBS_FILE" | xargs -P "$NUM_PARALLEL_JOBS" -I {} bash -c '
    IFS="|" read -r row_num projid library_id tsai_path backup_path <<< "{}"
    bash "'"$WORKER_SCRIPT"'" "$row_num" "$projid" "$library_id" "$tsai_path" "$backup_path" "'"$WORKER_OUTPUT_DIR"'/worker_${row_num}.csv"
'

echo "Parallel jobs completed. Merging results..."

# Merge all worker outputs
for worker_file in "$WORKER_OUTPUT_DIR"/worker_*.csv; do
    if [[ -f "$worker_file" ]]; then
        cat "$worker_file" >> "$OUTPUT_CSV"
    fi
done

# Run the Python script to finalize and validate
echo "Finalizing CSV with Python script..."
python3 "${SCRIPT_DIR}/merge_and_finalize_csv.py" \
    --input "$OUTPUT_CSV" \
    --output "$ALL_FASTQS_CSV" \
    --master-csv "$MASTER_CSV"

# Cleanup temp files
rm -rf "$RUN_TEMP"

echo "=============================================="
echo "Master FASTQ CSV built successfully!"
echo "Output: ${ALL_FASTQS_CSV}"
echo "Completed: $(date)"
echo "=============================================="

