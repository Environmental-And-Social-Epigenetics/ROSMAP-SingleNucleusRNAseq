#!/bin/bash
#SBATCH --job-name=organize_fastqs
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --output=/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/Tsai/Preprocessing/FASTQ_Transfer/New/Logs/organize_%j.out
#SBATCH --error=/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/Tsai/Preprocessing/FASTQ_Transfer/New/Logs/organize_%j.err
#SBATCH --mail-user=mabdel03@mit.edu
#SBATCH --mail-type=BEGIN,END,FAIL

#
# SLURM batch script to organize FASTQs by patient (projid)
# Creates symbolic links in FASTQs_By_Patient/{projid}/{Library_ID}_N/
#

set -e

# Define absolute paths (required for SLURM which copies the script)
PIPELINE_ROOT="/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/Tsai/Preprocessing/FASTQ_Transfer"
SCRIPT_DIR="${PIPELINE_ROOT}/New/Scripts/02_Organize_FASTQs"

# Load configuration
source "${PIPELINE_ROOT}/New/Scripts/Config/config.sh"

# Activate conda environment
source "$CONDA_INIT_SCRIPT"
conda activate "$CONDA_ENV"

# Ensure output directories exist
ensure_dirs

echo "=============================================="
echo "Organizing FASTQs by Patient"
echo "Started: $(date)"
echo "=============================================="

# Check that master FASTQ CSV exists
if [[ ! -f "$ALL_FASTQS_CSV" ]]; then
    echo "ERROR: Master FASTQ CSV not found: ${ALL_FASTQS_CSV}"
    echo "Please run run_build_csv.sh first."
    exit 1
fi

# Create a preprocessed CSV that groups by projid and Library_ID occurrence
# This handles the _1, _2 numbering per the original notebook convention
echo "Preprocessing FASTQ CSV to assign occurrence numbers..."

PREPROC_CSV="${TEMP_DIR}/fastqs_with_occurrence.csv"
mkdir -p "$TEMP_DIR"

python3 << 'PYTHON_SCRIPT' - "$ALL_FASTQS_CSV" "$PREPROC_CSV"
import sys
import pandas as pd
from collections import defaultdict

input_csv = sys.argv[1]
output_csv = sys.argv[2]

df = pd.read_csv(input_csv)

# Group by projid and Library_ID, then assign occurrence numbers based on source_dir
# This mimics the original notebook logic where each unique source_dir for a Library_ID gets _1, _2, etc.

occurrence_map = defaultdict(lambda: defaultdict(int))
occurrences = []

for idx, row in df.iterrows():
    projid = str(row['projid'])
    lib_id = row['Library_ID']
    source_dir = row['source_dir']
    
    key = f"{projid}|{lib_id}"
    
    # Get or assign occurrence number for this source_dir
    if source_dir not in occurrence_map[key]:
        occurrence_map[key][source_dir] = len(occurrence_map[key]) + 1
    
    occurrences.append(occurrence_map[key][source_dir])

df['occurrence'] = occurrences
df.to_csv(output_csv, index=False)

print(f"Preprocessed CSV saved with {len(df)} records")
print(f"Unique projids: {df['projid'].nunique()}")
PYTHON_SCRIPT

# Get unique projids
PROJIDS=$(tail -n +2 "$PREPROC_CSV" | cut -d',' -f1 | sort -u)
NUM_PROJIDS=$(echo "$PROJIDS" | wc -l)

echo "Found ${NUM_PROJIDS} unique projids to process"
echo "Output directory: ${FASTQ_OUTPUT_DIR}"

# Create jobs file
JOBS_FILE="${TEMP_DIR}/organize_jobs.txt"
echo "$PROJIDS" > "$JOBS_FILE"

# Worker script path
WORKER_SCRIPT="${SCRIPT_DIR}/organize_patient.sh"
chmod +x "$WORKER_SCRIPT"

# Run organization in parallel using xargs (GNU parallel not available)
echo "Running parallel organization with ${NUM_PARALLEL_JOBS} workers..."

cat "$JOBS_FILE" | xargs -P "$NUM_PARALLEL_JOBS" -I {} \
    bash "$WORKER_SCRIPT" {} "$PREPROC_CSV" "$FASTQ_OUTPUT_DIR"

echo "Parallel organization completed."

# Run validation
echo "Running validation..."
python3 "${SCRIPT_DIR}/validate_organization.py" \
    --fastq-csv "$ALL_FASTQS_CSV" \
    --output-dir "$FASTQ_OUTPUT_DIR" \
    --report "$VALIDATION_REPORT"

# Cleanup
rm -f "$PREPROC_CSV" "$JOBS_FILE"

echo "=============================================="
echo "FASTQ Organization Complete!"
echo "Output directory: ${FASTQ_OUTPUT_DIR}"
echo "Validation report: ${VALIDATION_REPORT}"
echo "Completed: $(date)"
echo "=============================================="

