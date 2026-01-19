#!/bin/bash
#
# Cell Ranger + CellBender Pipeline Configuration
# For processing all 480 Tsai patients in batches
#

# =============================================================================
# LOAD GLOBAL PATHS
# =============================================================================

# Get the directory where this script is located
CELLRANGER_CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source the main paths configuration
# Path: Config/ -> 02_Cellranger_Counts/ -> Tsai/ -> Preprocessing/ -> ROSMAP-SingleNucleusRNAseq/
source "${CELLRANGER_CONFIG_DIR}/../../../../config/paths.sh"

# =============================================================================
# PIPELINE DIRECTORIES
# =============================================================================

# Pipeline root
export PIPELINE_DIR="${REPO_ROOT}/Preprocessing/Tsai/02_Cellranger_Counts"

# Scripts directory
export SCRIPTS_DIR="${PIPELINE_DIR}/Scripts"

# Batch scripts output
export BATCH_SCRIPTS_DIR="${PIPELINE_DIR}/Batch_Scripts"

# Logs
export LOGS_DIR="${PIPELINE_DIR}/Logs"
export LOGS_OUT="${LOGS_DIR}/Outs"
export LOGS_ERR="${LOGS_DIR}/Errs"

# Tracking
export TRACKING_DIR="${PIPELINE_DIR}/Tracking"

# =============================================================================
# BATCH CONFIGURATION
# =============================================================================

# Number of patients per batch (conservative: 80 x 10GB = 800GB < 1TB)
export BATCH_SIZE=80

# Total number of batches (480 patients / 80 = 6 batches)
export NUM_BATCHES=6

# =============================================================================
# CELL RANGER PARAMETERS
# =============================================================================

# Cell Ranger version/path
export CELLRANGER_BIN="${CELLRANGER_PATH}/cellranger"

# Reference transcriptome
export TRANSCRIPTOME="${CELLRANGER_REF}"

# Cell Ranger options
export CR_CREATE_BAM="false"           # No BAM needed (patient assignment known)
export CR_INCLUDE_INTRONS="true"       # Include intronic reads (nuclear RNA)
export CR_NOSECONDARY="true"           # Skip secondary analysis
export CR_R1_LENGTH="26"               # Trim R1 to 26bp

# =============================================================================
# CELLBENDER PARAMETERS
# =============================================================================

# CellBender environment
export CELLBENDER_CONDA="${CELLBENDER_ENV}"

# CellBender options
export CB_EXPECTED_CELLS="5000"        # Expected number of cells
export CB_TOTAL_DROPLETS="20000"       # Total droplets to include
export CB_FPR="0.01"                   # False positive rate
export CB_EPOCHS="150"                 # Training epochs
export CB_CUDA="true"                  # Use GPU if available

# =============================================================================
# SLURM PARAMETERS
# =============================================================================

# Cell Ranger SLURM settings
export CR_SLURM_TIME="47:00:00"
export CR_SLURM_CPUS="32"
export CR_SLURM_MEM="128G"

# CellBender SLURM settings (GPU job)
export CB_SLURM_TIME="4:00:00"
export CB_SLURM_CPUS="4"
export CB_SLURM_MEM="64G"
export CB_SLURM_GPU="1"

# Email for notifications
export SLURM_MAIL_USER="mabdel03@mit.edu"

# =============================================================================
# INPUT/OUTPUT PATHS
# =============================================================================

# Input CSV with all FASTQ locations
export INPUT_CSV="${TSAI_FASTQS_CSV}"

# Scratch directories (temporary)
export CELLRANGER_OUTPUT="${TSAI_CELLRANGER_SCRATCH}"
export CELLBENDER_OUTPUT="${TSAI_CELLBENDER_SCRATCH}"

# Permanent storage for final CellBender outputs
export FINAL_OUTPUT="${TSAI_PREPROCESSED}"

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Ensure all directories exist
ensure_pipeline_dirs() {
    mkdir -p "${BATCH_SCRIPTS_DIR}"
    mkdir -p "${LOGS_OUT}"
    mkdir -p "${LOGS_ERR}"
    mkdir -p "${TRACKING_DIR}"
    mkdir -p "${CELLRANGER_OUTPUT}"
    mkdir -p "${CELLBENDER_OUTPUT}"
    mkdir -p "${FINAL_OUTPUT}"
}

# Get batch number for a patient index (1-indexed)
get_batch_number() {
    local patient_index=$1
    echo $(( (patient_index - 1) / BATCH_SIZE + 1 ))
}

# Get patients for a specific batch
get_batch_patients() {
    local batch_num=$1
    local start_idx=$(( (batch_num - 1) * BATCH_SIZE + 1 ))
    local end_idx=$(( batch_num * BATCH_SIZE ))
    echo "${start_idx}-${end_idx}"
}

# Check scratch usage
check_scratch_usage() {
    local usage=$(du -sh "${SCRATCH_ROOT}" 2>/dev/null | cut -f1)
    echo "Current scratch usage: ${usage}"
}

