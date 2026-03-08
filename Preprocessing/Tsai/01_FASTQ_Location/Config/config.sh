#!/bin/bash
#
# Configuration file for FASTQ Transfer Pipeline
#
# This config sources the main config/paths.sh and adds FASTQ-transfer-specific
# variables on top.
#

# =============================================================================
# LOAD GLOBAL PATHS
# =============================================================================

CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${CONFIG_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

# =============================================================================
# FASTQ TRANSFER PIPELINE PATHS
# =============================================================================

# Root directory for this pipeline
export PIPELINE_ROOT="${TSAI_FASTQ_TRANSFER_ROOT}"

# Old directory containing original CSVs
export OLD_DIR="${PIPELINE_ROOT}/Old"

# New directory for outputs
export NEW_DIR="${PIPELINE_ROOT}/New"

# =============================================================================
# INPUT CSVs (from original pipeline)
# =============================================================================

# Master mapping CSV: projid -> Library_ID -> Tsai_path
export MASTER_CSV="${OLD_DIR}/CSVs/Tsai_To_Openmind.csv"

# Backup locations CSV
export BACKUPS_CSV="${OLD_DIR}/CSVs/Reformatted_Backups.csv"

# =============================================================================
# OUTPUT DIRECTORIES
# =============================================================================

# Directory for generated CSVs
export CSV_OUTPUT_DIR="${NEW_DIR}/CSVs"

# Directory for organized FASTQs (symlinks by patient)
export FASTQ_OUTPUT_DIR="${NEW_DIR}/Output/FASTQs_By_Patient"

# Directory for logs
export LOG_DIR="${NEW_DIR}/Logs"

# Temporary directory for intermediate files
export TEMP_DIR="${NEW_DIR}/Logs/tmp"

# =============================================================================
# OUTPUT FILES
# =============================================================================

# Master CSV containing all FASTQ file paths
export ALL_FASTQS_CSV="${CSV_OUTPUT_DIR}/All_ROSMAP_FASTQs.csv"

# Validation report
export VALIDATION_REPORT="${CSV_OUTPUT_DIR}/validation_report.txt"

# =============================================================================
# PARALLELIZATION PARAMETERS
# =============================================================================

# Number of parallel jobs for GNU parallel
export NUM_PARALLEL_JOBS=32

# SLURM job parameters
export SLURM_TIME="10:00:00"
export SLURM_MEM="64G"
export SLURM_CPUS=32

# =============================================================================
# HELPER FUNCTION: Ensure directories exist
# =============================================================================
ensure_dirs() {
    mkdir -p "${CSV_OUTPUT_DIR}"
    mkdir -p "${FASTQ_OUTPUT_DIR}"
    mkdir -p "${LOG_DIR}"
    mkdir -p "${TEMP_DIR}"
}
