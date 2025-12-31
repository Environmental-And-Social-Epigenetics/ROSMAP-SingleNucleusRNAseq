#!/bin/bash
#
# Configuration file for FASTQ Transfer Pipeline
# All paths and parameters are defined here for easy modification
#

# =============================================================================
# CONDA CONFIGURATION
# =============================================================================

# Conda initialization script
export CONDA_INIT_SCRIPT="/orcd/data/lhtsai/001/om2/mabdel03/miniforge3/etc/profile.d/conda.sh"

# Conda environment with pandas for Python scripts
export CONDA_ENV="/orcd/data/lhtsai/001/om2/mabdel03/conda_envs/python_data_analysis"

# =============================================================================
# BASE PATHS
# =============================================================================

# Root directory for this pipeline
export PIPELINE_ROOT="/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/Tsai/Preprocessing/FASTQ_Transfer"

# Old directory containing original CSVs
export OLD_DIR="${PIPELINE_ROOT}/Old"

# New directory for outputs
export NEW_DIR="${PIPELINE_ROOT}/New"

# =============================================================================
# INPUT CSVs (from original pipeline)
# =============================================================================

# Master mapping CSV: projid -> Library_ID -> Tsai_path (797 rows, 483 projids)
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
export SLURM_MAIL_USER="mabdel03@mit.edu"

# =============================================================================
# HELPER FUNCTION: Ensure directories exist
# =============================================================================
ensure_dirs() {
    mkdir -p "${CSV_OUTPUT_DIR}"
    mkdir -p "${FASTQ_OUTPUT_DIR}"
    mkdir -p "${LOG_DIR}"
    mkdir -p "${TEMP_DIR}"
}

