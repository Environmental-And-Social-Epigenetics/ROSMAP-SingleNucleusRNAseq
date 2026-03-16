#!/bin/bash
#
# CellBender Standalone Pipeline Configuration
# For processing Tsai CellRanger outputs with CellBender
#

# =============================================================================
# LOAD GLOBAL PATHS
# =============================================================================

CELLBENDER_CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${CELLBENDER_CONFIG_DIR}/../../../../config/paths.sh"

# =============================================================================
# PIPELINE DIRECTORIES
# =============================================================================

export PIPELINE_DIR="${REPO_ROOT}/Preprocessing/Tsai/03_Cellbender"
export SCRIPTS_DIR="${PIPELINE_DIR}/Scripts"
export BATCH_SCRIPTS_DIR="${PIPELINE_DIR}/Batch_Scripts"
export LOGS_DIR="${PIPELINE_DIR}/Logs"
export LOGS_OUT="${LOGS_DIR}/Outs"
export LOGS_ERR="${LOGS_DIR}/Errs"
export TRACKING_DIR="${PIPELINE_DIR}/Tracking"

# =============================================================================
# INPUT/OUTPUT PATHS
# =============================================================================

export CELLRANGER_OUTPUT="${TSAI_CELLRANGER_OUTPUT}"
export CELLBENDER_OUTPUT="${TSAI_CELLBENDER_SCRATCH}"
export FINAL_OUTPUT="${TSAI_PREPROCESSED}"

# Tracking from CellRanger pipeline
export CELLRANGER_TRACKING_DIR="${REPO_ROOT}/Preprocessing/Tsai/02_Cellranger_Counts/Tracking"
export BATCH_ASSIGNMENTS_CSV="${CELLRANGER_TRACKING_DIR}/batch_assignments.csv"
export CELLRANGER_COMPLETED="${CELLRANGER_TRACKING_DIR}/cellranger_completed.txt"

# =============================================================================
# CELLBENDER PARAMETERS
# =============================================================================

export CONDA_INIT_SCRIPT="${CONDA_INIT_SCRIPT}"
export CELLBENDER_CONDA="${CELLBENDER_ENV}"
export CB_FPR="0"
export CB_CUDA="true"

# =============================================================================
# SLURM PARAMETERS
# =============================================================================

export CB_SLURM_TIME="47:00:00"
export CB_SLURM_CPUS="32"
export CB_SLURM_MEM="128G"
export CB_SLURM_GPU="1"
export SLURM_MAIL_USER="${SLURM_MAIL_USER}"

# =============================================================================
# HELPERS
# =============================================================================

ensure_pipeline_dirs() {
    mkdir -p "${BATCH_SCRIPTS_DIR}"
    mkdir -p "${LOGS_OUT}"
    mkdir -p "${LOGS_ERR}"
    mkdir -p "${TRACKING_DIR}"
}
