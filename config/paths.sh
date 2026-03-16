#!/bin/bash
#
# Path Configuration for ROSMAP snRNA-seq Pipeline
#
# This file is the SINGLE SOURCE OF TRUTH for all paths used across the
# pipeline.  Every shell wrapper should `source` it rather than hard-coding
# paths.  See config/README.md for setup instructions.
#
# HOW TO CONFIGURE:
#   1. Copy config/paths.local.sh.template to config/paths.local.sh
#   2. Edit config/paths.local.sh with your cluster-specific paths
#   3. Run: source config/paths.sh && check_paths
#
# Variables marked __UNCONFIGURED__ MUST be set in paths.local.sh before
# running the pipeline.  They will cause immediate, visible errors if used
# as-is (the string is an invalid filesystem path by design).
#

# =============================================================================
# AUTO-DETECTED PATHS (no edits needed)
# =============================================================================

# Repository root — derived from this file's location (config/ is one level
# below the repo root).  Works regardless of where the repo is cloned.
_CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export REPO_ROOT="$(cd "${_CONFIG_DIR}/.." && pwd)"

# Workspace root — the parent of the repo, where large data directories live
# (e.g. Tsai_Data/).
export WORKSPACE_ROOT="$(cd "${REPO_ROOT}/.." && pwd)"

# =============================================================================
# LOCAL OVERRIDES (loaded early so they take precedence)
# =============================================================================
# To customize paths without modifying this tracked file, create
# config/paths.local.sh and override any variables there.  That file
# is listed in .gitignore and will never be committed.
#
# A template is provided: config/paths.local.sh.template
if [[ -f "${_CONFIG_DIR}/paths.local.sh" ]]; then
    source "${_CONFIG_DIR}/paths.local.sh"
fi

# =============================================================================
# CONDA CONFIGURATION
# =============================================================================

# Conda initialization script — update in paths.local.sh to match your install:
#   miniconda:  ~/miniconda3/etc/profile.d/conda.sh
#   miniforge:  ~/miniforge3/etc/profile.d/conda.sh
#   anaconda:   ~/anaconda3/etc/profile.d/conda.sh
export CONDA_INIT_SCRIPT="${CONDA_INIT_SCRIPT:-${HOME}/miniforge3/etc/profile.d/conda.sh}"

# Base conda environment directory
export CONDA_ENV_BASE="${CONDA_ENV_BASE:-${HOME}/conda_envs}"

# Specific environments (all derive from CONDA_ENV_BASE by default)
export CELLBENDER_ENV="${CELLBENDER_ENV:-${CONDA_ENV_BASE}/Cellbender_env}"
export SYNAPSE_ENV="${SYNAPSE_ENV:-${CONDA_ENV_BASE}/synapse_env}"
export PYTHON_ENV="${PYTHON_ENV:-${CONDA_ENV_BASE}/python_env}"
export SINGLECELL_ENV="${SINGLECELL_ENV:-${CONDA_ENV_BASE}/single_cell_BP}"
export BATCHCORR_ENV="${BATCHCORR_ENV:-${CONDA_ENV_BASE}/BatchCorrection_SingleCell}"
export QC_ENV="${QC_ENV:-${CONDA_ENV_BASE}/qcEnv}"
export GLOBUS_ENV="${GLOBUS_ENV:-${CONDA_ENV_BASE}/globus_env}"

# =============================================================================
# BASE DATA PATHS
# =============================================================================

# Permanent storage (home/project directory) — MUST be set in paths.local.sh
export DATA_ROOT="${DATA_ROOT:-__UNCONFIGURED__set_DATA_ROOT_in_paths_local_sh}"

# Scratch/temporary storage — MUST be set in paths.local.sh
export SCRATCH_ROOT="${SCRATCH_ROOT:-__UNCONFIGURED__set_SCRATCH_ROOT_in_paths_local_sh}"

# =============================================================================
# REFERENCE DATA
# =============================================================================

# Cell Ranger reference transcriptome — MUST be set in paths.local.sh
export CELLRANGER_REF="${CELLRANGER_REF:-__UNCONFIGURED__set_CELLRANGER_REF_in_paths_local_sh}"

# Cell Ranger executable — MUST be set in paths.local.sh
export CELLRANGER_PATH="${CELLRANGER_PATH:-__UNCONFIGURED__set_CELLRANGER_PATH_in_paths_local_sh}"

# =============================================================================
# DEJAGER DATA PATHS
# =============================================================================

# Input FASTQs (downloaded from Synapse)
export DEJAGER_FASTQS="${DEJAGER_FASTQS:-${SCRATCH_ROOT}/FASTQs}"

# Cell Ranger output counts
export DEJAGER_COUNTS="${DEJAGER_COUNTS:-${SCRATCH_ROOT}/Counts}"

# Preprocessed counts (CellBender output) - permanent storage
export DEJAGER_PREPROCESSED="${DEJAGER_PREPROCESSED:-${DATA_ROOT}/Data/DeJager/Preprocessed_Counts}"

# =============================================================================
# DEJAGER DEMUXLET PATHS
# =============================================================================

# WGS data directory (filtered BAMs, pileups, demuxlet outputs)
# Only needed for DeJager preprocessing — set in paths.local.sh if running Demuxlet
export DEJAGER_WGS_DIR="${DEJAGER_WGS_DIR:-__UNCONFIGURED__set_DEJAGER_WGS_DIR_in_paths_local_sh}"

# SNP VCF for demuxlet (SNP-only, concatenated, lifted to GRCh38, ~84GB)
export DEJAGER_DEMUX_VCF="${DEJAGER_DEMUX_VCF:-${DEJAGER_WGS_DIR}/snp_fixedconcatenated_liftedROSMAP.vcf.gz}"

# Demuxafy Singularity container
export DEMUXAFY_SIF="${DEMUXAFY_SIF:-${DEJAGER_WGS_DIR}/Demuxafy.sif}"

# Per-library patient ID files
export DEJAGER_PATIENT_IDS_DIR="${DEJAGER_PATIENT_IDS_DIR:-${DEJAGER_WGS_DIR}/individ}"

# Bcftools conda environment (for BAM filtering)
export BCFTOOLS_ENV="${BCFTOOLS_ENV:-${CONDA_ENV_BASE}/bcftools_env}"

# Singularity module (cluster-specific; override in paths.local.sh if needed)
export SINGULARITY_MODULE="${SINGULARITY_MODULE:-openmind/singularity/3.10.4}"

# =============================================================================
# CANONICAL DATA PATHS (under repo Data/Transcriptomics/)
# =============================================================================
#
# These are the in-repo data locations. After cloning, populate them using the
# scripts in Data_Access/ (download from NAS, Globus transfer, or Synapse).
# Data is also copied here from legacy scattered locations via
# Data/Transcriptomics/copy_data.sbatch.
#
# If your data lives on separate scratch/fast storage, you can either:
#   (a) Symlink these directories to your actual data locations, or
#   (b) Override the WORKING PATHS below in paths.local.sh

export TSAI_FASTQS="${TSAI_FASTQS:-${REPO_ROOT}/Data/Transcriptomics/Tsai/FASTQs}"
export TSAI_CELLRANGER="${TSAI_CELLRANGER:-${REPO_ROOT}/Data/Transcriptomics/Tsai/Cellranger_Output}"
export TSAI_CELLBENDER="${TSAI_CELLBENDER:-${REPO_ROOT}/Data/Transcriptomics/Tsai/Cellbender_Output}"
export DEJAGER_FASTQS_DIR="${DEJAGER_FASTQS_DIR:-${REPO_ROOT}/Data/Transcriptomics/DeJager/FASTQs}"
export DEJAGER_CELLRANGER="${DEJAGER_CELLRANGER:-${REPO_ROOT}/Data/Transcriptomics/DeJager/Cellranger_Output}"
export DEJAGER_CELLBENDER="${DEJAGER_CELLBENDER:-${REPO_ROOT}/Data/Transcriptomics/DeJager/Cellbender_Output}"
export TSAI_PROCESSING="${TSAI_PROCESSING:-${REPO_ROOT}/Data/Transcriptomics/Tsai/Processing_Outputs}"
export DEJAGER_PROCESSING="${DEJAGER_PROCESSING:-${REPO_ROOT}/Data/Transcriptomics/DeJager/Processing_Outputs}"

# =============================================================================
# TSAI WORKING PATHS (for active processing on scratch/fast storage)
# =============================================================================
#
# These paths may point to faster scratch storage during active processing.
# Override in paths.local.sh if your scratch differs from the defaults.
# On a simple single-filesystem setup, these default to the canonical paths.

# Tsai FASTQs CSV (generated by Preprocessing/Tsai/01_FASTQ_Location pipeline)
export TSAI_FASTQS_CSV="${TSAI_FASTQS_CSV:-${DATA_ROOT}/Data/Tsai/Preprocessing/FASTQ_Transfer/New/CSVs/All_ROSMAP_FASTQs.csv}"

# Tsai FASTQs (Openmind transfer destination)
export TSAI_FASTQS_DIR="${TSAI_FASTQS_DIR:-${SCRATCH_ROOT}/Tsai_Data/FASTQs}"

# Cell Ranger output (per projid)
export TSAI_CELLRANGER_OUTPUT="${TSAI_CELLRANGER_OUTPUT:-${SCRATCH_ROOT}/Tsai_Data/Cellranger_Outputs}"

# Temporary Cell Ranger output (scratch - deleted after CellBender)
export TSAI_CELLRANGER_SCRATCH="${TSAI_CELLRANGER_SCRATCH:-${TSAI_CELLRANGER_OUTPUT}}"

# Temporary CellBender output (scratch - copied to permanent then deleted)
export TSAI_CELLBENDER_SCRATCH="${TSAI_CELLBENDER_SCRATCH:-${SCRATCH_ROOT}/Tsai/Cellbender_Output}"

# Preprocessed counts (CellBender output — used as input to Processing pipeline)
export TSAI_PREPROCESSED="${TSAI_PREPROCESSED:-${WORKSPACE_ROOT}/Tsai_Data/Cellbender_Outputs}"

# Legacy cohort-specific paths (DEPRECATED — kept for backward compatibility only)
export TSAI_COUNTS_ACE="${TSAI_COUNTS_ACE:-${SCRATCH_ROOT}/Tsai/ACE/Counts}"
export TSAI_COUNTS_RESILIENT="${TSAI_COUNTS_RESILIENT:-${SCRATCH_ROOT}/Tsai/Resilient/Counts}"
export TSAI_COUNTS_SOCISL="${TSAI_COUNTS_SOCISL:-${SCRATCH_ROOT}/Tsai/SocIsl/Counts}"

# =============================================================================
# TSAI PROCESSING PIPELINE PATHS
# =============================================================================

# Base output directory for the three-stage processing pipeline
export TSAI_PROCESSING_OUTPUTS="${TSAI_PROCESSING_OUTPUTS:-${WORKSPACE_ROOT}/Tsai_Data/Processing_Outputs}"

# Per-stage I/O directories
export TSAI_QC_FILTERED="${TSAI_QC_FILTERED:-${TSAI_PROCESSING_OUTPUTS}/01_QC_Filtered}"
export TSAI_DOUBLET_REMOVED="${TSAI_DOUBLET_REMOVED:-${TSAI_PROCESSING_OUTPUTS}/02_Doublet_Removed}"
export TSAI_INTEGRATED="${TSAI_INTEGRATED:-${TSAI_PROCESSING_OUTPUTS}/03_Integrated}"

# SLURM log directory
export TSAI_PROCESSING_LOGS="${TSAI_PROCESSING_LOGS:-${TSAI_PROCESSING_OUTPUTS}/Logs}"

# Patient metadata (tracked in repo)
export TSAI_METADATA_CSV="${TSAI_METADATA_CSV:-${REPO_ROOT}/Preprocessing/Tsai/02_Cellranger_Counts/Tracking/patient_metadata.csv}"

# Marker reference for ORA annotation (tracked in repo)
export TSAI_MARKERS_RDS="${TSAI_MARKERS_RDS:-${REPO_ROOT}/Processing/Tsai/Pipeline/Resources/Brain_Human_PFC_Markers_Mohammadi2020.rds}"

# =============================================================================
# DEJAGER PROCESSING PIPELINE PATHS
# =============================================================================

# Base output directory for the three-stage processing pipeline
export DEJAGER_PROCESSING_OUTPUTS="${DEJAGER_PROCESSING_OUTPUTS:-${WORKSPACE_ROOT}/DeJager_Data/Processing_Outputs}"

# Per-stage I/O directories
export DEJAGER_QC_FILTERED="${DEJAGER_QC_FILTERED:-${DEJAGER_PROCESSING_OUTPUTS}/01_QC_Filtered}"
export DEJAGER_DOUBLET_REMOVED="${DEJAGER_DOUBLET_REMOVED:-${DEJAGER_PROCESSING_OUTPUTS}/02_Doublet_Removed}"
export DEJAGER_INTEGRATED="${DEJAGER_INTEGRATED:-${DEJAGER_PROCESSING_OUTPUTS}/03_Integrated}"

# SLURM log directory
export DEJAGER_PROCESSING_LOGS="${DEJAGER_PROCESSING_LOGS:-${DEJAGER_PROCESSING_OUTPUTS}/Logs}"

# Patient barcode-to-ID mapping CSV (155MB, too large for git — see
# Processing/DeJager/Pipeline/README.md for how to obtain this file)
export DEJAGER_PATIENT_MAP="${DEJAGER_PATIENT_MAP:-${DEJAGER_WGS_DIR}/cell_to_patient_assignmentsFinal0.csv}"

# Patient ID overrides for "alone" libraries (JSON, tracked in repo)
export DEJAGER_PATIENT_ID_OVERRIDES="${DEJAGER_PATIENT_ID_OVERRIDES:-${REPO_ROOT}/Processing/DeJager/Pipeline/Resources/patient_id_overrides.json}"

# Marker reference for ORA annotation (shared with Tsai)
export DEJAGER_MARKERS_RDS="${DEJAGER_MARKERS_RDS:-${REPO_ROOT}/Processing/Tsai/Pipeline/Resources/Brain_Human_PFC_Markers_Mohammadi2020.rds}"

# =============================================================================
# PHENOTYPE / CLINICAL DATA (tracked in repo)
# =============================================================================

export PHENOTYPE_DIR="${PHENOTYPE_DIR:-${REPO_ROOT}/Data/Phenotypes}"
export ROSMAP_CLINICAL_CSV="${ROSMAP_CLINICAL_CSV:-${PHENOTYPE_DIR}/ROSMAP_clinical.csv}"
export ACE_SCORES_CSV="${ACE_SCORES_CSV:-${PHENOTYPE_DIR}/TSAI_DEJAGER_all_patients_wACEscores.csv}"
export DEJAGER_ID_MAP="${DEJAGER_ID_MAP:-${PHENOTYPE_DIR}/DeJager_ID_Map.csv}"

# =============================================================================
# SLURM SETTINGS
# =============================================================================

# Email for SLURM job notifications (leave empty to disable)
export SLURM_MAIL_USER="${SLURM_MAIL_USER:-}"

# Default SLURM partitions (leave empty for cluster default)
export SLURM_PARTITION="${SLURM_PARTITION:-}"
export SLURM_PARTITION_GPU="${SLURM_PARTITION_GPU:-}"

# =============================================================================
# NAS / TRANSFER SETTINGS
# =============================================================================

# NAS SFTP user (defaults to $USER; override in paths.local.sh if different)
export NAS_SFTP_USER="${NAS_SFTP_USER:-${USER}}"

# =============================================================================
# FASTQ TRANSFER PIPELINE PATHS (Tsai)
# =============================================================================

# Root directory for the FASTQ transfer pipeline outputs
export TSAI_FASTQ_TRANSFER_ROOT="${TSAI_FASTQ_TRANSFER_ROOT:-${DATA_ROOT}/Data/Tsai/Preprocessing/FASTQ_Transfer}"

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Initialize conda
init_conda() {
    source "${CONDA_INIT_SCRIPT}"
}

# Activate a conda environment by full path
activate_env() {
    local env_path="$1"
    init_conda
    conda activate "${env_path}"
}

# Verify paths exist and detect unconfigured sentinels
check_paths() {
    local errors=0
    local warnings=0
    echo "Checking required paths..."
    echo ""

    # Check for unconfigured sentinel values
    local sentinel_vars=(
        "DATA_ROOT"
        "SCRATCH_ROOT"
        "CELLRANGER_REF"
        "CELLRANGER_PATH"
    )
    for var in "${sentinel_vars[@]}"; do
        if [[ "${!var}" == *"__UNCONFIGURED__"* ]]; then
            echo "ERROR: ${var} is not configured."
            echo "       Create config/paths.local.sh (see config/paths.local.sh.template)"
            echo "       and set: export ${var}=\"/your/path/here\""
            echo ""
            errors=$((errors + 1))
        fi
    done

    # Check DEJAGER_WGS_DIR separately (only needed for DeJager preprocessing)
    if [[ "${DEJAGER_WGS_DIR}" == *"__UNCONFIGURED__"* ]]; then
        echo "NOTE:  DEJAGER_WGS_DIR is not configured (only needed for DeJager Demuxlet step)."
        echo ""
    fi

    # Validate configured paths exist
    if [[ "${DATA_ROOT}" != *"__UNCONFIGURED__"* && ! -d "${DATA_ROOT}" ]]; then
        echo "ERROR: DATA_ROOT does not exist: ${DATA_ROOT}"
        errors=$((errors + 1))
    fi

    if [[ "${SCRATCH_ROOT}" != *"__UNCONFIGURED__"* && ! -d "${SCRATCH_ROOT}" ]]; then
        echo "WARNING: SCRATCH_ROOT does not exist: ${SCRATCH_ROOT}"
        warnings=$((warnings + 1))
    fi

    if [[ ! -f "${CONDA_INIT_SCRIPT}" ]]; then
        echo "ERROR: Conda init script not found: ${CONDA_INIT_SCRIPT}"
        errors=$((errors + 1))
    fi

    if [[ "${CELLRANGER_PATH}" != *"__UNCONFIGURED__"* && ! -f "${CELLRANGER_PATH}/cellranger" ]]; then
        echo "WARNING: Cell Ranger executable not found: ${CELLRANGER_PATH}/cellranger"
        warnings=$((warnings + 1))
    fi

    if [[ "${CELLRANGER_REF}" != *"__UNCONFIGURED__"* && ! -d "${CELLRANGER_REF}" ]]; then
        echo "WARNING: Cell Ranger reference not found: ${CELLRANGER_REF}"
        warnings=$((warnings + 1))
    fi

    echo ""
    if [[ ${errors} -gt 0 ]]; then
        echo "FAILED: ${errors} critical error(s). See above for details."
        return 1
    elif [[ ${warnings} -gt 0 ]]; then
        echo "PASSED with ${warnings} warning(s)."
        return 0
    fi

    echo "All paths verified."
    return 0
}
