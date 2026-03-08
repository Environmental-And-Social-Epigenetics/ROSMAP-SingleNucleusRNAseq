#!/bin/bash
#
# Path Configuration for ROSMAP snRNA-seq Pipeline
# Edit these paths to match your cluster environment.
#
# This file is the SINGLE SOURCE OF TRUTH for all paths used across the
# pipeline.  Every shell wrapper should `source` it rather than hard-coding
# paths.  See config/README.md for setup instructions.
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
# CONDA CONFIGURATION
# =============================================================================

# Conda initialization script
export CONDA_INIT_SCRIPT="/om2/user/mabdel03/anaconda/etc/profile.d/conda.sh"

# Base conda environment directory
export CONDA_ENV_BASE="/orcd/data/lhtsai/001/om2/mabdel03/conda_envs"

# Specific environments
export CELLBENDER_ENV="/om2/user/mabdel03/conda_envs/Cellbender_env"
export SYNAPSE_ENV="${CONDA_ENV_BASE}/synapse_env"
export PYTHON_ENV="${CONDA_ENV_BASE}/python_env"
export SINGLECELL_ENV="${CONDA_ENV_BASE}/single_cell_BP"
export BATCHCORR_ENV="${CONDA_ENV_BASE}/BatchCorrection_SingleCell"
export QC_ENV="/om2/user/mabdel03/conda_envs/qcEnv"

# =============================================================================
# BASE DATA PATHS
# =============================================================================

# Permanent storage (home/project directory)
export DATA_ROOT="/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis"

# Scratch/temporary storage (Openmind scratch)
export SCRATCH_ROOT="/om/scratch/Mon/mabdel03"

# =============================================================================
# REFERENCE DATA
# =============================================================================

# Cell Ranger reference transcriptome
export CELLRANGER_REF="/om2/user/mabdel03/yard/references/human/refdata-gex-GRCh38-2020-A"

# Cell Ranger executable
export CELLRANGER_PATH="/om2/user/mabdel03/apps/yard/cellranger-8.0.0"

# =============================================================================
# DEJAGER DATA PATHS
# =============================================================================

# Input FASTQs (downloaded from Synapse)
export DEJAGER_FASTQS="${SCRATCH_ROOT}/FASTQs"

# Cell Ranger output counts
export DEJAGER_COUNTS="${SCRATCH_ROOT}/Counts"

# Preprocessed counts (CellBender output) - permanent storage
export DEJAGER_PREPROCESSED="${DATA_ROOT}/Data/DeJager/Preprocessed_Counts"

# =============================================================================
# TSAI DATA PATHS
# =============================================================================

# Tsai FASTQs CSV (self-contained in repo)
export TSAI_FASTQS_CSV="${REPO_ROOT}/Data/Tsai/All_ROSMAP_FASTQs.csv"

# Tsai FASTQs (Openmind transfer destination)
export TSAI_FASTQS_DIR="${SCRATCH_ROOT}/Tsai_Data/FASTQs"

# Cell Ranger output (per projid)
export TSAI_CELLRANGER_OUTPUT="${SCRATCH_ROOT}/Tsai_Data/Cellranger_Outputs"

# Temporary Cell Ranger output (scratch - deleted after CellBender)
export TSAI_CELLRANGER_SCRATCH="${TSAI_CELLRANGER_OUTPUT}"

# Temporary CellBender output (scratch - copied to permanent then deleted)
export TSAI_CELLBENDER_SCRATCH="${SCRATCH_ROOT}/Tsai/Cellbender_Output"

# Preprocessed counts (CellBender output)
export TSAI_PREPROCESSED="${SCRATCH_ROOT}/Tsai_Data/Cellbender_Outputs"

# Legacy cohort-specific paths (for reference)
export TSAI_COUNTS_ACE="${SCRATCH_ROOT}/Tsai/ACE/Counts"
export TSAI_COUNTS_RESILIENT="${SCRATCH_ROOT}/Tsai/Resilient/Counts"
export TSAI_COUNTS_SOCISL="${SCRATCH_ROOT}/Tsai/SocIsl/Counts"

# =============================================================================
# TSAI PROCESSING PIPELINE PATHS
# =============================================================================

# Base output directory for the three-stage processing pipeline
export TSAI_PROCESSING_OUTPUTS="${WORKSPACE_ROOT}/Tsai_Data/Processing_Outputs"

# Per-stage I/O directories
export TSAI_QC_FILTERED="${TSAI_PROCESSING_OUTPUTS}/01_QC_Filtered"
export TSAI_DOUBLET_REMOVED="${TSAI_PROCESSING_OUTPUTS}/02_Doublet_Removed"
export TSAI_INTEGRATED="${TSAI_PROCESSING_OUTPUTS}/03_Integrated"

# SLURM log directory
export TSAI_PROCESSING_LOGS="${TSAI_PROCESSING_OUTPUTS}/Logs"

# Patient metadata (tracked in repo)
export TSAI_METADATA_CSV="${REPO_ROOT}/Preprocessing/Tsai/02_Cellranger_Counts/Tracking/patient_metadata.csv"

# Marker reference for ORA annotation (tracked in repo)
export TSAI_MARKERS_RDS="${REPO_ROOT}/Processing/Tsai/Pipeline/Resources/Brain_Human_PFC_Markers_Mohammadi2020.rds"

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Initialize conda
init_conda() {
    source "${CONDA_INIT_SCRIPT}"
}

# Activate a conda environment
activate_env() {
    local env_name="$1"
    init_conda
    conda activate "${CONDA_ENV_BASE}/${env_name}"
}

# Verify paths exist
check_paths() {
    echo "Checking required paths..."
    
    if [[ ! -d "${DATA_ROOT}" ]]; then
        echo "ERROR: DATA_ROOT does not exist: ${DATA_ROOT}"
        return 1
    fi
    
    if [[ "${SCRATCH_ROOT}" == "/path/to/your/scratch" ]]; then
        echo "WARNING: SCRATCH_ROOT is not configured. Update config/paths.sh"
        return 1
    fi
    
    if [[ ! -f "${CONDA_INIT_SCRIPT}" ]]; then
        echo "ERROR: Conda init script not found: ${CONDA_INIT_SCRIPT}"
        return 1
    fi
    
    echo "All paths verified."
    return 0
}

