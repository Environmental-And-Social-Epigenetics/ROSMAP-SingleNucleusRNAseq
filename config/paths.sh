#!/bin/bash
#
# Path Configuration for ROSMAP snRNA-seq Pipeline
# Edit these paths to match your cluster environment
#

# =============================================================================
# CONDA CONFIGURATION
# =============================================================================

# Conda initialization script
export CONDA_INIT_SCRIPT="/orcd/data/lhtsai/001/om2/mabdel03/miniforge3/etc/profile.d/conda.sh"

# Base conda environment directory
export CONDA_ENV_BASE="/orcd/data/lhtsai/001/om2/mabdel03/conda_envs"

# Specific environments
export CELLBENDER_ENV="${CONDA_ENV_BASE}/Cellbender_env"
export SYNAPSE_ENV="${CONDA_ENV_BASE}/synapse_env"
export PYTHON_ENV="${CONDA_ENV_BASE}/python_env"
export SINGLECELL_ENV="${CONDA_ENV_BASE}/single_cell_BP"
export BATCHCORR_ENV="${CONDA_ENV_BASE}/BatchCorrection_SingleCell"
export QC_ENV="${CONDA_ENV_BASE}/qcEnv"

# =============================================================================
# BASE DATA PATHS
# =============================================================================

# Permanent storage (home/project directory)
export DATA_ROOT="/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis"

# Scratch/temporary storage - UPDATE THIS FOR YOUR CLUSTER
# Note: The old OpenMind paths were /om/scratch/Mon/mabdel03 or /om/scratch/Sun/mabdel03
# These need to be updated to your cluster's scratch filesystem
export SCRATCH_ROOT="/path/to/your/scratch"  # <-- UPDATE THIS

# =============================================================================
# REFERENCE DATA
# =============================================================================

# Cell Ranger reference transcriptome
export CELLRANGER_REF="/orcd/data/lhtsai/001/om2/mabdel03/yard/references/human/refdata-gex-GRCh38-2020-A"

# Cell Ranger executable
export CELLRANGER_PATH="/orcd/data/lhtsai/001/om2/mabdel03/apps/yard/cellranger-8.0.0"

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

# Tsai FASTQs are on Engaging filesystem - paths come from CSV
# Output counts (per cohort)
export TSAI_COUNTS_ACE="${SCRATCH_ROOT}/Tsai/ACE/Counts"
export TSAI_COUNTS_RESILIENT="${SCRATCH_ROOT}/Tsai/Resilient/Counts"
export TSAI_COUNTS_SOCISL="${SCRATCH_ROOT}/Tsai/SocIsl/Counts"

# Preprocessed counts (CellBender output) - permanent storage
export TSAI_PREPROCESSED="${DATA_ROOT}/Data/Tsai/Preprocessing/Preprocessed_Counts"

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

