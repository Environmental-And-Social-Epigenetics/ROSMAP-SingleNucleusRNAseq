#!/bin/bash
#
# Preflight Check System for ROSMAP snRNA-seq Pipeline
#
# Verifies that all prerequisites (conda environments, config variables, input
# files, software) are in place before submitting pipeline jobs.
#
# Usage (standalone):
#   bash config/preflight.sh tsai-stage1       # Check Tsai Processing Stage 1
#   bash config/preflight.sh dejager-04-demuxlet  # Check DeJager Demuxlet step
#   bash config/preflight.sh all               # Check everything
#
# Usage (sourced from submit_pipeline.sh):
#   source config/preflight.sh
#   preflight_check tsai-stage1 || exit 1
#
# Exit codes: 0=PASS, 1=FAIL (critical), 2=WARN (non-critical)
#

# =============================================================================
# INTERNAL STATE
# =============================================================================

_PF_ERRORS=0
_PF_WARNINGS=0

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

# Map dataset identifier to directory name (handles DeJager capitalization)
_pf_dataset_dir() {
    case "$1" in
        tsai)    echo "Tsai" ;;
        dejager) echo "DeJager" ;;
        *)       echo "$1" ;;
    esac
}

_pf_pass() {
    echo -e "  [\033[32mPASS\033[0m] $1"
}

_pf_fail() {
    echo -e "  [\033[31mFAIL\033[0m] $1"
    echo "         -> $2"
    _PF_ERRORS=$((_PF_ERRORS + 1))
}

_pf_warn() {
    echo -e "  [\033[33mWARN\033[0m] $1"
    echo "         -> $2"
    _PF_WARNINGS=$((_PF_WARNINGS + 1))
}

_pf_check_var() {
    local var_name="$1"
    if [[ -z "${!var_name:-}" ]]; then
        _pf_fail "Config variable ${var_name} is not set" \
                 "Check config/paths.sh or config/paths.local.sh"
    else
        _pf_pass "Config variable ${var_name} is set"
    fi
}

_pf_check_dir() {
    local dir_path="$1"
    local label="${2:-$1}"
    if [[ -d "${dir_path}" ]]; then
        _pf_pass "${label} exists"
    else
        _pf_fail "${label} not found: ${dir_path}" \
                 "Create the directory or update the path in config/paths.sh"
    fi
}

_pf_check_file() {
    local file_path="$1"
    local label="${2:-$1}"
    if [[ -f "${file_path}" ]]; then
        _pf_pass "${label} exists"
    else
        _pf_fail "${label} not found: ${file_path}" \
                 "Ensure the file exists or update the path in config/paths.sh"
    fi
}

_pf_check_conda_env() {
    local env_path="$1"
    local env_name="${2:-$1}"
    if [[ ! -d "${env_path}" ]]; then
        _pf_fail "Conda env '${env_name}' not found: ${env_path}" \
                 "Run: conda env create -f <spec.yml> -p ${env_path}"
        return 1
    fi
    if [[ ! -d "${env_path}/conda-meta" ]] && [[ ! -d "${env_path}/conda_meta" ]]; then
        _pf_fail "Directory exists but is not a valid conda env: ${env_path}" \
                 "Remove and recreate: conda env create -f <spec.yml> -p ${env_path}"
        return 1
    fi
    _pf_pass "Conda env '${env_name}' exists"
    return 0
}

_pf_check_python_pkg() {
    local env_path="$1"
    local pkg="$2"
    local python_bin="${env_path}/bin/python"

    if [[ ! -x "${python_bin}" ]]; then
        _pf_fail "Python not found in env: ${python_bin}" \
                 "Recreate the conda environment"
        return 1
    fi

    local version
    version=$("${python_bin}" -c "import ${pkg}; print(getattr(${pkg}, '__version__', 'unknown'))" 2>/dev/null)
    if [[ $? -ne 0 ]]; then
        _pf_fail "Python package '${pkg}' not importable in ${env_path}" \
                 "Run: ${env_path}/bin/pip install ${pkg}"
        return 1
    fi

    _pf_pass "Python package '${pkg}' (v${version})"
    return 0
}

_pf_check_r_pkg() {
    local env_path="$1"
    local pkg="$2"
    local extra_lib="${3:-}"
    local rscript_bin="${env_path}/bin/Rscript"

    if [[ ! -x "${rscript_bin}" ]]; then
        _pf_fail "Rscript not found in env: ${rscript_bin}" \
                 "Recreate the conda environment with r-base"
        return 1
    fi

    local r_cmd="library(${pkg})"
    if [[ -n "${extra_lib}" ]]; then
        r_cmd=".libPaths(c('${extra_lib}', .libPaths())); library(${pkg})"
    fi

    if "${rscript_bin}" -e "${r_cmd}" &>/dev/null; then
        _pf_pass "R package '${pkg}' loadable"
    else
        _pf_fail "R package '${pkg}' not loadable in ${env_path}" \
                 "Install via conda or R within the environment"
    fi
}

_pf_check_command() {
    local cmd="$1"
    local label="${2:-$1}"
    if [[ -x "${cmd}" ]] || command -v "${cmd}" &>/dev/null; then
        _pf_pass "${label} found"
    else
        _pf_fail "${label} not found: ${cmd}" \
                 "Install or update the path"
    fi
}

_pf_check_gpu() {
    if command -v nvidia-smi &>/dev/null; then
        _pf_pass "nvidia-smi available (GPU present)"
    else
        _pf_warn "nvidia-smi not found (GPU not available on this node)" \
                 "Expected on login nodes; GPU will be available on compute nodes"
    fi
}

_pf_check_module() {
    local mod="$1"
    # Module system may not be initialized on login nodes
    if ! command -v module &>/dev/null && [[ -f /etc/profile.d/modules.sh ]]; then
        source /etc/profile.d/modules.sh 2>/dev/null
    fi
    if command -v module &>/dev/null; then
        if module avail "${mod}" 2>&1 | grep -q "${mod}"; then
            _pf_pass "Module '${mod}' available"
        else
            _pf_fail "Module '${mod}' not available" \
                     "Check 'module avail' or contact system administrator"
        fi
    else
        _pf_warn "Module system not available on this node" \
                 "Module '${mod}' will be checked on compute nodes"
    fi
}

# =============================================================================
# PREPROCESSING — TSAI
# =============================================================================

preflight_tsai_01_fastq_location() {
    echo "--- Tsai 01: FASTQ Location ---"
    _pf_check_var CONDA_INIT_SCRIPT
    _pf_check_file "${CONDA_INIT_SCRIPT}" "Conda init script"
    _pf_check_var TSAI_FASTQ_TRANSFER_ROOT
    _pf_check_var TSAI_FASTQS_CSV
    _pf_check_command "python3" "Python 3"
    _pf_check_command "parallel" "GNU parallel"
    _pf_check_dir "${DATA_ROOT}" "DATA_ROOT"
}

preflight_tsai_02_cellranger_counts() {
    echo "--- Tsai 02: Cell Ranger + CellBender ---"
    _pf_check_var CONDA_INIT_SCRIPT
    _pf_check_file "${CONDA_INIT_SCRIPT}" "Conda init script"
    _pf_check_var CELLRANGER_PATH
    _pf_check_var CELLRANGER_REF
    _pf_check_var TSAI_FASTQS_DIR
    _pf_check_var TSAI_FASTQS_CSV
    _pf_check_var CELLBENDER_ENV

    # Cell Ranger executable and reference
    _pf_check_command "${CELLRANGER_PATH}/cellranger" "Cell Ranger executable"
    _pf_check_dir "${CELLRANGER_REF}" "Cell Ranger reference genome"

    # Cell Ranger version
    if [[ -x "${CELLRANGER_PATH}/cellranger" ]]; then
        local cr_version
        cr_version=$("${CELLRANGER_PATH}/cellranger" --version 2>&1 | grep -oP '\d+\.\d+\.\d+' | head -1)
        if [[ "${cr_version}" == "8.0.0" ]]; then
            _pf_pass "Cell Ranger version: ${cr_version}"
        else
            _pf_warn "Cell Ranger version ${cr_version} (expected 8.0.0)" \
                     "Pipeline was tested with v8.0.0"
        fi
    fi

    # Input FASTQs
    _pf_check_dir "${TSAI_FASTQS_DIR}" "Tsai FASTQs directory"
    _pf_check_file "${TSAI_FASTQS_CSV}" "Tsai FASTQs master CSV"
    if [[ -d "${TSAI_FASTQS_DIR}" ]]; then
        local n_patient_dirs
        n_patient_dirs=$(ls -1d "${TSAI_FASTQS_DIR}"/*/ 2>/dev/null | wc -l)
        if [[ "${n_patient_dirs}" -lt 10 ]]; then
            _pf_warn "Only ${n_patient_dirs} patient directories in TSAI_FASTQS_DIR" \
                     "Expected ~480 patient directories"
        else
            _pf_pass "Found ${n_patient_dirs} patient directories in FASTQs dir"
        fi
    fi

    # CellBender environment
    _pf_check_conda_env "${CELLBENDER_ENV}" "CellBender"

    # GPU (warn only)
    _pf_check_gpu
}

# =============================================================================
# PREPROCESSING — DEJAGER
# =============================================================================

preflight_dejager_01_fastq_download() {
    echo "--- DeJager 01: FASTQ Download ---"
    _pf_check_var CONDA_INIT_SCRIPT
    _pf_check_file "${CONDA_INIT_SCRIPT}" "Conda init script"
    _pf_check_var SYNAPSE_ENV
    _pf_check_var DEJAGER_FASTQS

    # Synapse environment and packages
    _pf_check_conda_env "${SYNAPSE_ENV}" "Synapse"
    if [[ -d "${SYNAPSE_ENV}" ]]; then
        _pf_check_python_pkg "${SYNAPSE_ENV}" "synapseclient"
        _pf_check_python_pkg "${SYNAPSE_ENV}" "pandas"
    fi

    # Synapse credentials
    if [[ -f "${HOME}/.synapseConfig" ]] || [[ -n "${SYNAPSE_AUTH_TOKEN:-}" ]]; then
        _pf_pass "Synapse credentials found"
    else
        _pf_fail "Synapse credentials not found" \
                 "Run 'synapse login --rememberMe' or set SYNAPSE_AUTH_TOKEN"
    fi

    # Output parent directory
    local fastq_parent
    fastq_parent="$(dirname "${DEJAGER_FASTQS}")"
    _pf_check_dir "${fastq_parent}" "DEJAGER_FASTQS parent directory"
}

preflight_dejager_02_cellranger_counts() {
    echo "--- DeJager 02: Cell Ranger Counts ---"
    _pf_check_var CONDA_INIT_SCRIPT
    _pf_check_file "${CONDA_INIT_SCRIPT}" "Conda init script"
    _pf_check_var CELLRANGER_PATH
    _pf_check_var CELLRANGER_REF
    _pf_check_var DEJAGER_FASTQS
    _pf_check_var PYTHON_ENV

    # Cell Ranger
    _pf_check_command "${CELLRANGER_PATH}/cellranger" "Cell Ranger executable"
    _pf_check_dir "${CELLRANGER_REF}" "Cell Ranger reference genome"

    # Python environment (used by Count_DeJager.py)
    _pf_check_conda_env "${PYTHON_ENV}" "Python"
    if [[ -d "${PYTHON_ENV}" ]]; then
        _pf_check_python_pkg "${PYTHON_ENV}" "pandas"
    fi

    # Input FASTQs
    _pf_check_dir "${DEJAGER_FASTQS}" "DeJager FASTQs directory"
    if [[ -d "${DEJAGER_FASTQS}" ]]; then
        local n_lib_dirs
        n_lib_dirs=$(ls -1d "${DEJAGER_FASTQS}"/*/ 2>/dev/null | wc -l)
        if [[ "${n_lib_dirs}" -eq 0 ]]; then
            _pf_fail "No library directories in DEJAGER_FASTQS" \
                     "Run step 01_FASTQ_Download first"
        else
            _pf_pass "Found ${n_lib_dirs} library directories in FASTQs dir"
        fi
    fi
}

preflight_dejager_03_cellbender() {
    echo "--- DeJager 03: CellBender ---"
    _pf_check_var CONDA_INIT_SCRIPT
    _pf_check_file "${CONDA_INIT_SCRIPT}" "Conda init script"
    _pf_check_var CELLBENDER_ENV
    _pf_check_var DEJAGER_COUNTS
    _pf_check_var DEJAGER_PREPROCESSED

    # CellBender environment
    _pf_check_conda_env "${CELLBENDER_ENV}" "CellBender"

    # Input: Cell Ranger outputs
    _pf_check_dir "${DEJAGER_COUNTS}" "DeJager Cell Ranger counts"
    if [[ -d "${DEJAGER_COUNTS}" ]]; then
        local n_counts
        n_counts=$(ls -1d "${DEJAGER_COUNTS}"/*/ 2>/dev/null | wc -l)
        if [[ "${n_counts}" -eq 0 ]]; then
            _pf_fail "No library directories in DEJAGER_COUNTS" \
                     "Run step 02_Cellranger_Counts first"
        else
            _pf_pass "Found ${n_counts} library directories in counts dir"
        fi
    fi

    # GPU (warn only)
    _pf_check_gpu
}

preflight_dejager_04_demuxlet() {
    echo "--- DeJager 04: Demuxlet/Freemuxlet ---"
    _pf_check_var CONDA_INIT_SCRIPT
    _pf_check_file "${CONDA_INIT_SCRIPT}" "Conda init script"
    _pf_check_var BCFTOOLS_ENV
    _pf_check_var SINGULARITY_MODULE
    _pf_check_var DEMUXAFY_SIF
    _pf_check_var DEJAGER_DEMUX_VCF
    _pf_check_var DEJAGER_PATIENT_IDS_DIR
    _pf_check_var DEJAGER_COUNTS
    _pf_check_var DEJAGER_PREPROCESSED

    # Bcftools environment
    _pf_check_conda_env "${BCFTOOLS_ENV}" "bcftools"

    # Singularity module
    _pf_check_module "${SINGULARITY_MODULE}"

    # Demuxafy container
    _pf_check_file "${DEMUXAFY_SIF}" "Demuxafy Singularity container"

    # WGS VCF
    _pf_check_file "${DEJAGER_DEMUX_VCF}" "WGS VCF for demuxlet"

    # Patient ID files
    _pf_check_dir "${DEJAGER_PATIENT_IDS_DIR}" "Patient IDs directory"
    if [[ -d "${DEJAGER_PATIENT_IDS_DIR}" ]]; then
        local n_id_files
        n_id_files=$(ls -1 "${DEJAGER_PATIENT_IDS_DIR}"/individPat*.txt 2>/dev/null | wc -l)
        if [[ "${n_id_files}" -eq 0 ]]; then
            _pf_fail "No individPat*.txt files in DEJAGER_PATIENT_IDS_DIR" \
                     "Patient ID files are required for demuxlet"
        else
            _pf_pass "Found ${n_id_files} patient ID files"
        fi
    fi

    # Cell Ranger BAMs
    _pf_check_dir "${DEJAGER_COUNTS}" "DeJager Cell Ranger counts (for BAMs)"

    # CellBender outputs
    _pf_check_dir "${DEJAGER_PREPROCESSED}" "DeJager CellBender outputs"

    # Helper tools
    local helper_dir="${REPO_ROOT}/Preprocessing/DeJager/04_Demuxlet_Freemuxlet/popscle_helper_tools"
    _pf_check_dir "${helper_dir}" "popscle helper tools"
    if [[ -d "${helper_dir}" ]]; then
        _pf_check_file "${helper_dir}/filter_bam_file_for_popscle_dsc_pileup.sh" "BAM filter script"
        _pf_check_file "${helper_dir}/filter_vcf_file_for_popscle.sh" "VCF filter script"
        _pf_check_file "${helper_dir}/sort_vcf_same_as_bam.sh" "VCF sort script"
    fi
}

# =============================================================================
# PROCESSING PIPELINE (shared between Tsai and DeJager)
# =============================================================================

preflight_processing_stage1() {
    local dataset="${1:?Usage: preflight_processing_stage1 tsai|dejager}"
    echo "--- $(_pf_dataset_dir "${dataset}") Processing Stage 1: QC Filter ---"

    # Common
    _pf_check_var CONDA_INIT_SCRIPT
    _pf_check_file "${CONDA_INIT_SCRIPT}" "Conda init script"
    _pf_check_var QC_ENV

    # QC conda environment + packages
    _pf_check_conda_env "${QC_ENV}" "QC (qcEnv)"
    if [[ -d "${QC_ENV}" ]]; then
        _pf_check_python_pkg "${QC_ENV}" "scanpy"
        _pf_check_python_pkg "${QC_ENV}" "anndata"
        _pf_check_python_pkg "${QC_ENV}" "pandas"
        _pf_check_python_pkg "${QC_ENV}" "scipy"
        _pf_check_python_pkg "${QC_ENV}" "h5py"
    fi

    # Dataset-specific paths
    if [[ "${dataset}" == "tsai" ]]; then
        _pf_check_var TSAI_PREPROCESSED
        _pf_check_dir "${TSAI_PREPROCESSED}" "Tsai CellBender outputs (input)"
        _pf_check_var TSAI_METADATA_CSV
        _pf_check_file "${TSAI_METADATA_CSV}" "Tsai patient metadata CSV"

        if [[ -d "${TSAI_PREPROCESSED}" ]]; then
            local n_inputs
            n_inputs=$(find "${TSAI_PREPROCESSED}" -maxdepth 2 -name "*.h5" 2>/dev/null | head -20 | wc -l)
            if [[ "${n_inputs}" -eq 0 ]]; then
                _pf_fail "No .h5 files found in TSAI_PREPROCESSED" \
                         "Run CellBender preprocessing first"
            else
                _pf_pass "Found CellBender .h5 input files"
            fi
        fi
    elif [[ "${dataset}" == "dejager" ]]; then
        _pf_check_var DEJAGER_PREPROCESSED
        _pf_check_dir "${DEJAGER_PREPROCESSED}" "DeJager CellBender outputs (input)"
        _pf_check_var DEJAGER_PATIENT_MAP
        _pf_check_file "${DEJAGER_PATIENT_MAP}" "DeJager patient map CSV"
        _pf_check_var DEJAGER_PATIENT_ID_OVERRIDES
        _pf_check_file "${DEJAGER_PATIENT_ID_OVERRIDES}" "DeJager patient ID overrides JSON"

        if [[ -d "${DEJAGER_PREPROCESSED}" ]]; then
            local n_inputs
            n_inputs=$(find "${DEJAGER_PREPROCESSED}" -maxdepth 2 -name "*.h5" 2>/dev/null | head -20 | wc -l)
            if [[ "${n_inputs}" -eq 0 ]]; then
                _pf_fail "No .h5 files found in DEJAGER_PREPROCESSED" \
                         "Run CellBender preprocessing first"
            else
                _pf_pass "Found CellBender .h5 input files"
            fi
        fi
    fi

    # Pipeline scripts
    local pipeline_dir="${REPO_ROOT}/Processing/$(_pf_dataset_dir "${dataset}")/Pipeline"
    _pf_check_file "${pipeline_dir}/01_qc_filter.py" "Stage 1 Python script"
    _pf_check_file "${pipeline_dir}/01_qc_filter.sh" "Stage 1 shell wrapper"
}

preflight_processing_stage2() {
    local dataset="${1:?Usage: preflight_processing_stage2 tsai|dejager}"
    echo "--- $(_pf_dataset_dir "${dataset}") Processing Stage 2: Doublet Removal ---"

    # Common
    _pf_check_var CONDA_INIT_SCRIPT
    _pf_check_file "${CONDA_INIT_SCRIPT}" "Conda init script"
    _pf_check_var SINGLECELL_ENV

    # SingleCell R environment + packages
    # The R script adds Pipeline/R_libs to .libPaths() at runtime, so we
    # must include it here to match the actual package resolution.
    local pipeline_dir="${REPO_ROOT}/Processing/$(_pf_dataset_dir "${dataset}")/Pipeline"
    local r_libs_extra="${pipeline_dir}/R_libs"

    _pf_check_conda_env "${SINGLECELL_ENV}" "SingleCell (single_cell_BP)"
    if [[ -d "${SINGLECELL_ENV}" ]]; then
        _pf_check_r_pkg "${SINGLECELL_ENV}" "scDblFinder" "${r_libs_extra}"
        _pf_check_r_pkg "${SINGLECELL_ENV}" "BiocParallel" "${r_libs_extra}"
        _pf_check_r_pkg "${SINGLECELL_ENV}" "zellkonverter" "${r_libs_extra}"
        _pf_check_r_pkg "${SINGLECELL_ENV}" "SingleCellExperiment" "${r_libs_extra}"
        _pf_check_r_pkg "${SINGLECELL_ENV}" "ggplot2" "${r_libs_extra}"
    fi

    # Dataset-specific: input from Stage 1
    if [[ "${dataset}" == "tsai" ]]; then
        _pf_check_var TSAI_QC_FILTERED
        _pf_check_dir "${TSAI_QC_FILTERED}" "Tsai Stage 1 outputs (input)"
        _pf_check_var TSAI_METADATA_CSV
        _pf_check_file "${TSAI_METADATA_CSV}" "Tsai patient metadata CSV"

        if [[ -d "${TSAI_QC_FILTERED}" ]]; then
            local n_inputs
            n_inputs=$(ls -1 "${TSAI_QC_FILTERED}"/*_qc.h5ad 2>/dev/null | wc -l)
            if [[ "${n_inputs}" -eq 0 ]]; then
                _pf_fail "No *_qc.h5ad files in TSAI_QC_FILTERED" \
                         "Run Stage 1 (QC filter) first"
            else
                _pf_pass "Found ${n_inputs} QC-filtered .h5ad files"
            fi
        fi
    elif [[ "${dataset}" == "dejager" ]]; then
        _pf_check_var DEJAGER_QC_FILTERED
        _pf_check_dir "${DEJAGER_QC_FILTERED}" "DeJager Stage 1 outputs (input)"

        if [[ -d "${DEJAGER_QC_FILTERED}" ]]; then
            local n_inputs
            n_inputs=$(ls -1 "${DEJAGER_QC_FILTERED}"/*_qc.h5ad 2>/dev/null | wc -l)
            if [[ "${n_inputs}" -eq 0 ]]; then
                _pf_fail "No *_qc.h5ad files in DEJAGER_QC_FILTERED" \
                         "Run Stage 1 (QC filter) first"
            else
                _pf_pass "Found ${n_inputs} QC-filtered .h5ad files"
            fi
        fi
    fi

    # Pipeline scripts
    _pf_check_file "${pipeline_dir}/02_doublet_removal.Rscript" "Stage 2 R script"
    _pf_check_file "${pipeline_dir}/02_doublet_removal.sh" "Stage 2 shell wrapper"

    # GenomeInfoDbData (optional local supplement)
    if [[ -d "${r_libs_extra}/GenomeInfoDbData" ]]; then
        _pf_pass "GenomeInfoDbData available in R_libs"
    else
        _pf_warn "GenomeInfoDbData not in ${r_libs_extra}" \
                 "Will be bootstrapped on first run, or install: R -e 'BiocManager::install(\"GenomeInfoDbData\")'"
    fi
}

preflight_processing_stage3() {
    local dataset="${1:?Usage: preflight_processing_stage3 tsai|dejager}"
    echo "--- $(_pf_dataset_dir "${dataset}") Processing Stage 3: Integration & Annotation ---"

    # Common
    _pf_check_var CONDA_INIT_SCRIPT
    _pf_check_file "${CONDA_INIT_SCRIPT}" "Conda init script"
    _pf_check_var BATCHCORR_ENV

    # BatchCorrection environment + packages
    _pf_check_conda_env "${BATCHCORR_ENV}" "BatchCorrection"
    if [[ -d "${BATCHCORR_ENV}" ]]; then
        _pf_check_python_pkg "${BATCHCORR_ENV}" "scanpy"
        _pf_check_python_pkg "${BATCHCORR_ENV}" "anndata"
        _pf_check_python_pkg "${BATCHCORR_ENV}" "harmonypy"
        _pf_check_python_pkg "${BATCHCORR_ENV}" "decoupler"
        _pf_check_python_pkg "${BATCHCORR_ENV}" "rpy2"
        _pf_check_python_pkg "${BATCHCORR_ENV}" "pandas"
    fi

    # Dataset-specific: input from Stage 2 + markers
    if [[ "${dataset}" == "tsai" ]]; then
        _pf_check_var TSAI_DOUBLET_REMOVED
        _pf_check_dir "${TSAI_DOUBLET_REMOVED}" "Tsai Stage 2 outputs (input)"
        _pf_check_var TSAI_METADATA_CSV
        _pf_check_file "${TSAI_METADATA_CSV}" "Tsai patient metadata CSV"
        _pf_check_var TSAI_MARKERS_RDS
        _pf_check_file "${TSAI_MARKERS_RDS}" "Tsai markers RDS file"

        if [[ -d "${TSAI_DOUBLET_REMOVED}" ]]; then
            local n_inputs
            n_inputs=$(ls -1 "${TSAI_DOUBLET_REMOVED}"/*_singlets.h5ad 2>/dev/null | wc -l)
            if [[ "${n_inputs}" -eq 0 ]]; then
                _pf_fail "No *_singlets.h5ad files in TSAI_DOUBLET_REMOVED" \
                         "Run Stage 2 (doublet removal) first"
            else
                _pf_pass "Found ${n_inputs} doublet-removed .h5ad files"
            fi
        fi
    elif [[ "${dataset}" == "dejager" ]]; then
        _pf_check_var DEJAGER_DOUBLET_REMOVED
        _pf_check_dir "${DEJAGER_DOUBLET_REMOVED}" "DeJager Stage 2 outputs (input)"
        _pf_check_var DEJAGER_MARKERS_RDS
        _pf_check_file "${DEJAGER_MARKERS_RDS}" "DeJager markers RDS file"

        if [[ -d "${DEJAGER_DOUBLET_REMOVED}" ]]; then
            local n_inputs
            n_inputs=$(ls -1 "${DEJAGER_DOUBLET_REMOVED}"/*_singlets.h5ad 2>/dev/null | wc -l)
            if [[ "${n_inputs}" -eq 0 ]]; then
                _pf_fail "No *_singlets.h5ad files in DEJAGER_DOUBLET_REMOVED" \
                         "Run Stage 2 (doublet removal) first"
            else
                _pf_pass "Found ${n_inputs} doublet-removed .h5ad files"
            fi
        fi
    fi

    # Pipeline scripts
    local pipeline_dir="${REPO_ROOT}/Processing/$(_pf_dataset_dir "${dataset}")/Pipeline"
    _pf_check_file "${pipeline_dir}/03_integration_annotation.py" "Stage 3 Python script"
    _pf_check_file "${pipeline_dir}/03_integration_annotation.sh" "Stage 3 shell wrapper"
}

# =============================================================================
# DISPATCHER
# =============================================================================

preflight_check() {
    local step="${1:?Usage: preflight_check STEP_NAME}"
    shift

    _PF_ERRORS=0
    _PF_WARNINGS=0

    echo ""
    echo "=========================================="
    echo "  PREFLIGHT CHECK: ${step}"
    echo "=========================================="
    echo ""

    # Common checks
    _pf_check_var REPO_ROOT
    _pf_check_dir "${REPO_ROOT}" "Repository root"
    _pf_check_file "${REPO_ROOT}/config/paths.sh" "Config file (paths.sh)"
    echo ""

    case "${step}" in
        # Preprocessing
        tsai-01-fastq)            preflight_tsai_01_fastq_location ;;
        tsai-02-cellranger)       preflight_tsai_02_cellranger_counts ;;
        dejager-01-download)      preflight_dejager_01_fastq_download ;;
        dejager-02-cellranger)    preflight_dejager_02_cellranger_counts ;;
        dejager-03-cellbender)    preflight_dejager_03_cellbender ;;
        dejager-04-demuxlet)      preflight_dejager_04_demuxlet ;;

        # Processing (with aliases)
        tsai-stage1|tsai-qc)           preflight_processing_stage1 tsai ;;
        tsai-stage2|tsai-doublets)     preflight_processing_stage2 tsai ;;
        tsai-stage3|tsai-integrate)    preflight_processing_stage3 tsai ;;
        dejager-stage1|dejager-qc)           preflight_processing_stage1 dejager ;;
        dejager-stage2|dejager-doublets)     preflight_processing_stage2 dejager ;;
        dejager-stage3|dejager-integrate)    preflight_processing_stage3 dejager ;;

        # Batch
        all-tsai)
            preflight_tsai_01_fastq_location; echo ""
            preflight_tsai_02_cellranger_counts; echo ""
            preflight_processing_stage1 tsai; echo ""
            preflight_processing_stage2 tsai; echo ""
            preflight_processing_stage3 tsai
            ;;
        all-dejager)
            preflight_dejager_01_fastq_download; echo ""
            preflight_dejager_02_cellranger_counts; echo ""
            preflight_dejager_03_cellbender; echo ""
            preflight_dejager_04_demuxlet; echo ""
            preflight_processing_stage1 dejager; echo ""
            preflight_processing_stage2 dejager; echo ""
            preflight_processing_stage3 dejager
            ;;
        all)
            preflight_tsai_01_fastq_location; echo ""
            preflight_tsai_02_cellranger_counts; echo ""
            preflight_dejager_01_fastq_download; echo ""
            preflight_dejager_02_cellranger_counts; echo ""
            preflight_dejager_03_cellbender; echo ""
            preflight_dejager_04_demuxlet; echo ""
            preflight_processing_stage1 tsai; echo ""
            preflight_processing_stage2 tsai; echo ""
            preflight_processing_stage3 tsai; echo ""
            preflight_processing_stage1 dejager; echo ""
            preflight_processing_stage2 dejager; echo ""
            preflight_processing_stage3 dejager
            ;;
        *)
            echo "ERROR: Unknown step '${step}'"
            echo ""
            echo "Valid steps:"
            echo "  Preprocessing:"
            echo "    tsai-01-fastq          Tsai FASTQ location/transfer"
            echo "    tsai-02-cellranger     Tsai Cell Ranger + CellBender"
            echo "    dejager-01-download    DeJager FASTQ download (Synapse)"
            echo "    dejager-02-cellranger  DeJager Cell Ranger"
            echo "    dejager-03-cellbender  DeJager CellBender"
            echo "    dejager-04-demuxlet    DeJager Demuxlet/Freemuxlet"
            echo ""
            echo "  Processing:"
            echo "    tsai-stage1            Tsai QC filtering"
            echo "    tsai-stage2            Tsai doublet removal"
            echo "    tsai-stage3            Tsai integration & annotation"
            echo "    dejager-stage1         DeJager QC filtering"
            echo "    dejager-stage2         DeJager doublet removal"
            echo "    dejager-stage3         DeJager integration & annotation"
            echo ""
            echo "  Batch:"
            echo "    all-tsai               All Tsai steps"
            echo "    all-dejager            All DeJager steps"
            echo "    all                    Everything"
            return 1
            ;;
    esac

    # Summary
    echo ""
    echo "=========================================="
    if [[ ${_PF_ERRORS} -gt 0 ]]; then
        echo -e "  RESULT: \033[31mFAIL\033[0m (${_PF_ERRORS} error(s), ${_PF_WARNINGS} warning(s))"
        echo "=========================================="
        return 1
    elif [[ ${_PF_WARNINGS} -gt 0 ]]; then
        echo -e "  RESULT: \033[33mPASS with warnings\033[0m (${_PF_WARNINGS} warning(s))"
        echo "=========================================="
        return 0
    else
        echo -e "  RESULT: \033[32mPASS\033[0m"
        echo "=========================================="
        return 0
    fi
}

# =============================================================================
# CLI MODE
# =============================================================================

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    _PREFLIGHT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    source "${_PREFLIGHT_DIR}/paths.sh"

    if [[ $# -eq 0 ]]; then
        echo "Usage: bash config/preflight.sh <step-name>"
        echo "Run 'bash config/preflight.sh all' to check everything."
        echo "Run 'bash config/preflight.sh --help' for valid step names."
        exit 1
    fi

    if [[ "$1" == "--help" ]] || [[ "$1" == "-h" ]]; then
        preflight_check "SHOW_HELP" 2>/dev/null || true
        exit 0
    fi

    preflight_check "$@"
    exit $?
fi
