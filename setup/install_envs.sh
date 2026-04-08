#!/bin/bash
#
# Create and validate conda environments for the ROSMAP snRNA-seq pipeline.
#
# Usage:
#   ./setup/install_envs.sh
#   ./setup/install_envs.sh --analysis
#   ./setup/install_envs.sh --preprocessing
#   ./setup/install_envs.sh --all
#   ./setup/install_envs.sh --method=requirements
#   ./setup/install_envs.sh --analysis --recreate
#   ./setup/install_envs.sh /path/to/envs
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

source "${REPO_ROOT}/config/paths.sh"

INSTALL_ANALYSIS=false
INSTALL_PREPROCESSING=false
INSTALL_METHOD="conda"
ENV_BASE_OVERRIDE=""
RECREATE=false

for arg in "$@"; do
    case "${arg}" in
        --analysis)       INSTALL_ANALYSIS=true ;;
        --preprocessing)  INSTALL_PREPROCESSING=true ;;
        --all)            INSTALL_ANALYSIS=true; INSTALL_PREPROCESSING=true ;;
        --method=conda)         INSTALL_METHOD="conda" ;;
        --method=requirements)  INSTALL_METHOD="requirements" ;;
        --recreate)       RECREATE=true ;;
        --help|-h)
            echo "Usage: $0 [OPTIONS] [/path/to/envs]"
            echo ""
            echo "Options:"
            echo "  --analysis                Also install Analysis environments"
            echo "  --preprocessing           Also install Preprocessing environments"
            echo "  --all                     Install all environments"
            echo "  --method=conda            Create envs from environment.yml (default)"
            echo "  --method=requirements     Bootstrap envs then install requirements.txt"
            echo "  --recreate                Delete and recreate existing envs before validation"
            echo "  --help                    Show this help"
            exit 0
            ;;
        -*)
            echo "Unknown flag: ${arg}. Use --help for usage."
            exit 1
            ;;
        *)
            ENV_BASE_OVERRIDE="${arg}"
            ;;
    esac
done

ENV_BASE="${ENV_BASE_OVERRIDE:-${CONDA_ENV_BASE}}"
PROC_ENVS_DIR="${REPO_ROOT}/Processing/Tsai/Pipeline/envs"
ANALYSIS_ENVS_DIR="${REPO_ROOT}/Analysis/envs"
PREPROC_ENVS_DIR="${REPO_ROOT}/Preprocessing/envs"

if [[ -f "${CONDA_INIT_SCRIPT}" ]]; then
    # shellcheck disable=SC1090
    source "${CONDA_INIT_SCRIPT}"
fi

if ! command -v conda &>/dev/null; then
    echo "ERROR: conda is not available. Check CONDA_INIT_SCRIPT or your PATH."
    exit 1
fi

echo "=============================================="
echo "ROSMAP snRNA-seq — Conda Environment Setup"
echo "=============================================="
echo ""
echo "Install location: ${ENV_BASE}"
echo "Install method:   ${INSTALL_METHOD}"
echo "Recreate envs:    ${RECREATE}"
echo ""

verify_python_imports() {
    local path="$1"
    shift
    local python_bin="${path}/bin/python"
    if [[ ! -x "${python_bin}" ]]; then
        echo "ERROR: Python not found in ${path}"
        exit 1
    fi
    local import_name
    for import_name in "$@"; do
        echo "    Verifying Python import: ${import_name}"
        "${python_bin}" -c "import ${import_name}" >/dev/null
    done
}

verify_r_packages() {
    local path="$1"
    shift
    local rscript_bin="${path}/bin/Rscript"
    if [[ ! -x "${rscript_bin}" ]]; then
        echo "ERROR: Rscript not found in ${path}"
        exit 1
    fi
    local pkg
    for pkg in "$@"; do
        echo "    Verifying R package: ${pkg}"
        "${rscript_bin}" -e "library(${pkg})" >/dev/null
    done
}

verify_commands() {
    local path="$1"
    shift
    local cmd
    for cmd in "$@"; do
        echo "    Verifying command: ${cmd}"
        if [[ ! -x "${path}/bin/${cmd}" ]]; then
            echo "ERROR: command '${cmd}' not found in ${path}/bin"
            exit 1
        fi
    done
}

run_post_install() {
    local env_id="$1"
    local path="$2"

    case "${env_id}" in
        sccomp)
            echo "    Installing CRAN package: sccomp"
            "${path}/bin/Rscript" -e "dir.create('${path}/lib/R/library', recursive=TRUE, showWarnings=FALSE); install.packages('sccomp', repos='https://cloud.r-project.org', lib='${path}/lib/R/library')"
            ;;
    esac
}

verify_env() {
    local env_id="$1"
    local path="$2"

    echo "    Running post-install verification for ${env_id}"
    case "${env_id}" in
        stage1_qc)
            verify_python_imports "${path}" scanpy anndata pandas h5py
            ;;
        stage2_doublets)
            verify_r_packages "${path}" scDblFinder BiocParallel zellkonverter SingleCellExperiment ggplot2
            ;;
        stage3_integration)
            verify_python_imports "${path}" scanpy anndata pandas h5py harmonypy decoupler rpy2 scib sklearn
            ;;
        cellbender)
            verify_python_imports "${path}" torch
            verify_commands "${path}" cellbender
            ;;
        synapse)
            verify_python_imports "${path}" pandas synapseclient
            ;;
        bcftools)
            verify_commands "${path}" bcftools samtools
            ;;
        globus)
            verify_commands "${path}" globus
            ;;
        deg)
            verify_r_packages "${path}" DESeq2 edgeR limma ggplot2 dplyr tidyr
            ;;
        scenic)
            verify_python_imports "${path}" pandas numpy loompy pyscenic scanpy anndata
            ;;
        compass)
            verify_python_imports "${path}" pandas numpy scipy scanpy anndata
            verify_commands "${path}" compass
            ;;
        gsea)
            verify_r_packages "${path}" clusterProfiler org.Hs.eg.db WebGestaltR ggplot2 dplyr
            ;;
        nebula)
            verify_python_imports "${path}" anndata h5py pandas numpy scipy
            verify_r_packages "${path}" DESeq2 zellkonverter SingleCellExperiment scran dplyr
            ;;
        sccomp)
            verify_r_packages "${path}" sccomp readr ggplot2 dplyr tidyr patchwork
            ;;
        *)
            echo "ERROR: Unknown environment id '${env_id}'"
            exit 1
            ;;
    esac
}

remove_env() {
    local path="$1"

    if [[ ! -d "${path}" ]]; then
        return
    fi

    echo "    Removing existing env at ${path}"
    if [[ -d "${path}/conda-meta" ]] || [[ -d "${path}/conda_meta" ]]; then
        if ! conda env remove -y -p "${path}" >/dev/null 2>&1; then
            rm -rf "${path}"
        fi
    else
        rm -rf "${path}"
    fi
}

create_env() {
    local env_id="$1"
    local env_dir="$2"
    local path="$3"
    local label="$4"
    local pattern="$5"
    shift 5
    local bootstrap_args=("$@")

    local env_yml="${env_dir}/environment.yml"
    local req_txt="${env_dir}/requirements.txt"

    echo "  ${label}"

    if [[ ! -f "${env_yml}" ]]; then
        echo "ERROR: Missing environment.yml: ${env_yml}"
        exit 1
    fi
    if [[ ! -f "${req_txt}" ]]; then
        echo "ERROR: Missing requirements.txt: ${req_txt}"
        exit 1
    fi

    if [[ -d "${path}" ]]; then
        if [[ "${RECREATE}" == "true" ]]; then
            remove_env "${path}"
        else
            echo "    Already exists at ${path}; validating current env"
            verify_env "${env_id}" "${path}"
            echo "    Validation passed: ${path}"
            return
        fi
    fi

    if [[ "${INSTALL_METHOD}" == "conda" ]]; then
        conda env create -f "${env_yml}" -p "${path}" -q
    else
        conda create -y -p "${path}" "${bootstrap_args[@]}"
        "${path}/bin/python" -m pip install -r "${req_txt}"
        if [[ "${pattern}" == "hybrid" ]]; then
            echo "    Hybrid environment bootstrapped. See ${env_dir}/README.md for any extra runtime steps."
        fi
    fi

    run_post_install "${env_id}" "${path}"
    verify_env "${env_id}" "${path}"
    echo "    Created and validated: ${path}"
}

echo "── Processing Pipeline Environments ──"
echo ""

QC_PATH="${ENV_BASE}/qcEnv"
create_env "stage1_qc" \
    "${PROC_ENVS_DIR}/stage1_qc" \
    "${QC_PATH}" \
    "Stage 1: QC filtering (qcEnv)" \
    "requirements-complete" \
    python=3.10 pip

SC_PATH="${ENV_BASE}/single_cell_BP"
create_env "stage2_doublets" \
    "${PROC_ENVS_DIR}/stage2_doublets" \
    "${SC_PATH}" \
    "Stage 2: Doublet removal (single_cell_BP)" \
    "hybrid" \
    -c conda-forge -c bioconda -c defaults \
    python=3.10 pip r-base=4.2 r-ggplot2 \
    bioconductor-biocparallel bioconductor-genomeinfodbdata \
    bioconductor-singlecellexperiment bioconductor-scdblfinder bioconductor-zellkonverter

BC_PATH="${ENV_BASE}/BatchCorrection_SingleCell"
create_env "stage3_integration" \
    "${PROC_ENVS_DIR}/stage3_integration" \
    "${BC_PATH}" \
    "Stage 3: Integration (BatchCorrection_SingleCell)" \
    "hybrid" \
    -c conda-forge -c defaults \
    python=3.10 pip r-base=4.2 rpy2>=3.5

echo ""

if [[ "${INSTALL_PREPROCESSING}" == "true" ]]; then
    echo "── Preprocessing Environments ──"
    echo ""

    CB_PATH="${ENV_BASE}/Cellbender_env"
    create_env "cellbender" \
        "${PREPROC_ENVS_DIR}/cellbender" \
        "${CB_PATH}" \
        "CellBender (GPU)" \
        "hybrid" \
        -c conda-forge -c pytorch -c defaults \
        python=3.10 pip pytorch>=2.0 pytorch-cuda>=11.8

    SYN_PATH="${ENV_BASE}/synapse_env"
    create_env "synapse" \
        "${PREPROC_ENVS_DIR}/synapse" \
        "${SYN_PATH}" \
        "Synapse download" \
        "requirements-complete" \
        python=3.10 pip

    BCF_PATH="${ENV_BASE}/bcftools_env"
    create_env "bcftools" \
        "${PREPROC_ENVS_DIR}/bcftools" \
        "${BCF_PATH}" \
        "bcftools" \
        "hybrid" \
        -c bioconda -c conda-forge -c defaults \
        python=3.10 pip bcftools>=1.17 samtools>=1.17 htslib>=1.17

    GLB_PATH="${ENV_BASE}/globus_env"
    create_env "globus" \
        "${PREPROC_ENVS_DIR}/globus" \
        "${GLB_PATH}" \
        "Globus CLI" \
        "requirements-complete" \
        python=3.10 pip

    echo ""
fi

if [[ "${INSTALL_ANALYSIS}" == "true" ]]; then
    echo "── Analysis Environments ──"
    echo ""

    DEG_PATH="${ENV_BASE}/deg_analysis"
    create_env "deg" \
        "${ANALYSIS_ENVS_DIR}/deg" \
        "${DEG_PATH}" \
        "DEG analysis (DESeq2/limma)" \
        "hybrid" \
        -c conda-forge -c bioconda -c defaults \
        python=3.10 pip r-base=4.2 \
        bioconductor-deseq2 bioconductor-edger bioconductor-limma \
        bioconductor-biocparallel r-ggplot2 r-dplyr r-tidyr r-pheatmap

    SCENIC_PATH="${ENV_BASE}/scenic_analysis"
    create_env "scenic" \
        "${ANALYSIS_ENVS_DIR}/scenic" \
        "${SCENIC_PATH}" \
        "SCENIC (pySCENIC)" \
        "requirements-complete" \
        python=3.10 pip

    COMPASS_PATH="${ENV_BASE}/compass_analysis"
    create_env "compass" \
        "${ANALYSIS_ENVS_DIR}/compass" \
        "${COMPASS_PATH}" \
        "COMPASS (metabolic)" \
        "requirements-complete" \
        python=3.10 pip

    GSEA_PATH="${ENV_BASE}/gsea_analysis"
    create_env "gsea" \
        "${ANALYSIS_ENVS_DIR}/gsea" \
        "${GSEA_PATH}" \
        "GSEA (WebGestaltR)" \
        "hybrid" \
        -c conda-forge -c bioconda -c defaults \
        python=3.10 pip r-base=4.2 r-ggplot2 r-dplyr \
        bioconductor-clusterprofiler bioconductor-org.hs.eg.db r-webgestaltr

    NEBULA_PATH="${ENV_BASE}/nebulaAnalysis7"
    create_env "nebula" \
        "${ANALYSIS_ENVS_DIR}/nebula" \
        "${NEBULA_PATH}" \
        "ACE DEG / pseudobulk (nebulaAnalysis7)" \
        "hybrid" \
        -c conda-forge -c bioconda -c defaults \
        python=3.10 pip r-base=4.2 r-dplyr bioconductor-zellkonverter \
        bioconductor-deseq2 bioconductor-singlecellexperiment bioconductor-scran \
        pandas>=2.0 numpy>=1.24 scipy>=1.10 h5py>=3.8

    SCCOMP_PATH="${ENV_BASE}/sccompAnalysis"
    create_env "sccomp" \
        "${ANALYSIS_ENVS_DIR}/sccomp" \
        "${SCCOMP_PATH}" \
        "ACE cell type proportion (sccompAnalysis)" \
        "hybrid" \
        -c conda-forge -c defaults \
        python=3.10 pip r-base=4.2 r-dplyr r-tidyr r-readr r-ggplot2 r-patchwork

    echo ""
fi

echo "=============================================="
echo "Done! Add the following to config/paths.local.sh"
echo "(or rely on CONDA_ENV_BASE defaults):"
echo "=============================================="
echo ""
echo "export CONDA_ENV_BASE=\"${ENV_BASE}\""
echo "export QC_ENV=\"${QC_PATH}\""
echo "export SINGLECELL_ENV=\"${SC_PATH}\""
echo "export BATCHCORR_ENV=\"${BC_PATH}\""

if [[ "${INSTALL_PREPROCESSING}" == "true" ]]; then
    echo "export CELLBENDER_ENV=\"${CB_PATH}\""
    echo "export SYNAPSE_ENV=\"${SYN_PATH}\""
    echo "export BCFTOOLS_ENV=\"${BCF_PATH}\""
    echo "export GLOBUS_ENV=\"${GLB_PATH}\""
fi

if [[ "${INSTALL_ANALYSIS}" == "true" ]]; then
    echo "export DEG_ANALYSIS_ENV=\"${DEG_PATH}\""
    echo "export SCENIC_ANALYSIS_ENV=\"${SCENIC_PATH}\""
    echo "export COMPASS_ANALYSIS_ENV=\"${COMPASS_PATH}\""
    echo "export GSEA_ANALYSIS_ENV=\"${GSEA_PATH}\""
    echo "export NEBULA_ENV=\"${NEBULA_PATH}\""
    echo "export SCCOMP_ENV=\"${SCCOMP_PATH}\""
fi

echo ""
