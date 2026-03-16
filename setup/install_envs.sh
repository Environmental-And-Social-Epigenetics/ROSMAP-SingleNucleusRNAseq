#!/bin/bash
#
# Create conda environments for the ROSMAP snRNA-seq pipeline.
#
# Usage:
#   ./setup/install_envs.sh                  # Processing environments only
#   ./setup/install_envs.sh --analysis       # Processing + Analysis environments
#   ./setup/install_envs.sh --preprocessing  # Processing + Preprocessing environments
#   ./setup/install_envs.sh --all            # All environments
#   ./setup/install_envs.sh /path/to/envs    # Override install location
#
# After running, the script prints the paths to put in config/paths.sh
# (or config/paths.local.sh) so the pipeline can find the environments.
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

source "${REPO_ROOT}/config/paths.sh"

# Parse flags
INSTALL_ANALYSIS=false
INSTALL_PREPROCESSING=false
ENV_BASE_OVERRIDE=""

for arg in "$@"; do
    case "${arg}" in
        --analysis)       INSTALL_ANALYSIS=true ;;
        --preprocessing)  INSTALL_PREPROCESSING=true ;;
        --all)            INSTALL_ANALYSIS=true; INSTALL_PREPROCESSING=true ;;
        --help|-h)
            echo "Usage: $0 [OPTIONS] [/path/to/envs]"
            echo ""
            echo "Options:"
            echo "  --analysis        Also install Analysis environments (DEG, SCENIC, COMPASS, GSEA)"
            echo "  --preprocessing   Also install Preprocessing environments (CellBender, Synapse, bcftools, Globus)"
            echo "  --all             Install all environments"
            echo "  --help            Show this help"
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

echo "=============================================="
echo "ROSMAP snRNA-seq — Conda Environment Setup"
echo "=============================================="
echo ""
echo "Install location: ${ENV_BASE}"
echo ""

create_env() {
    local yml="$1"
    local path="$2"
    local label="$3"

    echo "  Creating ${label} ..."
    if [[ -d "${path}" ]]; then
        echo "    Already exists at ${path} — skipping (delete to recreate)"
    else
        conda env create -f "${yml}" -p "${path}" -q
        echo "    Created: ${path}"
    fi
}

# =====================================================================
# PROCESSING ENVIRONMENTS (always installed)
# =====================================================================
echo "── Processing Pipeline Environments ──"
echo ""

QC_PATH="${ENV_BASE}/qcEnv"
create_env "${PROC_ENVS_DIR}/stage1_qc.yml" "${QC_PATH}" "Stage 1: QC filtering (qcEnv)"

SC_PATH="${ENV_BASE}/single_cell_BP"
create_env "${PROC_ENVS_DIR}/stage2_doublets.yml" "${SC_PATH}" "Stage 2: Doublet removal (single_cell_BP)"

BC_PATH="${ENV_BASE}/BatchCorrection_SingleCell"
create_env "${PROC_ENVS_DIR}/stage3_integration.yml" "${BC_PATH}" "Stage 3: Integration (BatchCorrection_SingleCell)"

echo ""

# =====================================================================
# PREPROCESSING ENVIRONMENTS (optional)
# =====================================================================
if [[ "${INSTALL_PREPROCESSING}" == "true" ]]; then
    echo "── Preprocessing Environments ──"
    echo ""

    CB_PATH="${ENV_BASE}/Cellbender_env"
    create_env "${PREPROC_ENVS_DIR}/cellbender.yml" "${CB_PATH}" "CellBender (GPU)"

    SYN_PATH="${ENV_BASE}/synapse_env"
    create_env "${PREPROC_ENVS_DIR}/synapse.yml" "${SYN_PATH}" "Synapse download"

    BCF_PATH="${ENV_BASE}/bcftools_env"
    create_env "${PREPROC_ENVS_DIR}/bcftools.yml" "${BCF_PATH}" "bcftools"

    GLB_PATH="${ENV_BASE}/globus_env"
    create_env "${PREPROC_ENVS_DIR}/globus.yml" "${GLB_PATH}" "Globus CLI"

    echo ""
fi

# =====================================================================
# ANALYSIS ENVIRONMENTS (optional)
# =====================================================================
if [[ "${INSTALL_ANALYSIS}" == "true" ]]; then
    echo "── Analysis Environments ──"
    echo ""

    DEG_PATH="${ENV_BASE}/deg_analysis"
    create_env "${ANALYSIS_ENVS_DIR}/deg.yml" "${DEG_PATH}" "DEG analysis (DESeq2/limma)"

    SCENIC_PATH="${ENV_BASE}/scenic_analysis"
    create_env "${ANALYSIS_ENVS_DIR}/scenic.yml" "${SCENIC_PATH}" "SCENIC (pySCENIC)"

    COMPASS_PATH="${ENV_BASE}/compass_analysis"
    create_env "${ANALYSIS_ENVS_DIR}/compass.yml" "${COMPASS_PATH}" "COMPASS (metabolic)"

    GSEA_PATH="${ENV_BASE}/gsea_analysis"
    create_env "${ANALYSIS_ENVS_DIR}/gsea.yml" "${GSEA_PATH}" "GSEA (WebGestaltR)"

    echo ""
fi

# =====================================================================
# SUMMARY
# =====================================================================
echo "=============================================="
echo "Done!  Add the following to config/paths.local.sh"
echo "(or update config/paths.sh directly):"
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

echo ""
