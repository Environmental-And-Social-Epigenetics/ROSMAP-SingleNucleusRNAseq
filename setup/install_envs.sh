#!/bin/bash
#
# Create conda environments for the ROSMAP snRNA-seq Processing pipeline.
#
# Usage:
#   ./setup/install_envs.sh                  # uses CONDA_ENV_BASE from paths.sh
#   ./setup/install_envs.sh /path/to/envs    # override install location
#
# After running, the script prints the paths to put in config/paths.sh
# (or config/paths.local.sh) so the pipeline can find the environments.
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

source "${REPO_ROOT}/config/paths.sh"

# Allow the user to override the install location via argument
ENV_BASE="${1:-${CONDA_ENV_BASE}}"

ENVS_DIR="${REPO_ROOT}/Processing/Tsai/Pipeline/envs"

echo "=============================================="
echo "ROSMAP snRNA-seq — Conda Environment Setup"
echo "=============================================="
echo ""
echo "Install location: ${ENV_BASE}"
echo "YAML specs:       ${ENVS_DIR}"
echo ""

# ── Stage 1: QC filtering ──────────────────────────────────────────────
QC_PATH="${ENV_BASE}/qcEnv"
echo "[1/3] Creating Stage 1 environment (qcEnv) ..."
if [[ -d "${QC_PATH}" ]]; then
    echo "  Already exists at ${QC_PATH} — skipping (delete to recreate)"
else
    conda env create -f "${ENVS_DIR}/stage1_qc.yml" -p "${QC_PATH}" -q
    echo "  Created: ${QC_PATH}"
fi

# ── Stage 2: Doublet removal ───────────────────────────────────────────
SC_PATH="${ENV_BASE}/single_cell_BP"
echo "[2/3] Creating Stage 2 environment (single_cell_BP) ..."
if [[ -d "${SC_PATH}" ]]; then
    echo "  Already exists at ${SC_PATH} — skipping (delete to recreate)"
else
    conda env create -f "${ENVS_DIR}/stage2_doublets.yml" -p "${SC_PATH}" -q
    echo "  Created: ${SC_PATH}"
fi

# ── Stage 3: Integration & annotation ─────────────────────────────────
BC_PATH="${ENV_BASE}/BatchCorrection_SingleCell"
echo "[3/3] Creating Stage 3 environment (BatchCorrection_SingleCell) ..."
if [[ -d "${BC_PATH}" ]]; then
    echo "  Already exists at ${BC_PATH} — skipping (delete to recreate)"
else
    conda env create -f "${ENVS_DIR}/stage3_integration.yml" -p "${BC_PATH}" -q
    echo "  Created: ${BC_PATH}"
fi

echo ""
echo "=============================================="
echo "Done!  Add the following to config/paths.local.sh"
echo "(or update config/paths.sh directly):"
echo "=============================================="
echo ""
echo "export CONDA_ENV_BASE=\"${ENV_BASE}\""
echo "export QC_ENV=\"${QC_PATH}\""
echo "export SINGLECELL_ENV=\"${SC_PATH}\""
echo "export BATCHCORR_ENV=\"${BC_PATH}\""
echo ""
