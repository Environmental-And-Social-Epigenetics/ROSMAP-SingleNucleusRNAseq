#!/bin/bash
# =============================================================================
# ACE Epigenomic Integration — Synapse Data Download
#
# Downloads ROSMAP H3K9ac and DNA methylation tables to a configurable
# location. Requires SYNAPSE_AUTH_TOKEN (Personal Access Token) in env.
#
# Usage:
#   export SYNAPSE_AUTH_TOKEN=eyJ...    # from synapse.org settings
#   bash download_epigenomic.sh [--dry-run]
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

DRY_RUN=0
if [[ "${1:-}" == "--dry-run" ]]; then
    DRY_RUN=1
fi

# Synapse IDs are configurable via env; defaults from README.md.
# These should be verified against the current ROSMAP data dictionary.
H3K9AC_SYN_ID="${ACE_EPI_H3K9AC_SYN:-syn7428143}"
METHYL_SYN_ID="${ACE_EPI_METHYL_SYN:-syn3157275}"

# Land downloads outside the repo (large files), under WORKSPACE_ROOT.
EPI_DATA_DIR="${ACE_OUTPUT_ROOT}/EpigenomicIntegration/Tsai/epigenomic_data"
H3K9AC_DIR="${EPI_DATA_DIR}/h3k9ac"
METHYL_DIR="${EPI_DATA_DIR}/methylation"

mkdir -p "${H3K9AC_DIR}" "${METHYL_DIR}"

echo "=== ACE Epigenomic Data Download ==="
echo "H3K9ac Synapse ID:    ${H3K9AC_SYN_ID}"
echo "Methylation Syn ID:   ${METHYL_SYN_ID}"
echo "Destination root:     ${EPI_DATA_DIR}"
echo ""

if [[ -z "${SYNAPSE_AUTH_TOKEN:-}" ]]; then
    echo "ERROR: SYNAPSE_AUTH_TOKEN is not set."
    echo "       Generate a Personal Access Token at https://www.synapse.org/PersonalAccessTokens"
    echo "       and export it before running."
    exit 1
fi

if [[ "${DRY_RUN}" == "1" ]]; then
    echo "[dry-run] Would download ${H3K9AC_SYN_ID} -> ${H3K9AC_DIR}"
    echo "[dry-run] Would download ${METHYL_SYN_ID} -> ${METHYL_DIR}"
    exit 0
fi

# Pick a python with synapseclient. Prefer the synapse_env if it exists,
# otherwise fall back to consortium (which the codebase already uses widely).
SYNAPSE_PY=""
for env_path in "${SYNAPSE_ENV:-}" "${CONDA_ENV_BASE:-}/synapse_env" "${CONDA_ENV_BASE:-}/consortium"; do
    if [[ -x "${env_path}/bin/python" ]]; then
        if "${env_path}/bin/python" -c "import synapseclient" 2>/dev/null; then
            SYNAPSE_PY="${env_path}/bin/python"
            break
        fi
    fi
done

if [[ -z "${SYNAPSE_PY}" ]]; then
    echo "ERROR: No conda env with synapseclient found."
    echo "       pip install synapseclient into one of your envs and rerun."
    exit 1
fi
echo "Using synapse client at: ${SYNAPSE_PY}"

# Recursive `synapse get -r SYN_ID` materializes the full file tree under cwd.
# We change directory per Synapse target so the layout matches our expectations.
download_synid() {
    local syn_id="$1"
    local dest_dir="$2"

    echo ""
    echo ">>> Downloading ${syn_id} to ${dest_dir}"
    "${SYNAPSE_PY}" -m synapseclient -p "${SYNAPSE_AUTH_TOKEN}" \
        get -r "${syn_id}" --downloadLocation "${dest_dir}"
}

download_synid "${H3K9AC_SYN_ID}" "${H3K9AC_DIR}"
download_synid "${METHYL_SYN_ID}" "${METHYL_DIR}"

echo ""
echo "=== Download complete ==="
echo "Inventory:"
find "${EPI_DATA_DIR}" -maxdepth 3 -type f -printf '  %p (%s bytes)\n' | head -50
