#!/bin/bash
# =============================================================================
# ACE Tsai end-to-end pipeline orchestrator.
#
# Submits every ACE workflow in dependency order. Per-workflow launchers
# remain independently usable; this is just the canonical "run everything"
# entrypoint.
#
# Usage:
#   bash run_all_ace_tsai.sh [INTEGRATION] [PHENOTYPE]
#     INTEGRATION  default: derived_batch
#     PHENOTYPE    default: all (loops over the three primary phenotypes)
#
# Optional env vars:
#   ACE_RUN_FLAGS=skip_existing  -- pass through to GSEA's skip-existing flag
#   ACE_SKIP_PREFLIGHT=1         -- skip preflight check
#   ACE_SKIP_EPI=1               -- skip EpigenomicIntegration (needs Synapse token)
#   SYNAPSE_AUTH_TOKEN=...       -- needed only if EpigenomicIntegration runs
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

INTEGRATION="${1:-derived_batch}"
PHENOTYPE="${2:-all}"

MANIFEST_DIR="${ACE_OUTPUT_ROOT}/_run_manifest"
mkdir -p "${MANIFEST_DIR}"
TIMESTAMP="$(date -u +%Y%m%dT%H%M%SZ)"
MANIFEST="${MANIFEST_DIR}/${TIMESTAMP}.tsv"
echo -e "workflow\tphenotype\tjobid\tdepends_on\tsubmit_time" > "${MANIFEST}"

log_jobs() {
    # parse `bash <launcher>` output, extract any "job NNNN" patterns and append
    # to the manifest. Argument 1: workflow name. The text comes on stdin.
    local wf="$1"
    while IFS= read -r line; do
        # Match "job 12345" or "job: 12345"
        if [[ "${line}" =~ job[:[:space:]]+([0-9]+) ]]; then
            local jobid="${BASH_REMATCH[1]}"
            local now="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
            echo -e "${wf}\t${PHENOTYPE}\t${jobid}\t${LAST_DEPS:-}\t${now}" >> "${MANIFEST}"
        fi
        echo "${line}"
    done
}

# ----------------------------------------------------------------------------
# Preflight
# ----------------------------------------------------------------------------
if [[ "${ACE_SKIP_PREFLIGHT:-0}" != "1" ]]; then
    echo "=== Preflight ==="
    bash "${REPO_ROOT}/config/preflight.sh" ace-tsai-full || {
        echo "ERROR: preflight failed. Re-run with ACE_SKIP_PREFLIGHT=1 to bypass."
        exit 1
    }
fi

echo ""
echo "=== ACE Tsai end-to-end run ==="
echo "Integration: ${INTEGRATION}"
echo "Phenotype:   ${PHENOTYPE}"
echo "Manifest:    ${MANIFEST}"
echo ""

# ----------------------------------------------------------------------------
# Stage 1: independent (DEG / CellTypeProportion). DEG produces the
# per-celltype splits and DEG result tables that downstream analyses depend
# on, so we can't formally make it parallel with downstream Stage 2. But it
# IS parallel with CellTypeProportion.
# ----------------------------------------------------------------------------
echo "--- Stage 1: DEG + CellTypeProportion (independent) ---"
LAST_DEPS=""

bash "${SCRIPT_DIR}/DEG/Tsai/aceDegT.sh" "${INTEGRATION}" 2>&1 | log_jobs "DEG"
bash "${SCRIPT_DIR}/CellTypeProportion/Tsai/acePropT.sh" "${INTEGRATION}" 2>&1 | log_jobs "CellTypeProportion"

# We don't currently track DEG job IDs deterministically, so downstream
# stages here just submit without an afterok dependency. That matches the
# current usage pattern where DEG is normally already complete before
# downstream jobs are submitted.
# ----------------------------------------------------------------------------
# Stage 2: workflows that read DEG results
# ----------------------------------------------------------------------------
echo ""
echo "--- Stage 2: GSEA / TFActivity / SCENIC / MicState / hdWGCNA / CellChat ---"

ACE_GSEA_SKIP_EXISTING="${ACE_GSEA_SKIP_EXISTING:-1}" \
    bash "${SCRIPT_DIR}/GSEA/Tsai/aceGseaT.sh" "${INTEGRATION}" 2>&1 | log_jobs "GSEA"

bash "${SCRIPT_DIR}/TFActivity/Tsai/aceTfActT.sh" "${INTEGRATION}" 2>&1 | log_jobs "TFActivity"

if [[ "${PHENOTYPE}" == "all" ]]; then
    PHENOTYPES=("tot_adverse_exp" "early_hh_ses" "ace_aggregate")
else
    PHENOTYPES=("${PHENOTYPE}")
fi

for ph in "${PHENOTYPES[@]}"; do
    bash "${SCRIPT_DIR}/SCENIC/Tsai/aceScenicT.sh" "${INTEGRATION}" "${ph}" 2>&1 | log_jobs "SCENIC"
done

bash "${SCRIPT_DIR}/MicState/Tsai/aceMicStateT.sh" "${INTEGRATION}" "${PHENOTYPE}" 2>&1 | log_jobs "MicState"
bash "${SCRIPT_DIR}/hdWGCNA/Tsai/aceWgcnaT.sh"   "${INTEGRATION}" "${PHENOTYPE}" 2>&1 | log_jobs "hdWGCNA"
bash "${SCRIPT_DIR}/CellChat/Tsai/aceCellChatT.sh" "${INTEGRATION}" "${PHENOTYPE}" 2>&1 | log_jobs "CellChat"

# ----------------------------------------------------------------------------
# Stage 3: epigenomic integration (optional)
# ----------------------------------------------------------------------------
if [[ "${ACE_SKIP_EPI:-0}" != "1" ]]; then
    echo ""
    echo "--- Stage 3: EpigenomicIntegration ---"
    if [[ -z "${SYNAPSE_AUTH_TOKEN:-}" && ! -d "${ACE_OUTPUT_ROOT}/EpigenomicIntegration/Tsai/epigenomic_data" ]]; then
        echo "Skipping EpigenomicIntegration: SYNAPSE_AUTH_TOKEN not set and data not yet downloaded."
        echo "  Run with SYNAPSE_AUTH_TOKEN=... or first 'bash EpigenomicIntegration/Tsai/download_epigenomic.sh'"
    else
        bash "${SCRIPT_DIR}/EpigenomicIntegration/Tsai/aceEpiIntT.sh" "${INTEGRATION}" "${PHENOTYPE}" 2>&1 | log_jobs "EpigenomicIntegration"
    fi
fi

echo ""
echo "=== Submission complete ==="
echo "Manifest written to ${MANIFEST}"
echo ""
echo "Track: squeue -u \$USER"
echo "Verify when complete: python ${SCRIPT_DIR}/_shared/verify_outputs.py"
