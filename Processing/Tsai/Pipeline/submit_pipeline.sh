#!/bin/bash
#
# Submit Processing pipeline stages to SLURM.
#
# This wrapper sources config/paths.sh and passes --output/--error flags on the
# sbatch command line so that SLURM logs land in the configured log directory.
# (Setting SLURM_OUTPUT/SLURM_ERROR as environment variables inside the job
# script does NOT work — SLURM reads log paths from #SBATCH directives or
# command-line flags at submission time.)
#
# Usage:
#   ./submit_pipeline.sh 1            # submit Stage 1 (QC filtering)
#   ./submit_pipeline.sh 2            # submit Stage 2 (doublet removal)
#   ./submit_pipeline.sh 3            # submit Stage 3 (integration)
#   ./submit_pipeline.sh 1 2          # submit Stages 1 and 2
#   ./submit_pipeline.sh all          # submit all three stages
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

source "${REPO_ROOT}/config/paths.sh"

LOG_DIR="${TSAI_PROCESSING_LOGS}"
mkdir -p "${LOG_DIR}"

# Build optional sbatch flags from config
EXTRA_SBATCH_FLAGS=()
if [[ -n "${SLURM_PARTITION:-}" ]]; then
    EXTRA_SBATCH_FLAGS+=(--partition="${SLURM_PARTITION}")
fi
if [[ -n "${SLURM_MAIL_USER:-}" ]]; then
    EXTRA_SBATCH_FLAGS+=(--mail-user="${SLURM_MAIL_USER}")
fi

submit_stage() {
    local stage="$1"
    case "${stage}" in
        1)
            echo "Submitting Stage 1 — QC filtering ..."
            sbatch \
                --output="${LOG_DIR}/tsai_qc_%A_%a.out" \
                --error="${LOG_DIR}/tsai_qc_%A_%a.err" \
                "${EXTRA_SBATCH_FLAGS[@]+"${EXTRA_SBATCH_FLAGS[@]}"}" \
                "${SCRIPT_DIR}/01_qc_filter.sh"
            ;;
        2)
            echo "Submitting Stage 2 — Doublet removal ..."
            sbatch \
                --output="${LOG_DIR}/tsai_doublets_%A_%a.out" \
                --error="${LOG_DIR}/tsai_doublets_%A_%a.err" \
                "${EXTRA_SBATCH_FLAGS[@]+"${EXTRA_SBATCH_FLAGS[@]}"}" \
                "${SCRIPT_DIR}/02_doublet_removal.sh"
            ;;
        3)
            echo "Submitting Stage 3 — Integration & annotation ..."
            sbatch \
                --output="${LOG_DIR}/tsai_integrate_%j.out" \
                --error="${LOG_DIR}/tsai_integrate_%j.err" \
                "${EXTRA_SBATCH_FLAGS[@]+"${EXTRA_SBATCH_FLAGS[@]}"}" \
                "${SCRIPT_DIR}/03_integration_annotation.sh"
            ;;
        *)
            echo "ERROR: Unknown stage '${stage}'. Use 1, 2, 3, or 'all'."
            exit 1
            ;;
    esac
}

if [[ $# -eq 0 ]]; then
    echo "Usage: $0 <stage> [stage ...]"
    echo "  Stages: 1 (QC), 2 (doublets), 3 (integration), all"
    exit 1
fi

for arg in "$@"; do
    if [[ "${arg}" == "all" ]]; then
        submit_stage 1
        submit_stage 2
        submit_stage 3
    else
        submit_stage "${arg}"
    fi
done
