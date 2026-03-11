#!/bin/bash
#
# Submit DeJager Processing pipeline stages to SLURM.
#
# This wrapper sources config/paths.sh and passes --output/--error flags on the
# sbatch command line so that SLURM logs land in the configured log directory.
# (Setting SLURM_OUTPUT/SLURM_ERROR as environment variables inside the job
# script does NOT work — SLURM reads log paths from #SBATCH directives or
# command-line flags at submission time.)
#
# When multiple stages are submitted together (e.g. "all"), each stage
# automatically depends on the previous one via --dependency=afterok.
#
# Usage:
#   ./submit_pipeline.sh 1            # submit Stage 1 (QC filtering)
#   ./submit_pipeline.sh 2            # submit Stage 2 (doublet removal)
#   ./submit_pipeline.sh 3            # submit Stage 3 (integration)
#   ./submit_pipeline.sh 1 2          # submit Stages 1 and 2 (2 waits for 1)
#   ./submit_pipeline.sh all          # submit all three stages (chained)
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

source "${REPO_ROOT}/config/paths.sh"
source "${REPO_ROOT}/config/preflight.sh"

LOG_DIR="${DEJAGER_PROCESSING_LOGS}"
mkdir -p "${LOG_DIR}"

# Build optional sbatch flags from config
EXTRA_SBATCH_FLAGS=()
if [[ -n "${SLURM_PARTITION:-}" ]]; then
    EXTRA_SBATCH_FLAGS+=(--partition="${SLURM_PARTITION}")
fi
if [[ -n "${SLURM_MAIL_USER:-}" ]]; then
    EXTRA_SBATCH_FLAGS+=(--mail-user="${SLURM_MAIL_USER}")
fi

# Parse --no-preflight flag
NO_PREFLIGHT=false
ARGS=()
for arg in "$@"; do
    if [[ "${arg}" == "--no-preflight" ]]; then
        NO_PREFLIGHT=true
    else
        ARGS+=("${arg}")
    fi
done
set -- "${ARGS[@]}"

# Track the last submitted job ID for dependency chaining
LAST_JOB_ID=""

parse_job_id() {
    # sbatch outputs: "Submitted batch job 12345"
    local output="$1"
    echo "${output}" | grep -oP '\d+$'
}

submit_stage() {
    local stage="$1"
    local dep_flags=()
    if [[ -n "${LAST_JOB_ID}" ]]; then
        dep_flags=(--dependency="afterok:${LAST_JOB_ID}")
        echo "  (depends on job ${LAST_JOB_ID})"
    fi

    # Run preflight check before submission
    if [[ "${NO_PREFLIGHT}" == "false" ]]; then
        preflight_check "dejager-stage${stage}" || \
            { echo "Preflight failed. Use --no-preflight to skip."; exit 1; }
        echo ""
    fi

    local sbatch_output
    case "${stage}" in
        1)
            echo "Submitting Stage 1 — QC filtering ..."
            sbatch_output=$(sbatch \
                --output="${LOG_DIR}/dej_qc_%A_%a.out" \
                --error="${LOG_DIR}/dej_qc_%A_%a.err" \
                "${dep_flags[@]+"${dep_flags[@]}"}" \
                "${EXTRA_SBATCH_FLAGS[@]+"${EXTRA_SBATCH_FLAGS[@]}"}" \
                "${SCRIPT_DIR}/01_qc_filter.sh")
            echo "  ${sbatch_output}"
            LAST_JOB_ID=$(parse_job_id "${sbatch_output}")

            echo "  Submitting Stage 1 aggregate plots ..."
            sbatch_output=$(sbatch \
                --dependency="afterok:${LAST_JOB_ID}" \
                --job-name=dej_qc_agg \
                --time=00:30:00 \
                --ntasks=1 \
                --mem=8G \
                --output="${LOG_DIR}/dej_qc_aggregate_%j.out" \
                --error="${LOG_DIR}/dej_qc_aggregate_%j.err" \
                "${EXTRA_SBATCH_FLAGS[@]+"${EXTRA_SBATCH_FLAGS[@]}"}" \
                --wrap="source ${REPO_ROOT}/config/paths.sh && set +u && source ${CONDA_INIT_SCRIPT} && conda activate ${QC_ENV} && set -u && python ${SCRIPT_DIR}/01_qc_filter.py --output-dir ${DEJAGER_QC_FILTERED} --aggregate-only")
            echo "  ${sbatch_output}"
            LAST_JOB_ID=$(parse_job_id "${sbatch_output}")
            ;;
        2)
            echo "Submitting Stage 2 — Doublet removal ..."
            sbatch_output=$(sbatch \
                --output="${LOG_DIR}/dej_doublets_%A_%a.out" \
                --error="${LOG_DIR}/dej_doublets_%A_%a.err" \
                "${dep_flags[@]+"${dep_flags[@]}"}" \
                "${EXTRA_SBATCH_FLAGS[@]+"${EXTRA_SBATCH_FLAGS[@]}"}" \
                "${SCRIPT_DIR}/02_doublet_removal.sh")
            echo "  ${sbatch_output}"
            LAST_JOB_ID=$(parse_job_id "${sbatch_output}")

            echo "  Submitting Stage 2 aggregate plots ..."
            sbatch_output=$(sbatch \
                --dependency="afterok:${LAST_JOB_ID}" \
                --job-name=dej_dbl_agg \
                --time=00:30:00 \
                --ntasks=1 \
                --mem=8G \
                --output="${LOG_DIR}/dej_doublets_aggregate_%j.out" \
                --error="${LOG_DIR}/dej_doublets_aggregate_%j.err" \
                "${EXTRA_SBATCH_FLAGS[@]+"${EXTRA_SBATCH_FLAGS[@]}"}" \
                --wrap="source ${REPO_ROOT}/config/paths.sh && set +u && source ${CONDA_INIT_SCRIPT} && conda activate ${SINGLECELL_ENV} && set -u && Rscript ${SCRIPT_DIR}/02_doublet_removal.Rscript --output-dir ${DEJAGER_DOUBLET_REMOVED} --input-dir ${DEJAGER_QC_FILTERED} --aggregate-only")
            echo "  ${sbatch_output}"
            LAST_JOB_ID=$(parse_job_id "${sbatch_output}")
            ;;
        3)
            echo "Submitting Stage 3 — Integration & annotation ..."
            sbatch_output=$(sbatch \
                --output="${LOG_DIR}/dej_integrate_%j.out" \
                --error="${LOG_DIR}/dej_integrate_%j.err" \
                "${dep_flags[@]+"${dep_flags[@]}"}" \
                "${EXTRA_SBATCH_FLAGS[@]+"${EXTRA_SBATCH_FLAGS[@]}"}" \
                "${SCRIPT_DIR}/03_integration_annotation.sh")
            echo "  ${sbatch_output}"
            LAST_JOB_ID=$(parse_job_id "${sbatch_output}")
            ;;
        *)
            echo "ERROR: Unknown stage '${stage}'. Use 1, 2, 3, or 'all'."
            exit 1
            ;;
    esac
}

if [[ $# -eq 0 ]]; then
    echo "Usage: $0 [--no-preflight] <stage> [stage ...]"
    echo "  Stages: 1 (QC), 2 (doublets), 3 (integration), all"
    echo "  --no-preflight   Skip preflight checks before submission"
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
