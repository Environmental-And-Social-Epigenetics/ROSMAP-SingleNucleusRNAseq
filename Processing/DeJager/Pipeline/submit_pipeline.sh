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
#   ./submit_pipeline.sh 3            # submit Stage 3 (integration, library_id)
#   ./submit_pipeline.sh 3b           # submit Stage 3b (integration, patient_id)
#   ./submit_pipeline.sh 3c           # submit Stage 3c (integration, pool_batch)
#   ./submit_pipeline.sh 3d           # submit Stage 3d (integration, derived_batch)
#   ./submit_pipeline.sh 3e           # submit Stage 3e (comprehensive comparison + report)
#   ./submit_pipeline.sh 1 2          # submit Stages 1 and 2 (2 waits for 1)
#   ./submit_pipeline.sh all          # submit all stages (3/3b/3c/3d in parallel, 3e after)
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
# Track all Stage 3 variant job IDs (they run in parallel; 3e depends on all)
STAGE3_JOB_ID=""
STAGE3B_JOB_ID=""
STAGE3C_JOB_ID=""
STAGE3D_JOB_ID=""

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
            echo "Submitting Stage 3 — Integration & annotation (library_id) ..."
            sbatch_output=$(sbatch \
                --output="${LOG_DIR}/dej_integrate_%j.out" \
                --error="${LOG_DIR}/dej_integrate_%j.err" \
                "${dep_flags[@]+"${dep_flags[@]}"}" \
                "${EXTRA_SBATCH_FLAGS[@]+"${EXTRA_SBATCH_FLAGS[@]}"}" \
                "${SCRIPT_DIR}/03_integration_annotation.sh")
            echo "  ${sbatch_output}"
            LAST_JOB_ID=$(parse_job_id "${sbatch_output}")
            STAGE3_JOB_ID="${LAST_JOB_ID}"
            ;;
        3b)
            echo "Submitting Stage 3b — Integration & annotation (patient_id) ..."
            sbatch_output=$(sbatch \
                --output="${LOG_DIR}/dej_integrate_patid_%j.out" \
                --error="${LOG_DIR}/dej_integrate_patid_%j.err" \
                "${dep_flags[@]+"${dep_flags[@]}"}" \
                "${EXTRA_SBATCH_FLAGS[@]+"${EXTRA_SBATCH_FLAGS[@]}"}" \
                "${SCRIPT_DIR}/03_integration_annotation_patient_id.sh")
            echo "  ${sbatch_output}"
            LAST_JOB_ID=$(parse_job_id "${sbatch_output}")
            STAGE3B_JOB_ID="${LAST_JOB_ID}"
            ;;
        3c)
            echo "Submitting Stage 3c — Integration & annotation (pool_batch) ..."
            sbatch_output=$(sbatch \
                --output="${LOG_DIR}/dej_integrate_pool_%j.out" \
                --error="${LOG_DIR}/dej_integrate_pool_%j.err" \
                "${dep_flags[@]+"${dep_flags[@]}"}" \
                "${EXTRA_SBATCH_FLAGS[@]+"${EXTRA_SBATCH_FLAGS[@]}"}" \
                "${SCRIPT_DIR}/03_integration_annotation_pool_batch.sh")
            echo "  ${sbatch_output}"
            LAST_JOB_ID=$(parse_job_id "${sbatch_output}")
            STAGE3C_JOB_ID="${LAST_JOB_ID}"
            ;;
        3d)
            echo "Submitting Stage 3d — Integration & annotation (derived_batch) ..."
            sbatch_output=$(sbatch \
                --output="${LOG_DIR}/dej_integrate_derived_%j.out" \
                --error="${LOG_DIR}/dej_integrate_derived_%j.err" \
                "${dep_flags[@]+"${dep_flags[@]}"}" \
                "${EXTRA_SBATCH_FLAGS[@]+"${EXTRA_SBATCH_FLAGS[@]}"}" \
                "${SCRIPT_DIR}/03_integration_annotation_derived_batch.sh")
            echo "  ${sbatch_output}"
            LAST_JOB_ID=$(parse_job_id "${sbatch_output}")
            STAGE3D_JOB_ID="${LAST_JOB_ID}"
            ;;
        3e)
            # Stage 3e depends on ALL four Stage 3 variants completing.
            local eval_dep_ids=()
            for jid in "${STAGE3_JOB_ID}" "${STAGE3B_JOB_ID}" "${STAGE3C_JOB_ID}" "${STAGE3D_JOB_ID}"; do
                [[ -n "${jid}" ]] && eval_dep_ids+=("${jid}")
            done
            local eval_dep_flags=()
            if [[ ${#eval_dep_ids[@]} -gt 0 ]]; then
                local dep_str
                dep_str=$(IFS=:; echo "${eval_dep_ids[*]}")
                eval_dep_flags=(--dependency="afterok:${dep_str}")
                echo "  (depends on jobs ${dep_str})"
            elif [[ -n "${LAST_JOB_ID}" ]]; then
                eval_dep_flags=(--dependency="afterok:${LAST_JOB_ID}")
                echo "  (depends on job ${LAST_JOB_ID})"
            fi

            echo "Submitting Stage 3e — Comprehensive batch correction comparison ..."
            sbatch_output=$(sbatch \
                --job-name=dej_compare \
                --time=04:00:00 \
                --ntasks=16 \
                --mem=256G \
                --output="${LOG_DIR}/dej_compare_%j.out" \
                --error="${LOG_DIR}/dej_compare_%j.err" \
                "${eval_dep_flags[@]+"${eval_dep_flags[@]}"}" \
                "${EXTRA_SBATCH_FLAGS[@]+"${EXTRA_SBATCH_FLAGS[@]}"}" \
                --wrap="export HDF5_USE_FILE_LOCKING=FALSE && export PYTHONUNBUFFERED=1 && source ${REPO_ROOT}/config/paths.sh && set +u && source ${CONDA_INIT_SCRIPT} && conda activate ${BATCHCORR_ENV} && set -u && export PATH=\${CONDA_PREFIX}/bin:\${PATH} && python ${SCRIPT_DIR}/03c_compare_corrections.py --input ${DEJAGER_INTEGRATED}/dejager_integrated.h5ad ${DEJAGER_INTEGRATED_PATIENT_ID}/dejager_integrated.h5ad ${DEJAGER_INTEGRATED_POOL_BATCH}/dejager_integrated.h5ad ${DEJAGER_INTEGRATED_DERIVED_BATCH}/dejager_integrated.h5ad --labels 'library_id (canonical)' 'patient_id (439)' 'pool_batch (60)' 'derived_batch (24)' --clinical-csv ${ROSMAP_CLINICAL_CSV} --output-dir ${DEJAGER_PROCESSING_OUTPUTS}/03_Evaluation")
            echo "  ${sbatch_output}"
            LAST_JOB_ID=$(parse_job_id "${sbatch_output}")
            ;;
        *)
            echo "ERROR: Unknown stage '${stage}'. Use 1, 2, 3, 3b, 3c, 3d, 3e, or 'all'."
            exit 1
            ;;
    esac
}

if [[ $# -eq 0 ]]; then
    echo "Usage: $0 [--no-preflight] <stage> [stage ...]"
    echo "  Stages: 1 (QC), 2 (doublets), 3 (integration, library_id),"
    echo "          3b (integration, patient_id), 3c (integration, pool_batch),"
    echo "          3d (integration, derived_batch), 3e (comparison + report), all"
    echo "  --no-preflight   Skip preflight checks before submission"
    exit 1
fi

for arg in "$@"; do
    if [[ "${arg}" == "all" ]]; then
        submit_stage 1
        submit_stage 2
        # All four Stage 3 variants depend on Stage 2 and run in parallel.
        # Note: all jobs load the same 122 singlet files simultaneously,
        # which may cause I/O contention. To avoid this, submit them
        # separately (e.g., ./submit_pipeline.sh 3 then 3b after it finishes).
        STAGE2_JOB_ID="${LAST_JOB_ID}"
        submit_stage 3
        LAST_JOB_ID="${STAGE2_JOB_ID}"
        submit_stage 3b
        LAST_JOB_ID="${STAGE2_JOB_ID}"
        submit_stage 3c
        LAST_JOB_ID="${STAGE2_JOB_ID}"
        submit_stage 3d
        # Stage 3e depends on all four variants (handled inside submit_stage 3e).
        submit_stage 3e
    else
        submit_stage "${arg}"
    fi
done
