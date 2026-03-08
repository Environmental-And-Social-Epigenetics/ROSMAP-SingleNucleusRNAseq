#!/bin/bash
#
# Submit all Cell Ranger scripts, then all CellBender scripts.
# No batching or scratch checks.
#

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../Config/cellranger_config.sh"

ensure_pipeline_dirs

submit_all() {
    local job_type="$1"
    local scripts=("$@")
    scripts=("${scripts[@]:1}")
    local job_ids=()

    for script in "${scripts[@]}"; do
        if [[ -f "$script" ]]; then
            job_id=$(sbatch --parsable "$script")
            job_ids+=("$job_id")
        fi
    done

    echo "${job_ids[*]}"
}

wait_for_jobs() {
    local job_type="$1"
    shift
    local job_ids=("$@")

    if [[ ${#job_ids[@]} -eq 0 ]]; then
        echo "No ${job_type} jobs to wait for."
        return
    fi

    echo "Waiting for ${#job_ids[@]} ${job_type} jobs..."
    while true; do
        local running=0
        for job_id in "${job_ids[@]}"; do
            if squeue -j "$job_id" &>/dev/null; then
                state=$(squeue -j "$job_id" -h -o "%t" 2>/dev/null || echo "")
                if [[ -n "$state" ]]; then
                    running=$((running + 1))
                fi
            fi
        done

        if [[ $running -eq 0 ]]; then
            echo "All ${job_type} jobs completed."
            break
        fi

        echo "$(date '+%H:%M:%S') - ${running} ${job_type} jobs still running..."
        sleep 300
    done
}

CR_SCRIPTS=("${BATCH_SCRIPTS_DIR}"/batch_*/cellranger/*.sh)
CB_SCRIPTS=("${BATCH_SCRIPTS_DIR}"/batch_*/cellbender/*.sh)

echo "Submitting ${#CR_SCRIPTS[@]} Cell Ranger jobs..."
cr_job_ids_str=$(submit_all "Cell Ranger" "${CR_SCRIPTS[@]}")
read -ra CR_JOB_IDS <<< "$cr_job_ids_str"

wait_for_jobs "Cell Ranger" "${CR_JOB_IDS[@]}"

echo "Submitting ${#CB_SCRIPTS[@]} CellBender jobs..."
cb_job_ids_str=$(submit_all "CellBender" "${CB_SCRIPTS[@]}")
read -ra CB_JOB_IDS <<< "$cb_job_ids_str"

wait_for_jobs "CellBender" "${CB_JOB_IDS[@]}"

echo "All Cell Ranger and CellBender jobs complete."
