#!/bin/bash
#
# Run a single batch of Cell Ranger + CellBender jobs
#
# Usage: ./run_batch.sh <batch_number> [--cellranger-only] [--cellbender-only] [--cleanup-only]
#
# This script:
# 1. Submits all Cell Ranger jobs for the batch
# 2. Waits for Cell Ranger to complete
# 3. Submits all CellBender jobs
# 4. Waits for CellBender to complete
# 5. Cleans up scratch space
#

set -e

# Parse arguments
BATCH_NUM="${1:-}"
CELLRANGER_ONLY=false
CELLBENDER_ONLY=false
CLEANUP_ONLY=false

for arg in "$@"; do
    case $arg in
        --cellranger-only) CELLRANGER_ONLY=true ;;
        --cellbender-only) CELLBENDER_ONLY=true ;;
        --cleanup-only) CLEANUP_ONLY=true ;;
    esac
done

if [[ -z "$BATCH_NUM" ]]; then
    echo "Usage: $0 <batch_number> [--cellranger-only] [--cellbender-only] [--cleanup-only]"
    echo ""
    echo "Options:"
    echo "  --cellranger-only   Only run Cell Ranger (skip CellBender and cleanup)"
    echo "  --cellbender-only   Only run CellBender (assumes Cell Ranger completed)"
    echo "  --cleanup-only      Only run cleanup (assumes CellBender completed)"
    exit 1
fi

# Get script directory and load config
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../Config/cellranger_config.sh"

# Batch scripts directory
BATCH_DIR="${BATCH_SCRIPTS_DIR}/batch_${BATCH_NUM}"

if [[ ! -d "$BATCH_DIR" ]]; then
    echo "ERROR: Batch directory not found: $BATCH_DIR"
    echo "Run generate_batch_scripts.py first."
    exit 1
fi

echo "=============================================="
echo "Running Batch ${BATCH_NUM}"
echo "Started: $(date)"
echo "=============================================="

# Ensure directories exist
ensure_pipeline_dirs

# Function to submit jobs and collect job IDs
submit_jobs() {
    local job_type="$1"
    local scripts_dir="$2"
    local job_ids=()
    
    echo ""
    echo "Submitting ${job_type} jobs..."
    
    for script in "$scripts_dir"/*.sh; do
        if [[ -f "$script" ]]; then
            job_id=$(sbatch --parsable "$script")
            job_ids+=("$job_id")
            echo "  Submitted: $(basename "$script") -> Job $job_id"
        fi
    done
    
    echo "Submitted ${#job_ids[@]} ${job_type} jobs"
    
    # Return job IDs as space-separated string
    echo "${job_ids[*]}"
}

# Function to wait for jobs to complete
wait_for_jobs() {
    local job_type="$1"
    shift
    local job_ids=("$@")
    
    if [[ ${#job_ids[@]} -eq 0 ]]; then
        echo "No jobs to wait for"
        return
    fi
    
    echo ""
    echo "Waiting for ${#job_ids[@]} ${job_type} jobs to complete..."
    echo "Job IDs: ${job_ids[*]}"
    
    while true; do
        # Check how many jobs are still running
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
            echo "All ${job_type} jobs completed!"
            break
        fi
        
        echo "$(date '+%H:%M:%S') - $running jobs still running..."
        sleep 300  # Check every 5 minutes
    done
}

# Function to run cleanup
run_cleanup() {
    echo ""
    echo "Running cleanup for batch ${BATCH_NUM}..."
    bash "${SCRIPT_DIR}/cleanup_batch.sh" "$BATCH_NUM"
}

# Check scratch usage
check_scratch_usage

# Run Cell Ranger
if [[ "$CELLBENDER_ONLY" == "false" ]] && [[ "$CLEANUP_ONLY" == "false" ]]; then
    cr_job_ids_str=$(submit_jobs "Cell Ranger" "${BATCH_DIR}/cellranger")
    read -ra CR_JOB_IDS <<< "$cr_job_ids_str"
    
    if [[ "$CELLRANGER_ONLY" == "false" ]]; then
        wait_for_jobs "Cell Ranger" "${CR_JOB_IDS[@]}"
    else
        echo ""
        echo "Cell Ranger jobs submitted. Use 'squeue -u \$USER' to monitor."
        echo "Run with --cellbender-only after completion."
        exit 0
    fi
fi

# Run CellBender
if [[ "$CELLRANGER_ONLY" == "false" ]] && [[ "$CLEANUP_ONLY" == "false" ]]; then
    cb_job_ids_str=$(submit_jobs "CellBender" "${BATCH_DIR}/cellbender")
    read -ra CB_JOB_IDS <<< "$cb_job_ids_str"
    
    wait_for_jobs "CellBender" "${CB_JOB_IDS[@]}"
fi

# Run cleanup
if [[ "$CELLRANGER_ONLY" == "false" ]] && [[ "$CELLBENDER_ONLY" == "false" ]]; then
    run_cleanup
fi

if [[ "$CLEANUP_ONLY" == "true" ]]; then
    run_cleanup
fi

# Final status
echo ""
echo "=============================================="
echo "Batch ${BATCH_NUM} Complete!"
echo "Finished: $(date)"
echo "=============================================="

# Check tracking files
CR_COMPLETED=$(wc -l < "${TRACKING_DIR}/cellranger_completed.txt" 2>/dev/null || echo 0)
CB_COMPLETED=$(wc -l < "${TRACKING_DIR}/cellbender_completed.txt" 2>/dev/null || echo 0)
CR_FAILED=$(wc -l < "${TRACKING_DIR}/cellranger_failed.txt" 2>/dev/null || echo 0)
CB_FAILED=$(wc -l < "${TRACKING_DIR}/cellbender_failed.txt" 2>/dev/null || echo 0)

echo ""
echo "Overall Progress:"
echo "  Cell Ranger completed: $CR_COMPLETED"
echo "  CellBender completed: $CB_COMPLETED"
echo "  Cell Ranger failed: $CR_FAILED"
echo "  CellBender failed: $CB_FAILED"

check_scratch_usage

