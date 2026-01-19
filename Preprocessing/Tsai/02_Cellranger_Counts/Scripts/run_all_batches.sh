#!/bin/bash
#
# Master orchestrator: Run Cell Ranger + CellBender for all batches
#
# Usage:
#   ./run_all_batches.sh [--start-batch N] [--end-batch M] [--skip-cellranger-batch N]
#
# Options:
#   --start-batch N        Start from batch N (default: 1)
#   --end-batch M          End at batch M (default: 16)
#   --skip-cellranger-batch N  Skip Cell Ranger for batch N (already complete)
#
# Example:
#   ./run_all_batches.sh --skip-cellranger-batch 1   # Start from batch 1, skip CR for batch 1
#

set -e

# =============================================================================
# PARSE ARGUMENTS
# =============================================================================

START_BATCH=1
END_BATCH=16
SKIP_CR_BATCHES=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --start-batch)
            START_BATCH="$2"
            shift 2
            ;;
        --end-batch)
            END_BATCH="$2"
            shift 2
            ;;
        --skip-cellranger-batch)
            SKIP_CR_BATCHES="${SKIP_CR_BATCHES} $2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# =============================================================================
# SETUP
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../Config/cellranger_config.sh"

# Create timestamped log file
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
MASTER_LOG="${LOGS_DIR}/pipeline_master_${TIMESTAMP}.log"

# Ensure directories exist
ensure_pipeline_dirs

# Logging function
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$MASTER_LOG"
}

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

submit_jobs() {
    local job_type="$1"
    local scripts_dir="$2"
    local job_ids=()
    
    for script in "$scripts_dir"/*.sh; do
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
        return
    fi
    
    log "Waiting for ${#job_ids[@]} ${job_type} jobs..."
    
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
            log "All ${job_type} jobs completed"
            break
        fi
        
        log "  $running jobs still running..."
        sleep 300  # Check every 5 minutes
    done
}

run_cleanup() {
    local batch_num="$1"
    log "Running cleanup for batch ${batch_num}..."
    bash "${SCRIPT_DIR}/cleanup_batch.sh" "$batch_num" --force
}

should_skip_cellranger() {
    local batch="$1"
    for skip_batch in $SKIP_CR_BATCHES; do
        if [[ "$batch" == "$skip_batch" ]]; then
            return 0
        fi
    done
    return 1
}

# =============================================================================
# MAIN PIPELINE
# =============================================================================

log "=============================================="
log "Starting Full Pipeline"
log "=============================================="
log "Start batch: ${START_BATCH}"
log "End batch: ${END_BATCH}"
log "Skip Cell Ranger for batches: ${SKIP_CR_BATCHES:-none}"
log "Master log: ${MASTER_LOG}"
log ""

check_scratch_usage

for BATCH_NUM in $(seq "$START_BATCH" "$END_BATCH"); do
    log ""
    log "=============================================="
    log "BATCH ${BATCH_NUM} of ${END_BATCH}"
    log "Started: $(date)"
    log "=============================================="
    
    BATCH_DIR="${BATCH_SCRIPTS_DIR}/batch_${BATCH_NUM}"
    
    if [[ ! -d "$BATCH_DIR" ]]; then
        log "ERROR: Batch directory not found: $BATCH_DIR"
        continue
    fi
    
    # -------------------------------------------
    # CELL RANGER
    # -------------------------------------------
    if should_skip_cellranger "$BATCH_NUM"; then
        log "Skipping Cell Ranger for batch ${BATCH_NUM} (already complete)"
    else
        log "Submitting Cell Ranger jobs..."
        cr_job_ids_str=$(submit_jobs "Cell Ranger" "${BATCH_DIR}/cellranger")
        read -ra CR_JOB_IDS <<< "$cr_job_ids_str"
        log "Submitted ${#CR_JOB_IDS[@]} Cell Ranger jobs"
        
        wait_for_jobs "Cell Ranger" "${CR_JOB_IDS[@]}"
    fi
    
    # Check scratch before CellBender
    check_scratch_usage
    
    # -------------------------------------------
    # CELLBENDER
    # -------------------------------------------
    log "Submitting CellBender jobs..."
    cb_job_ids_str=$(submit_jobs "CellBender" "${BATCH_DIR}/cellbender")
    read -ra CB_JOB_IDS <<< "$cb_job_ids_str"
    log "Submitted ${#CB_JOB_IDS[@]} CellBender jobs"
    
    wait_for_jobs "CellBender" "${CB_JOB_IDS[@]}"
    
    # -------------------------------------------
    # CLEANUP
    # -------------------------------------------
    run_cleanup "$BATCH_NUM"
    
    # Report progress
    log ""
    log "Batch ${BATCH_NUM} complete!"
    CR_DONE=$(wc -l < "${TRACKING_DIR}/cellranger_completed.txt" 2>/dev/null || echo 0)
    CB_DONE=$(wc -l < "${TRACKING_DIR}/cellbender_completed.txt" 2>/dev/null || echo 0)
    CR_FAIL=$(wc -l < "${TRACKING_DIR}/cellranger_failed.txt" 2>/dev/null || echo 0)
    CB_FAIL=$(wc -l < "${TRACKING_DIR}/cellbender_failed.txt" 2>/dev/null || echo 0)
    log "Progress: CR completed=${CR_DONE}, CB completed=${CB_DONE}, CR failed=${CR_FAIL}, CB failed=${CB_FAIL}"
    check_scratch_usage
done

# =============================================================================
# FINAL SUMMARY
# =============================================================================

log ""
log "=============================================="
log "PIPELINE COMPLETE"
log "Finished: $(date)"
log "=============================================="

CR_DONE=$(wc -l < "${TRACKING_DIR}/cellranger_completed.txt" 2>/dev/null || echo 0)
CB_DONE=$(wc -l < "${TRACKING_DIR}/cellbender_completed.txt" 2>/dev/null || echo 0)
CR_FAIL=$(wc -l < "${TRACKING_DIR}/cellranger_failed.txt" 2>/dev/null || echo 0)
CB_FAIL=$(wc -l < "${TRACKING_DIR}/cellbender_failed.txt" 2>/dev/null || echo 0)

log ""
log "Final Statistics:"
log "  Cell Ranger completed: ${CR_DONE}"
log "  CellBender completed: ${CB_DONE}"
log "  Cell Ranger failed: ${CR_FAIL}"
log "  CellBender failed: ${CB_FAIL}"
log ""
log "Master log saved to: ${MASTER_LOG}"

if [[ "$CR_FAIL" -gt 0 ]] || [[ "$CB_FAIL" -gt 0 ]]; then
    log ""
    log "ATTENTION: Some jobs failed. Run retry_failed.sh to reprocess them."
fi
