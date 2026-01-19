#!/bin/bash
#
# Retry failed Cell Ranger and CellBender jobs
#
# Usage: ./retry_failed.sh [--cellranger-only] [--cellbender-only] [--dry-run]
#
# This script:
# 1. Reads failed patient IDs from tracking files
# 2. Finds and resubmits their SLURM scripts
# 3. Reports what was resubmitted
#

set -e

# =============================================================================
# PARSE ARGUMENTS
# =============================================================================

RETRY_CR=true
RETRY_CB=true
DRY_RUN=false

for arg in "$@"; do
    case $arg in
        --cellranger-only)
            RETRY_CR=true
            RETRY_CB=false
            ;;
        --cellbender-only)
            RETRY_CR=false
            RETRY_CB=true
            ;;
        --dry-run)
            DRY_RUN=true
            ;;
        *)
            echo "Unknown option: $arg"
            echo "Usage: $0 [--cellranger-only] [--cellbender-only] [--dry-run]"
            exit 1
            ;;
    esac
done

# =============================================================================
# SETUP
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../Config/cellranger_config.sh"

echo "=============================================="
echo "Retrying Failed Jobs"
echo "Started: $(date)"
echo "=============================================="

if [[ "$DRY_RUN" == "true" ]]; then
    echo "[DRY RUN MODE - No jobs will be submitted]"
fi

# =============================================================================
# FIND FAILED JOBS
# =============================================================================

CR_FAILED_FILE="${TRACKING_DIR}/cellranger_failed.txt"
CB_FAILED_FILE="${TRACKING_DIR}/cellbender_failed.txt"
BATCH_ASSIGNMENTS="${TRACKING_DIR}/batch_assignments.csv"

# Function to find batch number for a projid
get_batch_for_projid() {
    local projid="$1"
    grep "^${projid}," "$BATCH_ASSIGNMENTS" | cut -d',' -f3
}

# =============================================================================
# RETRY CELL RANGER
# =============================================================================

if [[ "$RETRY_CR" == "true" ]]; then
    echo ""
    echo "=== Cell Ranger Failed Jobs ==="
    
    if [[ -f "$CR_FAILED_FILE" ]] && [[ -s "$CR_FAILED_FILE" ]]; then
        # Get unique projids
        CR_FAILED=$(sort -u "$CR_FAILED_FILE")
        N_CR_FAILED=$(echo "$CR_FAILED" | wc -l)
        echo "Found ${N_CR_FAILED} failed Cell Ranger jobs"
        
        CR_SUBMITTED=0
        for projid in $CR_FAILED; do
            batch_num=$(get_batch_for_projid "$projid")
            script="${BATCH_SCRIPTS_DIR}/batch_${batch_num}/cellranger/${projid}_cellranger.sh"
            
            if [[ -f "$script" ]]; then
                if [[ "$DRY_RUN" == "true" ]]; then
                    echo "  [DRY RUN] Would submit: $script"
                else
                    job_id=$(sbatch --parsable "$script")
                    echo "  Submitted ${projid} -> Job ${job_id}"
                    CR_SUBMITTED=$((CR_SUBMITTED + 1))
                fi
            else
                echo "  WARNING: Script not found for ${projid}: ${script}"
            fi
        done
        
        if [[ "$DRY_RUN" == "false" ]]; then
            echo "Submitted ${CR_SUBMITTED} Cell Ranger retry jobs"
            
            # Clear the failed file (they'll be re-added if they fail again)
            > "$CR_FAILED_FILE"
            echo "Cleared ${CR_FAILED_FILE}"
        fi
    else
        echo "No failed Cell Ranger jobs found"
    fi
fi

# =============================================================================
# RETRY CELLBENDER
# =============================================================================

if [[ "$RETRY_CB" == "true" ]]; then
    echo ""
    echo "=== CellBender Failed Jobs ==="
    
    if [[ -f "$CB_FAILED_FILE" ]] && [[ -s "$CB_FAILED_FILE" ]]; then
        # Get unique projids
        CB_FAILED=$(sort -u "$CB_FAILED_FILE")
        N_CB_FAILED=$(echo "$CB_FAILED" | wc -l)
        echo "Found ${N_CB_FAILED} failed CellBender jobs"
        
        CB_SUBMITTED=0
        for projid in $CB_FAILED; do
            batch_num=$(get_batch_for_projid "$projid")
            script="${BATCH_SCRIPTS_DIR}/batch_${batch_num}/cellbender/${projid}_cellbender.sh"
            
            if [[ -f "$script" ]]; then
                if [[ "$DRY_RUN" == "true" ]]; then
                    echo "  [DRY RUN] Would submit: $script"
                else
                    job_id=$(sbatch --parsable "$script")
                    echo "  Submitted ${projid} -> Job ${job_id}"
                    CB_SUBMITTED=$((CB_SUBMITTED + 1))
                fi
            else
                echo "  WARNING: Script not found for ${projid}: ${script}"
            fi
        done
        
        if [[ "$DRY_RUN" == "false" ]]; then
            echo "Submitted ${CB_SUBMITTED} CellBender retry jobs"
            
            # Clear the failed file
            > "$CB_FAILED_FILE"
            echo "Cleared ${CB_FAILED_FILE}"
        fi
    else
        echo "No failed CellBender jobs found"
    fi
fi

# =============================================================================
# SUMMARY
# =============================================================================

echo ""
echo "=============================================="
echo "Retry Complete"
echo "=============================================="

if [[ "$DRY_RUN" == "true" ]]; then
    echo "[DRY RUN] No jobs were actually submitted"
else
    echo "Monitor with: squeue -u \$USER"
fi
