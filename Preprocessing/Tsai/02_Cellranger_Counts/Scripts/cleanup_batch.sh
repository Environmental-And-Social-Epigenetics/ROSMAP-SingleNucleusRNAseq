#!/bin/bash
#
# Cleanup scratch space after a batch completes
#
# Usage: ./cleanup_batch.sh <batch_number> [--force]
#
# Options:
#   --force    Skip confirmation prompt for missing outputs (for automated runs)
#
# This script:
# 1. Verifies CellBender outputs are copied to permanent storage
# 2. Deletes Cell Ranger temporary files
# 3. Deletes CellBender temporary files
# 4. Reports space recovered
#

set -e

BATCH_NUM=""
FORCE_MODE=false

# Parse arguments
for arg in "$@"; do
    case $arg in
        --force)
            FORCE_MODE=true
            ;;
        *)
            if [[ -z "$BATCH_NUM" ]]; then
                BATCH_NUM="$arg"
            fi
            ;;
    esac
done

if [[ -z "$BATCH_NUM" ]]; then
    echo "Usage: $0 <batch_number> [--force]"
    exit 1
fi

# Get script directory and load config
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../Config/cellranger_config.sh"

echo "=============================================="
echo "Cleaning up Batch ${BATCH_NUM}"
echo "Started: $(date)"
echo "=============================================="

# Get patients in this batch
BATCH_ASSIGNMENTS="${TRACKING_DIR}/batch_assignments.csv"
if [[ ! -f "$BATCH_ASSIGNMENTS" ]]; then
    echo "ERROR: Batch assignments file not found: $BATCH_ASSIGNMENTS"
    exit 1
fi

# Extract projids for this batch
PROJIDS=$(tail -n +2 "$BATCH_ASSIGNMENTS" | awk -F',' -v batch="$BATCH_NUM" '$3 == batch {print $1}')
N_PATIENTS=$(echo "$PROJIDS" | wc -w)

echo "Patients in batch ${BATCH_NUM}: ${N_PATIENTS}"

# Check space before cleanup
echo ""
echo "Space usage before cleanup:"
du -sh "${CELLRANGER_OUTPUT}" 2>/dev/null || echo "  Cell Ranger scratch: (not found)"
du -sh "${CELLBENDER_OUTPUT}" 2>/dev/null || echo "  CellBender scratch: (not found)"

# Verify CellBender outputs are in permanent storage
echo ""
echo "Verifying CellBender outputs in permanent storage..."
MISSING=0
VERIFIED=0

for projid in $PROJIDS; do
    final_h5="${FINAL_OUTPUT}/${projid}/cellbender_output.h5"
    if [[ -f "$final_h5" ]]; then
        VERIFIED=$((VERIFIED + 1))
    else
        echo "  WARNING: Missing final output for ${projid}"
        MISSING=$((MISSING + 1))
    fi
done

echo "Verified: ${VERIFIED}/${N_PATIENTS}"
if [[ $MISSING -gt 0 ]]; then
    echo "Missing: ${MISSING}"
    if [[ "$FORCE_MODE" == "true" ]]; then
        echo "Force mode enabled - continuing with cleanup despite missing outputs"
    else
        echo ""
        read -p "Some outputs are missing. Continue with cleanup anyway? (y/N): " confirm
        if [[ "$confirm" != "y" ]] && [[ "$confirm" != "Y" ]]; then
            echo "Cleanup aborted."
            exit 1
        fi
    fi
fi

# Cleanup Cell Ranger outputs
echo ""
echo "Removing Cell Ranger temporary files..."
CR_REMOVED=0
for projid in $PROJIDS; do
    cr_dir="${CELLRANGER_OUTPUT}/${projid}"
    if [[ -d "$cr_dir" ]]; then
        rm -rf "$cr_dir"
        CR_REMOVED=$((CR_REMOVED + 1))
    fi
done
echo "Removed ${CR_REMOVED} Cell Ranger directories"

# Cleanup CellBender outputs
echo ""
echo "Removing CellBender temporary files..."
CB_REMOVED=0
for projid in $PROJIDS; do
    cb_dir="${CELLBENDER_OUTPUT}/${projid}"
    if [[ -d "$cb_dir" ]]; then
        rm -rf "$cb_dir"
        CB_REMOVED=$((CB_REMOVED + 1))
    fi
done
echo "Removed ${CB_REMOVED} CellBender directories"

# Check space after cleanup
echo ""
echo "Space usage after cleanup:"
du -sh "${CELLRANGER_OUTPUT}" 2>/dev/null || echo "  Cell Ranger scratch: (empty/not found)"
du -sh "${CELLBENDER_OUTPUT}" 2>/dev/null || echo "  CellBender scratch: (empty/not found)"
du -sh "${SCRATCH_ROOT}" 2>/dev/null || echo "  Total scratch: (empty/not found)"

# Record cleanup completion
echo "batch_${BATCH_NUM}: $(date)" >> "${TRACKING_DIR}/cleanup_completed.txt"

echo ""
echo "=============================================="
echo "Cleanup for Batch ${BATCH_NUM} Complete!"
echo "=============================================="

