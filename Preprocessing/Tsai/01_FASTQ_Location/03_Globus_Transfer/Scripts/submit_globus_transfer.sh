#!/bin/bash
#
# Submit Globus batch transfer for Tsai FASTQ files
#
# This script activates the Globus CLI conda environment and submits
# the batch transfer from Engaging to Openmind.
#
# Usage:
#   ./submit_globus_transfer.sh [--check-only]
#
# Options:
#   --check-only    Only check login status and validate batch file, don't submit
#

set -e

# =============================================================================
# Configuration
# =============================================================================

# Globus endpoint IDs
# Note: Use MIT ORCD Engaging Collection, not the older mithpc#engaging
ENGAGING_ENDPOINT="ec54b570-cac5-47f7-b2a1-100c2078686f"  # MIT ORCD Engaging Collection
OPENMIND_ENDPOINT="cbc6f8da-d37e-11eb-bde9-5111456017d9"

# Conda environment with Globus CLI
GLOBUS_CONDA_ENV="/home/mabdel03/conda_envs/globus_env"

# Script directory (resolve relative paths)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BASE_DIR="$(dirname "$SCRIPT_DIR")"

# Batch file location
BATCH_FILE="${BASE_DIR}/Batch_Files/globus_batch.txt"
LOG_DIR="${BASE_DIR}/Logs"

# Transfer label
TRANSFER_LABEL="Tsai_FASTQ_Transfer_$(date +%Y%m%d_%H%M%S)"

# =============================================================================
# Functions
# =============================================================================

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2
}

# =============================================================================
# Main Script
# =============================================================================

log "Globus Transfer Submission Script"
log "=================================="
echo

# Parse arguments
CHECK_ONLY=false
if [[ "$1" == "--check-only" ]]; then
    CHECK_ONLY=true
    log "Running in check-only mode"
fi

# Check if batch file exists
if [[ ! -f "$BATCH_FILE" ]]; then
    error "Batch file not found: $BATCH_FILE"
    error "Please run generate_globus_batch.py first"
    exit 1
fi

# Count entries in batch file
BATCH_LINES=$(wc -l < "$BATCH_FILE")
log "Batch file: $BATCH_FILE"
log "Entries to transfer: $BATCH_LINES"
echo

# Activate conda environment
log "Activating Globus conda environment..."
if [[ -f "/home/mabdel03/miniforge3/etc/profile.d/conda.sh" ]]; then
    source /home/mabdel03/miniforge3/etc/profile.d/conda.sh
elif [[ -f "$HOME/.bashrc" ]]; then
    source "$HOME/.bashrc"
fi

conda activate "$GLOBUS_CONDA_ENV"
log "Conda environment activated: $GLOBUS_CONDA_ENV"
echo

# Check Globus login status
log "Checking Globus login status..."
if globus whoami &>/dev/null; then
    GLOBUS_USER=$(globus whoami)
    log "Logged in as: $GLOBUS_USER"
else
    log "Not logged in. Initiating Globus login..."
    globus login
    GLOBUS_USER=$(globus whoami)
    log "Now logged in as: $GLOBUS_USER"
fi
echo

# Verify endpoint access
log "Verifying endpoint access..."
log "  Source (Engaging): $ENGAGING_ENDPOINT"
log "  Destination (Openmind): $OPENMIND_ENDPOINT"

# Check if endpoints are accessible
if ! globus endpoint show "$ENGAGING_ENDPOINT" &>/dev/null; then
    error "Cannot access source endpoint: $ENGAGING_ENDPOINT"
    error "You may need to activate the endpoint or check permissions"
    exit 1
fi

if ! globus endpoint show "$OPENMIND_ENDPOINT" &>/dev/null; then
    error "Cannot access destination endpoint: $OPENMIND_ENDPOINT"
    error "You may need to activate the endpoint or check permissions"
    exit 1
fi

log "Both endpoints are accessible"
echo

# Show sample entries from batch file
log "Sample batch entries (first 3):"
head -3 "$BATCH_FILE" | while read -r line; do
    echo "  $line"
done
echo

if [[ "$CHECK_ONLY" == true ]]; then
    log "Check-only mode: Skipping transfer submission"
    log "To submit the transfer, run without --check-only"
    exit 0
fi

# Confirm submission
echo "=========================================="
echo "Ready to submit Globus transfer:"
echo "  From: $ENGAGING_ENDPOINT (Engaging)"
echo "  To:   $OPENMIND_ENDPOINT (Openmind)"
echo "  Files: $BATCH_LINES"
echo "  Label: $TRANSFER_LABEL"
echo "=========================================="
echo

read -p "Submit transfer? [y/N] " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    log "Transfer cancelled by user"
    exit 0
fi

# Submit the transfer
log "Submitting Globus batch transfer..."
mkdir -p "$LOG_DIR"

TRANSFER_OUTPUT=$(globus transfer "$ENGAGING_ENDPOINT" "$OPENMIND_ENDPOINT" \
    --batch "$BATCH_FILE" \
    --label "$TRANSFER_LABEL" \
    --notify on \
    2>&1)

echo "$TRANSFER_OUTPUT"

# Extract task ID from output
TASK_ID=$(echo "$TRANSFER_OUTPUT" | grep -oP 'Task ID: \K[a-f0-9-]+' || true)

if [[ -n "$TASK_ID" ]]; then
    log "Transfer submitted successfully!"
    log "Task ID: $TASK_ID"
    
    # Save task ID to log file
    TASK_LOG="${LOG_DIR}/transfer_task_${TRANSFER_LABEL}.log"
    {
        echo "Globus Transfer Log"
        echo "==================="
        echo "Timestamp: $(date)"
        echo "Task ID: $TASK_ID"
        echo "Label: $TRANSFER_LABEL"
        echo "Source: $ENGAGING_ENDPOINT"
        echo "Destination: $OPENMIND_ENDPOINT"
        echo "Batch file: $BATCH_FILE"
        echo "Entries: $BATCH_LINES"
        echo ""
        echo "To check status:"
        echo "  globus task show $TASK_ID"
        echo ""
        echo "To monitor:"
        echo "  globus task wait $TASK_ID"
        echo ""
        echo "Web interface:"
        echo "  https://app.globus.org/activity/$TASK_ID"
    } > "$TASK_LOG"
    
    log "Task log saved: $TASK_LOG"
    echo
    log "Monitor transfer at: https://app.globus.org/activity/$TASK_ID"
    log "Or use: globus task show $TASK_ID"
else
    error "Could not extract Task ID from output"
    error "Check output above for details"
    exit 1
fi

echo
log "Transfer submission complete!"
