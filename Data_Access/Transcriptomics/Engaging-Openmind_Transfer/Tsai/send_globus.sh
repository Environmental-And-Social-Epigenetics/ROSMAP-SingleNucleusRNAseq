#!/bin/bash
#
# Send Tsai transcriptomics data from the local cluster to the remote cluster
# via Globus
#
# Auto-detects whether this is running on Openmind or Engaging and sets the
# source (local) and destination (remote) endpoints accordingly.
#
# Usage:
#   ./send_globus.sh [OPTIONS]
#
# Options:
#   --batch FILE     Use a pre-generated batch file instead of generating one
#   --check-only     Verify login/endpoints without submitting a transfer
#   --data-type TYPE Filter to a specific data type (FASTQs, Cellranger_Output,
#                    Cellbender_Output) when generating a batch file on the fly
#

set -euo pipefail

# =============================================================================
# Source project paths
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

# =============================================================================
# Globus endpoint IDs
# =============================================================================

OPENMIND_ENDPOINT="cbc6f8da-d37e-11eb-bde9-5111456017d9"
ENGAGING_ENDPOINT="c52fcff2-761c-11eb-8cfc-cd623f92e1c0"

# =============================================================================
# Auto-detect cluster
# =============================================================================

if [[ -d "/om2" ]]; then
    CURRENT_CLUSTER="openmind"
    CONDA_INIT="/om2/user/mabdel03/anaconda/etc/profile.d/conda.sh"
    GLOBUS_ENV="/om2/user/mabdel03/conda_envs/globus_env"
    SOURCE_ENDPOINT="$OPENMIND_ENDPOINT"
    DEST_ENDPOINT="$ENGAGING_ENDPOINT"
    SOURCE_LABEL="Openmind"
    DEST_LABEL="Engaging"
else
    CURRENT_CLUSTER="engaging"
    CONDA_INIT="/home/mabdel03/miniforge3/etc/profile.d/conda.sh"
    GLOBUS_ENV="/home/mabdel03/conda_envs/globus_env"
    SOURCE_ENDPOINT="$ENGAGING_ENDPOINT"
    DEST_ENDPOINT="$OPENMIND_ENDPOINT"
    SOURCE_LABEL="Engaging"
    DEST_LABEL="Openmind"
fi

# =============================================================================
# Configuration
# =============================================================================

COHORT="Tsai"
LOG_DIR="${SCRIPT_DIR}/Logs"
TRANSFER_LABEL="${COHORT}_Send_${SOURCE_LABEL}_to_${DEST_LABEL}_$(date +%Y%m%d_%H%M%S)"

# =============================================================================
# Functions
# =============================================================================

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2
}

usage() {
    echo "Usage: $0 [--batch FILE] [--check-only] [--data-type TYPE]"
    echo ""
    echo "Options:"
    echo "  --batch FILE       Use a pre-generated batch file"
    echo "  --check-only       Verify login and endpoints without submitting"
    echo "  --data-type TYPE   Filter data type (FASTQs, Cellranger_Output, Cellbender_Output)"
}

# =============================================================================
# Parse arguments
# =============================================================================

BATCH_FILE=""
CHECK_ONLY=false
DATA_TYPE="all"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --batch)
            BATCH_FILE="$2"
            shift 2
            ;;
        --check-only)
            CHECK_ONLY=true
            shift
            ;;
        --data-type)
            DATA_TYPE="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            error "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# =============================================================================
# Main
# =============================================================================

log "Globus Transfer: SEND ${COHORT} Data"
log "======================================"
log "Cluster detected: ${CURRENT_CLUSTER}"
log "Direction: ${SOURCE_LABEL} -> ${DEST_LABEL}"
echo

# ---- Activate Globus conda environment ----
log "Activating Globus conda environment..."
source "$CONDA_INIT"
conda activate "$GLOBUS_ENV"
log "Conda environment activated: ${GLOBUS_ENV}"
echo

# ---- Check Globus login status ----
log "Checking Globus login status..."
if globus whoami &>/dev/null; then
    GLOBUS_USER=$(globus whoami)
    log "Logged in as: ${GLOBUS_USER}"
else
    log "Not logged in. Initiating Globus login..."
    globus login
    GLOBUS_USER=$(globus whoami)
    log "Now logged in as: ${GLOBUS_USER}"
fi
echo

# ---- Verify endpoints ----
log "Verifying endpoint access..."
log "  Source (${SOURCE_LABEL}): ${SOURCE_ENDPOINT}"
log "  Destination (${DEST_LABEL}): ${DEST_ENDPOINT}"

if ! globus endpoint show "$SOURCE_ENDPOINT" &>/dev/null; then
    error "Cannot access source endpoint: ${SOURCE_ENDPOINT}"
    error "You may need to activate the endpoint or check permissions."
    exit 1
fi

if ! globus endpoint show "$DEST_ENDPOINT" &>/dev/null; then
    error "Cannot access destination endpoint: ${DEST_ENDPOINT}"
    error "You may need to activate the endpoint or check permissions."
    exit 1
fi

log "Both endpoints are accessible."
echo

# ---- Prepare batch file ----
if [[ -z "$BATCH_FILE" ]]; then
    log "No batch file specified; generating one on the fly..."
    BATCH_FILE=$("${SCRIPT_DIR}/generate_batch.sh" "$DATA_TYPE" | tail -1 | grep -oP 'Batch file written: \K.*' || true)

    # Fallback: use the latest symlink
    if [[ -z "$BATCH_FILE" || ! -f "$BATCH_FILE" ]]; then
        BATCH_FILE="${SCRIPT_DIR}/Batch_Files/tsai_latest_batch.txt"
    fi
fi

if [[ ! -f "$BATCH_FILE" ]]; then
    error "Batch file not found: ${BATCH_FILE}"
    exit 1
fi

BATCH_LINES=$(wc -l < "$BATCH_FILE")
log "Batch file: ${BATCH_FILE}"
log "Entries to transfer: ${BATCH_LINES}"
echo

# ---- Show sample entries ----
log "Sample batch entries (first 3):"
head -3 "$BATCH_FILE" | while IFS= read -r line; do
    echo "  ${line}"
done
echo

# ---- Check-only mode ----
if [[ "$CHECK_ONLY" == true ]]; then
    log "Check-only mode: skipping transfer submission."
    log "To submit, run without --check-only."
    exit 0
fi

# ---- Confirm submission ----
echo "=========================================="
echo "Ready to submit Globus transfer:"
echo "  From: ${SOURCE_ENDPOINT} (${SOURCE_LABEL})"
echo "  To:   ${DEST_ENDPOINT} (${DEST_LABEL})"
echo "  Entries: ${BATCH_LINES}"
echo "  Label: ${TRANSFER_LABEL}"
echo "=========================================="
echo

read -p "Submit transfer? [y/N] " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    log "Transfer cancelled by user."
    exit 0
fi

# ---- Submit transfer ----
log "Submitting Globus batch transfer..."
mkdir -p "$LOG_DIR"

TRANSFER_OUTPUT=$(globus transfer "$SOURCE_ENDPOINT" "$DEST_ENDPOINT" \
    --batch "$BATCH_FILE" \
    --label "$TRANSFER_LABEL" \
    --notify on \
    2>&1)

echo "$TRANSFER_OUTPUT"

# ---- Extract and save task ID ----
TASK_ID=$(echo "$TRANSFER_OUTPUT" | grep -oP 'Task ID: \K[a-f0-9-]+' || true)

if [[ -n "$TASK_ID" ]]; then
    log "Transfer submitted successfully!"
    log "Task ID: ${TASK_ID}"

    TASK_LOG="${LOG_DIR}/send_task_${TRANSFER_LABEL}.log"
    {
        echo "Globus Transfer Log"
        echo "==================="
        echo "Timestamp: $(date)"
        echo "Direction: SEND (${SOURCE_LABEL} -> ${DEST_LABEL})"
        echo "Task ID: ${TASK_ID}"
        echo "Label: ${TRANSFER_LABEL}"
        echo "Source: ${SOURCE_ENDPOINT} (${SOURCE_LABEL})"
        echo "Destination: ${DEST_ENDPOINT} (${DEST_LABEL})"
        echo "Batch file: ${BATCH_FILE}"
        echo "Entries: ${BATCH_LINES}"
        echo ""
        echo "To check status:"
        echo "  globus task show ${TASK_ID}"
        echo ""
        echo "To monitor:"
        echo "  globus task wait ${TASK_ID}"
        echo ""
        echo "Web interface:"
        echo "  https://app.globus.org/activity/${TASK_ID}"
    } > "$TASK_LOG"

    log "Task log saved: ${TASK_LOG}"
    echo
    log "Monitor transfer at: https://app.globus.org/activity/${TASK_ID}"
    log "Or use: globus task show ${TASK_ID}"
else
    error "Could not extract Task ID from output."
    error "Check the output above for details."
    exit 1
fi

echo
log "Transfer submission complete!"
