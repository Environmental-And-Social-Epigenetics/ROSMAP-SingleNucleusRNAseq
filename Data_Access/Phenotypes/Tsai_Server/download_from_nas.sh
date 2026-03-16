#!/bin/bash
# =============================================================================
# download_from_nas.sh — Download Phenotype data from Tsai Lab NAS
# =============================================================================
# Usage:
#   ./download_from_nas.sh              Download Phenotypes from NAS
#   ./download_from_nas.sh verify       Verify local files against NAS
#
# Phenotypes is a single flat directory of CSV files — no sharding needed.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---------------------------------------------------------------------------
# Source configuration
# ---------------------------------------------------------------------------
source "$(cd "${SCRIPT_DIR}/../../../config" && pwd)/paths.sh"

# ---------------------------------------------------------------------------
# Connection settings
# ---------------------------------------------------------------------------
REMOTE_HOST="tsailabnas.mit.edu"
SFTP_PORT=22
SFTP_USER="mabdel03"
SFTP_AUTHFILE="${HOME}/.smb_tsailabnas"
SFTP_NUM_REQUESTS=64
SFTP_CIPHER="aes128-gcm@openssh.com"
SFTP_BUFFER_SIZE=""

LOCAL_BASE="${REPO_ROOT}/Data/Phenotypes"
SFTP_REMOTE_BASE="/LabShared/mabdel03/ROSMAP/Data/Phenotypes"
STATE_DIR="${HOME}/.nas_transfer_state"

# ---------------------------------------------------------------------------
# Resolve mode
# ---------------------------------------------------------------------------
MODE="${1:-download}"
if [[ "$MODE" != "download" && "$MODE" != "verify" ]]; then
    echo "ERROR: Unknown mode '$MODE'. Use: download | verify"
    exit 1
fi

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

get_password() {
    local pass
    pass=$(grep -i '^password' "$SFTP_AUTHFILE" | sed 's/^[^=]*= *//')
    if [[ -z "$pass" ]]; then
        log "ERROR: Could not extract password from $SFTP_AUTHFILE"
        exit 1
    fi
    echo "$pass"
}

sftp_cmd() {
    local pass="$1"
    local extra_args=("${@:2}")
    local cmd=(sshpass -p "$pass" sftp
        -oStrictHostKeyChecking=no
        -oBatchMode=no
        -oPort="${SFTP_PORT}"
        -R "${SFTP_NUM_REQUESTS}"
        -c "${SFTP_CIPHER}")
    if [[ -n "${SFTP_BUFFER_SIZE}" ]]; then
        cmd+=(-B "${SFTP_BUFFER_SIZE}")
    fi
    cmd+=("${extra_args[@]}" "${SFTP_USER}@${REMOTE_HOST}")
    echo "${cmd[@]}"
}

# ---------------------------------------------------------------------------
# Build local manifest: relative_path<TAB>size_bytes
# ---------------------------------------------------------------------------
build_local_manifest() {
    local dir="$1"
    (cd "$dir" && find . -type f -exec stat --format='%n\t%s' {} + | sort)
}

# ---------------------------------------------------------------------------
# Build remote manifest (handles WIN\Domain Users group name with space)
# ---------------------------------------------------------------------------
build_remote_manifest() {
    local pass="$1"
    local remote_dir="$2"
    local tmpfile
    tmpfile=$(mktemp)

    local batchfile
    batchfile=$(mktemp)
    echo "ls -lR ${remote_dir}" > "$batchfile"

    eval "$(sftp_cmd "$pass" -b "$batchfile")" > "$tmpfile" 2>/dev/null || true
    rm -f "$batchfile"

    awk -v base="$remote_dir" '
    /^sftp>/ { next }
    /^Listing/ { next }
    /^$/ { next }
    /^\// {
        sub(/:$/, "")
        current_dir = $0
        sub("^" base "/?", "./", current_dir)
        next
    }
    /^d/ { next }
    /^total/ { next }
    /^-/ {
        size = ""
        for (i = 2; i <= NF; i++) {
            if ($i ~ /^[0-9]+$/) {
                size = $i
                break
            }
        }
        if (size == "") next
        fname = $NF
        if (current_dir == "") {
            printf "./%s\t%s\n", fname, size
        } else {
            printf "%s/%s\t%s\n", current_dir, fname, size
        }
    }
    ' "$tmpfile" | sort

    rm -f "$tmpfile"
}

# ---------------------------------------------------------------------------
# Verify: compare local vs remote manifests
# ---------------------------------------------------------------------------
verify_phenotypes() {
    local pass="$1"

    local local_manifest remote_manifest
    local_manifest=$(mktemp)
    remote_manifest=$(mktemp)

    build_local_manifest "$LOCAL_BASE" > "$local_manifest"
    build_remote_manifest "$pass" "$SFTP_REMOTE_BASE" > "$remote_manifest"

    if diff -q "$local_manifest" "$remote_manifest" > /dev/null 2>&1; then
        log "VERIFY OK: Phenotypes ($(wc -l < "$local_manifest") files match)"
        rm -f "$local_manifest" "$remote_manifest"
        return 0
    else
        log "VERIFY FAILED: Phenotypes"
        log "  Local files:  $(wc -l < "$local_manifest")"
        log "  Remote files: $(wc -l < "$remote_manifest")"
        diff "$local_manifest" "$remote_manifest" | head -20 || true
        rm -f "$local_manifest" "$remote_manifest"
        return 1
    fi
}

# ===========================================================================
# Main
# ===========================================================================
log "=== NAS Download: Phenotypes ==="
log "Mode: ${MODE}"
log "Local: ${LOCAL_BASE}"
log "Remote: ${SFTP_REMOTE_BASE}"
log ""

PASS=$(get_password)
mkdir -p "${STATE_DIR}"

completed_file="${STATE_DIR}/Phenotypes_dl_completed.txt"
failed_file="${STATE_DIR}/Phenotypes_dl_failed.txt"

if [[ "$MODE" == "download" ]]; then
    # Check if already completed
    if [[ -f "$completed_file" ]] && grep -qxF "Phenotypes" "$completed_file" 2>/dev/null; then
        log "SKIP: Phenotypes already marked as completed"
        log "  (Delete $completed_file to re-download)"
        exit 0
    fi

    log "Downloading Phenotypes directory..."

    # Ensure local directory exists
    mkdir -p "$LOCAL_BASE"

    max_attempts=3
    for attempt in $(seq 1 $max_attempts); do
        log "  Attempt $attempt/$max_attempts"

        batchfile=$(mktemp)
        echo "get -Pr ${SFTP_REMOTE_BASE}/* ${LOCAL_BASE}/" > "$batchfile"

        eval "$(sftp_cmd "$PASS" -b "$batchfile")" 2>&1 || true
        rm -f "$batchfile"

        log "  Download attempt $attempt finished"

        if [[ $attempt -lt $max_attempts ]]; then
            log "  Retrying in 10 seconds..."
            sleep 10
        fi
    done

    # Verify
    if verify_phenotypes "$PASS"; then
        echo "Phenotypes" >> "$completed_file"
        log "Download complete and verified."
    else
        echo "Phenotypes" >> "$failed_file"
        log "Download completed but verification FAILED."
        exit 1
    fi
fi

if [[ "$MODE" == "verify" ]]; then
    if verify_phenotypes "$PASS"; then
        log "Verification passed."
    else
        log "Verification FAILED."
        exit 1
    fi
fi

log "=== Done ==="
