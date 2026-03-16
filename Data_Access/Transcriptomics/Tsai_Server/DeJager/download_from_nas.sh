#!/bin/bash
# =============================================================================
# download_from_nas.sh — Download DeJager Transcriptomics data from Tsai Lab NAS
# =============================================================================
# Usage:
#   ./download_from_nas.sh              Download all remote folders not yet local
#   ./download_from_nas.sh verify       Verify downloaded data against NAS
#
# Downloads each folder within each DATA_TYPE from the NAS to the local
# Data/Transcriptomics/DeJager/ tree.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COHORT="DeJager"

# ---------------------------------------------------------------------------
# Source configuration
# ---------------------------------------------------------------------------
source "$(cd "${SCRIPT_DIR}/../../../../config" && pwd)/paths.sh"
source "${SCRIPT_DIR}/upload_to_nas.env"

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

in_list() {
    local name="$1" file="$2"
    [[ -f "$file" ]] && grep -qxF "$name" "$file"
}

append_state() {
    local name="$1" file="$2"
    mkdir -p "$(dirname "$file")"
    in_list "$name" "$file" 2>/dev/null || echo "$name" >> "$file"
}

build_local_manifest() {
    local dir="$1"
    (cd "$dir" && find . -type f -exec stat --format='%n\t%s' {} + | sort)
}

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

list_remote_folders() {
    local pass="$1"
    local remote_dir="$2"
    local tmpfile
    tmpfile=$(mktemp)

    local batchfile
    batchfile=$(mktemp)
    echo "ls -1 ${remote_dir}" > "$batchfile"

    eval "$(sftp_cmd "$pass" -b "$batchfile")" > "$tmpfile" 2>/dev/null || true
    rm -f "$batchfile"

    grep -v '^sftp>' "$tmpfile" | grep -v '^\.' | sed '/^$/d' | while read -r name; do
        echo "$name"
    done

    rm -f "$tmpfile"
}

verify_folder() {
    local pass="$1"
    local local_dir="$2"
    local remote_dir="$3"
    local folder_name="$4"

    local local_manifest remote_manifest
    local_manifest=$(mktemp)
    remote_manifest=$(mktemp)

    build_local_manifest "$local_dir" > "$local_manifest"
    build_remote_manifest "$pass" "$remote_dir" > "$remote_manifest"

    if diff -q "$local_manifest" "$remote_manifest" > /dev/null 2>&1; then
        log "  VERIFY OK: $folder_name ($(wc -l < "$local_manifest") files match)"
        rm -f "$local_manifest" "$remote_manifest"
        return 0
    else
        log "  VERIFY FAILED: $folder_name"
        log "  Local files:  $(wc -l < "$local_manifest")"
        log "  Remote files: $(wc -l < "$remote_manifest")"
        diff "$local_manifest" "$remote_manifest" | head -20 || true
        rm -f "$local_manifest" "$remote_manifest"
        return 1
    fi
}

download_folder() {
    local pass="$1"
    local remote_path="$2"
    local local_base="$3"
    local folder_name="$4"
    local max_attempts=3

    mkdir -p "${local_base}"

    for attempt in $(seq 1 $max_attempts); do
        log "  Attempt $attempt/$max_attempts: downloading $folder_name"

        local batchfile
        batchfile=$(mktemp)
        echo "get -Pr ${remote_path} ${local_base}/" > "$batchfile"

        eval "$(sftp_cmd "$pass" -b "$batchfile")" 2>&1 || true
        rm -f "$batchfile"

        log "  Download attempt $attempt finished for $folder_name"

        if [[ $attempt -lt $max_attempts ]]; then
            log "  Retrying in 10 seconds..."
            sleep 10
        fi
    done

    return 0
}

# ===========================================================================
# Main
# ===========================================================================
log "=== NAS Download: ${COHORT} Transcriptomics ==="
log "Mode: ${MODE}"
log "Local base: ${LOCAL_BASE}"
log "Remote base: ${SFTP_REMOTE_BASE}"
log ""

PASS=$(get_password)
mkdir -p "${STATE_DIR}"

total=0; skipped=0; succeeded=0; failed=0; verified=0

for dtype in ${DATA_TYPES}; do
    local_dtype_dir="${LOCAL_BASE}/${dtype}"
    remote_dtype_dir="${SFTP_REMOTE_BASE}/${dtype}"

    completed_file="${STATE_DIR}/${COHORT}_${dtype}_dl_completed.txt"
    failed_file="${STATE_DIR}/${COHORT}_${dtype}_dl_failed.txt"

    log "--- Processing data type: ${dtype} ---"

    remote_folders=$(list_remote_folders "$PASS" "$remote_dtype_dir")

    if [[ -z "$remote_folders" ]]; then
        log "  No remote folders found in ${remote_dtype_dir}"
        continue
    fi

    while IFS= read -r folder_name; do
        [[ -z "$folder_name" ]] && continue
        total=$((total + 1))

        # --- DOWNLOAD mode ---
        if [[ "$MODE" == "download" ]]; then
            if in_list "$folder_name" "$completed_file" 2>/dev/null; then
                skipped=$((skipped + 1))
                log "  SKIP (already completed): ${dtype}/${folder_name}"
                continue
            fi

            log "DOWNLOAD: ${dtype}/${folder_name}"
            download_folder "$PASS" "${remote_dtype_dir}/${folder_name}" "$local_dtype_dir" "$folder_name"

            if [[ -d "${local_dtype_dir}/${folder_name}" ]] && \
               verify_folder "$PASS" "${local_dtype_dir}/${folder_name}" "${remote_dtype_dir}/${folder_name}" "$folder_name"; then
                append_state "$folder_name" "$completed_file"
                succeeded=$((succeeded + 1))
            else
                append_state "$folder_name" "$failed_file"
                failed=$((failed + 1))
                log "  FAILED verification: ${dtype}/${folder_name}"
            fi
        fi

        # --- VERIFY mode ---
        if [[ "$MODE" == "verify" ]]; then
            log "VERIFY: ${dtype}/${folder_name}"
            if [[ -d "${local_dtype_dir}/${folder_name}" ]] && \
               verify_folder "$PASS" "${local_dtype_dir}/${folder_name}" "${remote_dtype_dir}/${folder_name}" "$folder_name"; then
                append_state "$folder_name" "$completed_file"
                verified=$((verified + 1))
            else
                append_state "$folder_name" "$failed_file"
                failed=$((failed + 1))
            fi
        fi
    done <<< "$remote_folders"
done

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
log ""
log "=== Summary ==="
log "Total remote folders found: ${total}"
if [[ "$MODE" == "download" ]]; then
    log "Skipped (already done): ${skipped}"
    log "Downloaded + verified:  ${succeeded}"
    log "Failed:                 ${failed}"
elif [[ "$MODE" == "verify" ]]; then
    log "Verified OK:            ${verified}"
    log "Failed verification:    ${failed}"
fi
log "=== Done ==="

exit $( [[ $failed -eq 0 ]] && echo 0 || echo 1 )
