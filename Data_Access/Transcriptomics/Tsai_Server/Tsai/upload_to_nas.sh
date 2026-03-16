#!/bin/bash
# =============================================================================
# upload_to_nas.sh — Upload Tsai Transcriptomics data to Tsai Lab NAS
# =============================================================================
# Usage:
#   ./upload_to_nas.sh upload    Upload all pending folders to NAS
#   ./upload_to_nas.sh verify    Verify already-uploaded folders (manifest check)
#   ./upload_to_nas.sh probe     List folders and their transfer status; no I/O
#
# Supports sharded parallelism via SFTP_SHARD_TOTAL / SFTP_SHARD_INDEX env vars
# (see slurm_upload.sh for the SLURM array wrapper).
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COHORT="Tsai"

# ---------------------------------------------------------------------------
# Source configuration
# ---------------------------------------------------------------------------
source "$(cd "${SCRIPT_DIR}/../../../../config" && pwd)/paths.sh"
source "${SCRIPT_DIR}/upload_to_nas.env"

# ---------------------------------------------------------------------------
# Resolve mode
# ---------------------------------------------------------------------------
MODE="${1:-upload}"
if [[ "$MODE" != "upload" && "$MODE" != "verify" && "$MODE" != "probe" ]]; then
    echo "ERROR: Unknown mode '$MODE'. Use: upload | verify | probe"
    exit 1
fi

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

# Extract password from SMB-format auth file
get_password() {
    local pass
    pass=$(grep -i '^password' "$SFTP_AUTHFILE" | sed 's/^[^=]*= *//')
    if [[ -z "$pass" ]]; then
        log "ERROR: Could not extract password from $SFTP_AUTHFILE"
        exit 1
    fi
    echo "$pass"
}

# Build the base sftp command (without batch flags)
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

# Deterministic hash of a string -> integer (portable)
string_hash() {
    echo "$1" | cksum | awk '{print $1}'
}

# Check if a folder name is in a file (one name per line)
in_list() {
    local name="$1" file="$2"
    [[ -f "$file" ]] && grep -qxF "$name" "$file"
}

# Append a folder name to a state file (idempotent)
append_state() {
    local name="$1" file="$2"
    mkdir -p "$(dirname "$file")"
    in_list "$name" "$file" 2>/dev/null || echo "$name" >> "$file"
}

# ---------------------------------------------------------------------------
# Build local manifest: relative_path<TAB>size_bytes
# ---------------------------------------------------------------------------
build_local_manifest() {
    local dir="$1"
    # Using find + stat; output sorted for stable diffing
    (cd "$dir" && find . -type f -exec stat --format='%n\t%s' {} + | sort)
}

# ---------------------------------------------------------------------------
# Build remote manifest from sftp "ls -l" output.
# The NAS reports WIN\Domain Users (group name with a space), so we cannot
# rely on a fixed column index for the size field.  Instead we find the
# column dynamically: size is the first purely-numeric field after the
# permission string (column 1).
# ---------------------------------------------------------------------------
build_remote_manifest() {
    local pass="$1"
    local remote_dir="$2"
    local tmpfile
    tmpfile=$(mktemp)

    # Batch: recursive ls
    local batchfile
    batchfile=$(mktemp)
    echo "ls -lR ${remote_dir}" > "$batchfile"

    eval "$(sftp_cmd "$pass" -b "$batchfile")" > "$tmpfile" 2>/dev/null || true
    rm -f "$batchfile"

    # Parse: extract relative paths and sizes
    awk -v base="$remote_dir" '
    /^sftp>/ { next }
    /^Listing/ { next }
    /^$/ { next }
    # Directory header lines from ls -lR
    /^\// {
        sub(/:$/, "")
        current_dir = $0
        # Make relative to base
        sub("^" base "/?", "./", current_dir)
        next
    }
    # Skip directory entries (first char is d) and total lines
    /^d/ { next }
    /^total/ { next }
    /^-/ {
        # Find size: first purely numeric field after permissions (field 1)
        size = ""
        for (i = 2; i <= NF; i++) {
            if ($i ~ /^[0-9]+$/) {
                size = $i
                break
            }
        }
        if (size == "") next
        # Filename is the last field
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
# Verify a single folder: compare local vs remote manifests
# Returns 0 on match, 1 on mismatch
# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
# Upload a single folder with retry
# ---------------------------------------------------------------------------
upload_folder() {
    local pass="$1"
    local local_path="$2"     # full path to local folder
    local remote_base="$3"    # remote directory to upload INTO
    local folder_name="$4"
    local max_attempts=3

    for attempt in $(seq 1 $max_attempts); do
        log "  Attempt $attempt/$max_attempts: uploading $folder_name"

        # Build batch file
        local batchfile
        batchfile=$(mktemp)
        # Ensure full remote directory tree exists (walk each path component)
        local _accum=""
        IFS='/' read -ra _parts <<< "${remote_base}/${folder_name}"
        for _p in "${_parts[@]}"; do
            [[ -z "$_p" ]] && continue
            _accum="${_accum}/${_p}"
            echo "-mkdir ${_accum}" >> "$batchfile"
        done
        echo "put -Pr ${local_path} ${remote_base}/" >> "$batchfile"

        local rc=0
        eval "$(sftp_cmd "$pass" -b "$batchfile")" 2>&1 || rc=$?
        rm -f "$batchfile"

        if [[ $rc -eq 0 ]]; then
            log "  Upload completed for $folder_name (attempt $attempt)"
            return 0
        fi

        log "  Upload attempt $attempt failed for $folder_name (exit code $rc)"

        if [[ $attempt -lt $max_attempts ]]; then
            log "  Retrying in 10 seconds..."
            sleep 10
        fi
    done

    log "  All $max_attempts attempts failed for $folder_name"
    return 1
}

# ===========================================================================
# Main
# ===========================================================================
log "=== NAS Transfer: ${COHORT} Transcriptomics ==="
log "Mode: ${MODE}"
log "Local base: ${LOCAL_BASE}"
log "Remote base: ${SFTP_REMOTE_BASE}"
log "Shard: ${SFTP_SHARD_INDEX} of ${SFTP_SHARD_TOTAL}"
log ""

# Get password once
PASS=$(get_password)

# Ensure state directory exists
mkdir -p "${STATE_DIR}"

# Counters
total=0; skipped=0; succeeded=0; failed=0; verified=0

for dtype in ${DATA_TYPES}; do
    local_dtype_dir="${LOCAL_BASE}/${dtype}"
    remote_dtype_dir="${SFTP_REMOTE_BASE}/${dtype}"

    completed_file="${STATE_DIR}/${COHORT}_${dtype}_completed.txt"
    failed_file="${STATE_DIR}/${COHORT}_${dtype}_failed.txt"

    if [[ ! -d "$local_dtype_dir" ]]; then
        log "SKIP: Local directory does not exist: $local_dtype_dir"
        continue
    fi

    log "--- Processing data type: ${dtype} ---"

    # Iterate over subdirectories (each is a library/sample folder)
    for folder_path in "${local_dtype_dir}"/*/; do
        # Skip if not a directory
        [[ -d "$folder_path" ]] || continue

        folder_name=$(basename "$folder_path")
        total=$((total + 1))

        # Sharding: skip folders not assigned to this worker
        if [[ "${SFTP_SHARD_TOTAL}" -gt 1 ]]; then
            hash_val=$(string_hash "$folder_name")
            if [[ $(( hash_val % SFTP_SHARD_TOTAL )) -ne ${SFTP_SHARD_INDEX} ]]; then
                continue
            fi
        fi

        # --- PROBE mode ---
        if [[ "$MODE" == "probe" ]]; then
            status="PENDING"
            in_list "$folder_name" "$completed_file" 2>/dev/null && status="COMPLETED"
            in_list "$folder_name" "$failed_file" 2>/dev/null && status="FAILED"
            log "  [$status] ${dtype}/${folder_name}"
            continue
        fi

        # Skip already-completed folders (upload mode only)
        if [[ "$MODE" == "upload" ]] && in_list "$folder_name" "$completed_file" 2>/dev/null; then
            skipped=$((skipped + 1))
            log "  SKIP (already completed): ${dtype}/${folder_name}"
            continue
        fi

        # --- UPLOAD mode ---
        if [[ "$MODE" == "upload" ]]; then
            log "UPLOAD: ${dtype}/${folder_name}"
            if upload_folder "$PASS" "$folder_path" "$remote_dtype_dir" "$folder_name"; then
                append_state "$folder_name" "$completed_file"
                succeeded=$((succeeded + 1))
            else
                append_state "$folder_name" "$failed_file"
                failed=$((failed + 1))
                log "  FAILED: ${dtype}/${folder_name}"
            fi
        fi

        # --- VERIFY mode ---
        if [[ "$MODE" == "verify" ]]; then
            log "VERIFY: ${dtype}/${folder_name}"
            if verify_folder "$PASS" "$folder_path" "${remote_dtype_dir}/${folder_name}" "$folder_name"; then
                append_state "$folder_name" "$completed_file"
                verified=$((verified + 1))
            else
                append_state "$folder_name" "$failed_file"
                failed=$((failed + 1))
            fi
        fi
    done
done

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
log ""
log "=== Summary ==="
log "Total folders scanned: ${total}"
if [[ "$MODE" == "upload" ]]; then
    log "Skipped (already done): ${skipped}"
    log "Uploaded + verified:    ${succeeded}"
    log "Failed:                 ${failed}"
elif [[ "$MODE" == "verify" ]]; then
    log "Verified OK:            ${verified}"
    log "Failed verification:    ${failed}"
fi
log "=== Done ==="

exit $( [[ $failed -eq 0 ]] && echo 0 || echo 1 )
