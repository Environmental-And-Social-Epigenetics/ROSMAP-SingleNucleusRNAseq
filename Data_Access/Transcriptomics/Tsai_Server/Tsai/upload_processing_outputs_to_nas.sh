#!/bin/bash
# =============================================================================
# upload_processing_outputs_to_nas.sh — Upload Tsai Processing Outputs to NAS
# =============================================================================
# Usage:
#   ./upload_processing_outputs_to_nas.sh upload    Upload all pending folders
#   ./upload_processing_outputs_to_nas.sh verify    Verify uploads (manifest check)
#   ./upload_processing_outputs_to_nas.sh probe     List folders and status; no I/O
#
# Transfers Processing_Outputs (QC-filtered, doublet-removed, and integrated/
# annotated h5ad files) to the Tsai Lab NAS for backup and distribution.
#
# Supports sharded parallelism via SFTP_SHARD_TOTAL / SFTP_SHARD_INDEX env vars.
# =============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COHORT="Tsai"

# ---------------------------------------------------------------------------
# Source configuration
# ---------------------------------------------------------------------------
source "$(cd "${SCRIPT_DIR}/../../../../config" && pwd)/paths.sh"
source "${SCRIPT_DIR}/upload_processing_outputs_to_nas.env"

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

string_hash() {
    echo "$1" | cksum | awk '{print $1}'
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

upload_folder() {
    local pass="$1"
    local local_path="$2"
    local remote_base="$3"
    local folder_name="$4"
    local max_attempts=3

    for attempt in $(seq 1 $max_attempts); do
        log "  Attempt $attempt/$max_attempts: uploading $folder_name"

        local batchfile
        batchfile=$(mktemp)
        local _accum=""
        IFS='/' read -ra _parts <<< "${remote_base}/${folder_name}"
        for _p in "${_parts[@]}"; do
            [[ -z "$_p" ]] && continue
            _accum="${_accum}/${_p}"
            echo "-mkdir ${_accum}" >> "$batchfile"
        done
        echo "put -Pr ${local_path} ${remote_base}/" >> "$batchfile"

        if eval "$(sftp_cmd "$pass" -b "$batchfile")" 2>&1 | tee /dev/stderr | grep -q "Uploading\|100%"; then
            rm -f "$batchfile"
            log "  Upload completed for $folder_name (attempt $attempt)"
            return 0
        fi

        rm -f "$batchfile"
        log "  Upload attempt $attempt finished for $folder_name (will verify)"

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
log "=== NAS Transfer: ${COHORT} Processing Outputs ==="
log "Mode: ${MODE}"
log "Local base: ${LOCAL_BASE}"
log "Remote base: ${SFTP_REMOTE_BASE}"
log "Shard: ${SFTP_SHARD_INDEX} of ${SFTP_SHARD_TOTAL}"
log ""

PASS=$(get_password)
mkdir -p "${STATE_DIR}"

total=0; skipped=0; succeeded=0; failed=0; verified=0

for dtype in ${DATA_TYPES}; do
    local_dtype_dir="${LOCAL_BASE}/${dtype}"
    remote_dtype_dir="${SFTP_REMOTE_BASE}/${dtype}"

    completed_file="${STATE_DIR}/${COHORT}_processing_${dtype}_completed.txt"
    failed_file="${STATE_DIR}/${COHORT}_processing_${dtype}_failed.txt"

    if [[ ! -d "$local_dtype_dir" ]]; then
        log "SKIP: Local directory does not exist: $local_dtype_dir"
        continue
    fi

    log "--- Processing stage: ${dtype} ---"

    # Check if directory has subdirectories (per-sample stages) or just files
    # (integrated stage). Handle both patterns.
    has_subdirs=false
    for sub in "${local_dtype_dir}"/*/; do
        if [[ -d "$sub" ]]; then
            has_subdirs=true
            break
        fi
    done

    if [[ "$has_subdirs" == "true" ]]; then
        # Per-sample directory pattern (01_QC_Filtered, 02_Doublet_Removed)
        for folder_path in "${local_dtype_dir}"/*/; do
            [[ -d "$folder_path" ]] || continue

            folder_name=$(basename "$folder_path")
            total=$((total + 1))

            if [[ "${SFTP_SHARD_TOTAL}" -gt 1 ]]; then
                hash_val=$(string_hash "$folder_name")
                if [[ $(( hash_val % SFTP_SHARD_TOTAL )) -ne ${SFTP_SHARD_INDEX} ]]; then
                    continue
                fi
            fi

            if [[ "$MODE" == "probe" ]]; then
                status="PENDING"
                in_list "$folder_name" "$completed_file" 2>/dev/null && status="COMPLETED"
                in_list "$folder_name" "$failed_file" 2>/dev/null && status="FAILED"
                log "  [$status] ${dtype}/${folder_name}"
                continue
            fi

            if [[ "$MODE" == "upload" ]] && in_list "$folder_name" "$completed_file" 2>/dev/null; then
                skipped=$((skipped + 1))
                log "  SKIP (already completed): ${dtype}/${folder_name}"
                continue
            fi

            if [[ "$MODE" == "upload" ]]; then
                log "UPLOAD: ${dtype}/${folder_name}"
                upload_folder "$PASS" "$folder_path" "$remote_dtype_dir" "$folder_name"
                if verify_folder "$PASS" "$folder_path" "${remote_dtype_dir}/${folder_name}" "$folder_name"; then
                    append_state "$folder_name" "$completed_file"
                    succeeded=$((succeeded + 1))
                else
                    append_state "$folder_name" "$failed_file"
                    failed=$((failed + 1))
                fi
            fi

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
    else
        # Flat file pattern (03_Integrated — contains h5ad + CSVs directly)
        total=$((total + 1))

        if [[ "$MODE" == "probe" ]]; then
            status="PENDING"
            in_list "$dtype" "$completed_file" 2>/dev/null && status="COMPLETED"
            in_list "$dtype" "$failed_file" 2>/dev/null && status="FAILED"
            log "  [$status] ${dtype}/"
            continue
        fi

        if [[ "$MODE" == "upload" ]] && in_list "$dtype" "$completed_file" 2>/dev/null; then
            skipped=$((skipped + 1))
            log "  SKIP (already completed): ${dtype}/"
            continue
        fi

        if [[ "$MODE" == "upload" ]]; then
            log "UPLOAD: ${dtype}/ (entire directory)"
            upload_folder "$PASS" "$local_dtype_dir" "${SFTP_REMOTE_BASE}" "$dtype"
            if verify_folder "$PASS" "$local_dtype_dir" "${SFTP_REMOTE_BASE}/${dtype}" "$dtype"; then
                append_state "$dtype" "$completed_file"
                succeeded=$((succeeded + 1))
            else
                append_state "$dtype" "$failed_file"
                failed=$((failed + 1))
            fi
        fi

        if [[ "$MODE" == "verify" ]]; then
            log "VERIFY: ${dtype}/"
            if verify_folder "$PASS" "$local_dtype_dir" "${SFTP_REMOTE_BASE}/${dtype}" "$dtype"; then
                append_state "$dtype" "$completed_file"
                verified=$((verified + 1))
            else
                append_state "$dtype" "$failed_file"
                failed=$((failed + 1))
            fi
        fi
    fi
done

log ""
log "=== Summary ==="
log "Total items scanned: ${total}"
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
