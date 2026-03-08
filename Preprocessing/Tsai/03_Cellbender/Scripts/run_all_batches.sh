#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../Config/cellbender_config.sh"

ensure_pipeline_dirs

shopt -s nullglob
batch_dirs=("${BATCH_SCRIPTS_DIR}"/batch_*)
shopt -u nullglob

if [[ ${#batch_dirs[@]} -eq 0 ]]; then
    echo "ERROR: No batch directories found in ${BATCH_SCRIPTS_DIR}"
    exit 1
fi

for batch_dir in "${batch_dirs[@]}"; do
    batch_name="$(basename "${batch_dir}")"
    batch_num="${batch_name#batch_}"
    "${SCRIPT_DIR}/run_batch.sh" "${batch_num}"
done
