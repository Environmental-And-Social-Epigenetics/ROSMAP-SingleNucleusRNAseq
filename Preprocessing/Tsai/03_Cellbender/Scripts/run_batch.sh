#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../Config/cellbender_config.sh"

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <batch_num>"
    exit 1
fi

ensure_pipeline_dirs

batch_num="$1"
batch_dir="${BATCH_SCRIPTS_DIR}/batch_${batch_num}/cellbender"

if [[ ! -d "${batch_dir}" ]]; then
    echo "ERROR: Batch directory not found: ${batch_dir}"
    exit 1
fi

shopt -s nullglob
scripts=("${batch_dir}"/*.sh)
shopt -u nullglob

if [[ ${#scripts[@]} -eq 0 ]]; then
    echo "ERROR: No CellBender scripts found in ${batch_dir}"
    exit 1
fi

echo "Submitting ${#scripts[@]} CellBender jobs for batch ${batch_num}..."
for script in "${scripts[@]}"; do
    sbatch "${script}"
done
