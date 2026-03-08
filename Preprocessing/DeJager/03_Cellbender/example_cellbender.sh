#!/bin/bash
#
# Example CellBender command for a single DeJager library.
# Paths are resolved from config/paths.sh environment variables.
#
#SBATCH -n 32
#SBATCH -t 47:00:00
#SBATCH --gres=gpu:a100:1
#SBATCH --mem=128G
#SBATCH --mail-type=FAIL

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

source "${CONDA_INIT_SCRIPT}"
conda activate "${CELLBENDER_ENV}"

LIBRARY_ID="190403-B4-A"
OUTPUT_DIR="${DEJAGER_PREPROCESSED}/${LIBRARY_ID}"
mkdir -p "${OUTPUT_DIR}"

cellbender remove-background \
    --cuda \
    --input "${DEJAGER_COUNTS}/${LIBRARY_ID}/outs/raw_feature_bc_matrix.h5" \
    --fpr 0 \
    --output "${OUTPUT_DIR}/processed_feature_bc_matrix.h5"
