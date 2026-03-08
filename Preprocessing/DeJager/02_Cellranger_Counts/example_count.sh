#!/bin/bash
#
# Example Cell Ranger count command for a single DeJager library.
# Paths are resolved from config/paths.sh environment variables.
#
#SBATCH -t 47:00:00
#SBATCH -n 32
#SBATCH --mem=128G
#SBATCH --mail-type=FAIL

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

export PATH="${CELLRANGER_PATH}:${PATH}"

LIBRARY_ID="190403-B4-A"

cellranger count \
    --create-bam true \
    --include-introns true \
    --nosecondary \
    --r1-length 26 \
    --id "${LIBRARY_ID}" \
    --transcriptome="${CELLRANGER_REF}" \
    --sample "${LIBRARY_ID}_Broad" \
    --fastqs "${DEJAGER_FASTQS}/${LIBRARY_ID}" \
    --output-dir="${DEJAGER_COUNTS}/${LIBRARY_ID}"
