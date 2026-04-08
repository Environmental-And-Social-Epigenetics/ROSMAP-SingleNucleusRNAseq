#!/bin/bash
#SBATCH -n 45
#SBATCH -t 5:00:00
#SBATCH --mem=100G
#SBATCH -o gsea_%j.out
#SBATCH -e gsea_%j.err
#SBATCH -p pi_lhtsai,pi_manoli

set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  SCRIPT_DIR="${SLURM_SUBMIT_DIR}"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

activate_env "${GSEA_ANALYSIS_ENV}"

OUTPUT_DIR="${ACE_OUTPUT_ROOT}/GSEA/Tsai"
mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}"

Rscript "${SCRIPT_DIR}/gseaResults.Rscript"
