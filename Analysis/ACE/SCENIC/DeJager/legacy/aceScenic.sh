#!/bin/bash
#SBATCH -n 40
#SBATCH -t 48:00:00
#SBATCH --mem=600G
#SBATCH -o aceScenic_%j.out
#SBATCH -e aceScenic_%j.err
#SBATCH -p pi_lhtsai,pi_manoli

set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  SCRIPT_DIR="${SLURM_SUBMIT_DIR}"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

activate_env "${SCENIC_ANALYSIS_ENV}"

export HDF5_USE_FILE_LOCKING=FALSE

mkdir -p "${ACE_OUTPUT_ROOT}/SCENIC/DeJager"
cd "${ACE_OUTPUT_ROOT}/SCENIC/DeJager"

python3 "${SCRIPT_DIR}/aceScenic.py"
