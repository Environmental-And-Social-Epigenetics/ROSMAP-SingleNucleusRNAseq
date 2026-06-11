#!/bin/bash
# DEPRECATED: superseded by the modular scenic_analysis.py pipeline; retained
# for reference only.
#
# The supported entry point is Analysis/ACE/SCENIC/DeJager/aceScenicDJ.sh ->
# run_scenic.sh, which delegates to the shared cohort-aware
# Analysis/ACE/SCENIC/Tsai/scenic_analysis.py (--cohort dejager, single
# pool_size). Both this wrapper and aceScenic.py were relocated into legacy/,
# so the REPO_ROOT/aceScenic.py paths below are no longer correct relative to
# this directory. Do not run this script.
#SBATCH -n 40
#SBATCH -t 48:00:00
#SBATCH --mem=600G
#SBATCH -o aceScenic_%j.out
#SBATCH -e aceScenic_%j.err
#SBATCH -p pi_lhtsai,pi_manoli

echo "DEPRECATED: use Analysis/ACE/SCENIC/DeJager/aceScenicDJ.sh (modular pipeline)." >&2
exit 1

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
