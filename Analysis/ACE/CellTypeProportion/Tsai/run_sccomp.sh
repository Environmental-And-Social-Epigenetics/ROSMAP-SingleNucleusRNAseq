#!/bin/bash
#SBATCH -p pi_lhtsai,pi_manoli
#SBATCH -n 8
#SBATCH -t 08:00:00
#SBATCH --mem=64G
#SBATCH -o %j.out
#SBATCH -e %j.err

set -euo pipefail

# Under SLURM, BASH_SOURCE points to a temp copy in /var/spool/slurmd/;
# use SLURM_SUBMIT_DIR (the directory where sbatch was run) instead.
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  SCRIPT_DIR="${SLURM_SUBMIT_DIR}"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

INTEGRATION="${1:?ERROR: integration argument required (derived_batch or projid)}"
SEX="${2:?ERROR: sex argument required (all, male, or female)}"
RESOLUTION="${3:?ERROR: resolution argument required (fine or broad)}"
shift 3 || true

OUTPUT_ROOT="${ACE_PROP_OUTPUT_ROOT:-${ANALYSIS_OUTPUT_ROOT}/ACE/CellTypeProportion/Tsai}"

set +u
activate_env "${SCCOMP_ENV}"
set -u
export HDF5_USE_FILE_LOCKING=FALSE

"${SCCOMP_ENV}/bin/Rscript" "${SCRIPT_DIR}/sccomp_analysis.R" \
    --integration "${INTEGRATION}" \
    --sex "${SEX}" \
    --resolution "${RESOLUTION}" \
    --output-root "${OUTPUT_ROOT}" \
    "$@"
