#!/bin/bash
#SBATCH -p pi_lhtsai,pi_manoli
#SBATCH -n 4
#SBATCH -t 01:00:00
#SBATCH --mem=16G
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
shift || true

OUTPUT_ROOT="${ACE_PROP_OUTPUT_ROOT:-${ANALYSIS_OUTPUT_ROOT}/ACE/CellTypeProportion/Tsai}"

set +u
activate_env "${SCCOMP_ENV}"
set -u

"${SCCOMP_ENV}/bin/Rscript" "${SCRIPT_DIR}/sccomp_visualize.R" \
    --integration "${INTEGRATION}" \
    --output-root "${OUTPUT_ROOT}" \
    "$@"
