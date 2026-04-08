#!/bin/bash
#SBATCH -p pi_lhtsai,pi_manoli
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH --mem=16G
#SBATCH -o %j.out
#SBATCH -e %j.err

set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  SCRIPT_DIR="${SLURM_SUBMIT_DIR}"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

INTEGRATION="${1:?ERROR: integration argument required (derived_batch or projid)}"
shift || true

OUTPUT_ROOT="${RESILIENT_OUTPUT_ROOT}/CellTypeProportion/Tsai"

set +u
activate_env "${NEBULA_ENV}"
set -u
export HDF5_USE_FILE_LOCKING=FALSE

"${NEBULA_ENV}/bin/python" "${SCRIPT_DIR}/prep_counts.py" \
    --integration "${INTEGRATION}" \
    --output-root "${OUTPUT_ROOT}" \
    "$@"
