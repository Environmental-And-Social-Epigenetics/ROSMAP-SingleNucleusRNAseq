#!/bin/bash
#SBATCH -p pi_lhtsai,pi_manoli
#SBATCH -n 4
#SBATCH --mem=150G
#SBATCH -t 06:00:00
#SBATCH -o logs/%j_donor_qc.out
#SBATCH -e logs/%j_donor_qc.err

set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  SCRIPT_DIR="${SLURM_SUBMIT_DIR}"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

set +u
activate_env "${NEBULA_ENV}"
set -u

export HDF5_USE_FILE_LOCKING=FALSE

mkdir -p "${SCRIPT_DIR}/figures" "${SCRIPT_DIR}/tables" "${SCRIPT_DIR}/logs"

EXTRA_ARGS=""
if [[ "${SMOKE_FLAG:-}" == "--smoke" ]]; then
  EXTRA_ARGS="${EXTRA_ARGS} --smoke"
fi
if [[ -n "${CELLTYPE:-}" ]]; then
  EXTRA_ARGS="${EXTRA_ARGS} --celltype ${CELLTYPE}"
fi

echo "Input H5AD:   ${TSAI_INTEGRATED}/tsai_annotated.h5ad"
echo "Output dir:   ${SCRIPT_DIR}"
echo "Extra args:   ${EXTRA_ARGS:-(none)}"

python "${SCRIPT_DIR}/donor_qc.py" \
  --output-dir "${SCRIPT_DIR}" \
  ${EXTRA_ARGS}
