#!/bin/bash
#SBATCH -p pi_lhtsai,pi_manoli
#SBATCH -n 4
#SBATCH --mem=100G
#SBATCH -t 04:00:00
#SBATCH -o logs/%j_marker_qc.out
#SBATCH -e logs/%j_marker_qc.err

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

INPUT_H5AD="${INPUT_H5AD:-${TSAI_INTEGRATED}/tsai_annotated.h5ad}"
PER_CT_CAP="${PER_CT_CAP:-20000}"
THRESHOLD="${THRESHOLD:-0.30}"
SEED="${SEED:-0}"

echo "Input H5AD:       ${INPUT_H5AD}"
echo "Output dir:       ${SCRIPT_DIR}"
echo "Per-celltype cap: ${PER_CT_CAP}"
echo "Threshold:        ${THRESHOLD}"
echo "Seed:             ${SEED}"

python "${SCRIPT_DIR}/make_marker_dotplots.py" \
  --input-h5ad "${INPUT_H5AD}" \
  --output-dir "${SCRIPT_DIR}" \
  --per-celltype-cap "${PER_CT_CAP}" \
  --threshold "${THRESHOLD}" \
  --seed "${SEED}"
