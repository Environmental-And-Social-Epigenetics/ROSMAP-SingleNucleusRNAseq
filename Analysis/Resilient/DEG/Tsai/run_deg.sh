#!/bin/bash
#SBATCH -p pi_lhtsai,pi_manoli
#SBATCH -n 4
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

normalize_integration() {
  case "$1" in
    batch|derived_batch) echo "derived_batch" ;;
    projid) echo "projid" ;;
    *)
      echo "ERROR: integration argument must be one of: derived_batch, projid" >&2
      exit 1
      ;;
  esac
}

activate_env "${NEBULA_ENV}"

export HDF5_USE_FILE_LOCKING=FALSE

INTEGRATION_RAW="${1:?ERROR: integration argument required (derived_batch or projid)}"
PHENOTYPE="${2:?ERROR: phenotype argument required (resilience_group)}"
CELLTYPE="${3:-}"
INTEGRATION="$(normalize_integration "${INTEGRATION_RAW}")"

OUTPUT_ROOT="${RESILIENT_OUTPUT_ROOT}/DEG/Tsai"
INPUT_DIR="${OUTPUT_ROOT}/celltype_splits_${INTEGRATION}"
RESULTS_DIR="${OUTPUT_ROOT}/results_${INTEGRATION}/${PHENOTYPE}"
mkdir -p "${RESULTS_DIR}"

echo "Integration: ${INTEGRATION}"
echo "Phenotype: ${PHENOTYPE}"
echo "Input dir: ${INPUT_DIR}"
echo "Output dir: ${RESULTS_DIR}"

CMD=(
  Rscript "${SCRIPT_DIR}/resilDegT.Rscript"
  --integration "${INTEGRATION}"
  --phenotype "${PHENOTYPE}"
  --input-dir "${INPUT_DIR}"
  --output-dir "${RESULTS_DIR}"
  --pheno-csv "${RESILIENT_PHENOTYPE_CSV}"
)

if [[ -n "${CELLTYPE}" ]]; then
  CMD+=(--celltype "${CELLTYPE}")
fi

"${CMD[@]}"
