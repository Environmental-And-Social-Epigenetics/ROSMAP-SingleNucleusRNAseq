#!/bin/bash
#SBATCH -p pi_lhtsai,pi_manoli
#SBATCH -n 4
#SBATCH -o %j.out
#SBATCH -e %j.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

activate_env "${NEBULA_ENV}"

export HDF5_USE_FILE_LOCKING=FALSE

INTEGRATION="${1:?ERROR: integration argument required}"
PHENOTYPE="${2:?ERROR: phenotype argument required}"
CELLTYPE="${3:-}"

OUTPUT_ROOT="${ANALYSIS_OUTPUT_ROOT}/ACE/DEG/DeJager"
INPUT_DIR="${OUTPUT_ROOT}/celltype_splits_${INTEGRATION}"
RESULTS_DIR="${OUTPUT_ROOT}/results_${INTEGRATION}/${PHENOTYPE}"
mkdir -p "${RESULTS_DIR}"

CMD=(
  Rscript "${SCRIPT_DIR}/aceDegDJ.Rscript"
  --integration "${INTEGRATION}"
  --phenotype "${PHENOTYPE}"
  --input-dir "${INPUT_DIR}"
  --output-dir "${RESULTS_DIR}"
  --pheno-csv "${ACE_SCORES_CSV}"
)

if [[ -n "${CELLTYPE}" ]]; then
  CMD+=(--celltype "${CELLTYPE}")
fi

"${CMD[@]}"
