#!/bin/bash
#SBATCH -p pi_lhtsai,pi_manoli
#SBATCH -n 4
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
PHENOTYPE="${2:?ERROR: phenotype argument required (tot_adverse_exp, early_hh_ses, or ace_aggregate)}"
CELLTYPE="${3:-}"
INTEGRATION="$(normalize_integration "${INTEGRATION_RAW}")"

OUTPUT_ROOT="${ANALYSIS_OUTPUT_ROOT}/ACE/DEG/Tsai"
INPUT_DIR="${OUTPUT_ROOT}/celltype_splits_${INTEGRATION}"
RESULTS_DIR="${OUTPUT_ROOT}/results_${INTEGRATION}/${PHENOTYPE}"
mkdir -p "${RESULTS_DIR}"

echo "Integration: ${INTEGRATION}"
echo "Phenotype: ${PHENOTYPE}"
echo "Input dir: ${INPUT_DIR}"
echo "Output dir: ${RESULTS_DIR}"

CMD=(
  Rscript "${SCRIPT_DIR}/aceDegT.Rscript"
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
