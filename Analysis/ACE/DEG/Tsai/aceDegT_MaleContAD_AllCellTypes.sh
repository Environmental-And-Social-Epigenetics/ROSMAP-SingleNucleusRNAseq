#!/bin/bash
#SBATCH -p pi_lhtsai,pi_manoli
#SBATCH -n 4
#SBATCH --mem=200G
#SBATCH -t 12:00:00
#SBATCH -o %j_MaleContAD_AllCellTypes.out
#SBATCH -e %j_MaleContAD_AllCellTypes.err

set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  SCRIPT_DIR="${SLURM_SUBMIT_DIR}"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

set +u
activate_env "${NEBULA_ENV}"
set -u

export HDF5_USE_FILE_LOCKING=FALSE

INTEGRATION="derived_batch"
PHENOTYPE="tot_adverse_exp"
OUTPUT_ROOT="${ANALYSIS_OUTPUT_ROOT}/ACE/DEG/Tsai"
INPUT_DIR="${OUTPUT_ROOT}/celltype_splits_${INTEGRATION}"
RESULTS_DIR="${OUTPUT_ROOT}/results_${INTEGRATION}_MaleContAD_AllCellTypes/${PHENOTYPE}"
mkdir -p "${RESULTS_DIR}"

echo "Integration:   ${INTEGRATION}"
echo "Phenotype:     ${PHENOTYPE}"
echo "Cohort:        Males only, continuous AD adjustment"
echo "Design:        ~ age_death + pmi + amylsqrt + tangsqrt + tot_adverse_exp"
echo "Cell types:    ALL (broad_Exc, broad_Inh, + 17 subtypes)"
echo "Input dir:     ${INPUT_DIR}"
echo "Output dir:    ${RESULTS_DIR}"
echo "Smoke flag:    ${SMOKE_FLAG:-(none)}"
echo "Celltype arg:  ${CELLTYPE:-(none, all)}"

CMD=(
  Rscript "${SCRIPT_DIR}/aceDegT_MaleContAD_AllCellTypes.Rscript"
  --integration "${INTEGRATION}"
  --phenotype "${PHENOTYPE}"
  --input-dir "${INPUT_DIR}"
  --output-dir "${RESULTS_DIR}"
  --pheno-csv "${ACE_SCORES_CSV}"
)
if [[ -n "${CELLTYPE:-}" ]]; then
  CMD+=(--celltype "${CELLTYPE}")
fi
if [[ "${SMOKE_FLAG:-}" == "--smoke" ]]; then
  CMD+=(--smoke)
fi

"${CMD[@]}"
