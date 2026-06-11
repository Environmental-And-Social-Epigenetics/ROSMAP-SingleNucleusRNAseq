#!/bin/bash
#SBATCH -p pi_lhtsai,pi_manoli,ou_bcs_low,mit_preemptable
#SBATCH -n 20
#SBATCH --mem=200G
#SBATCH -t 24:00:00
#SBATCH -o %j.out
#SBATCH -e %j.err

set -euo pipefail

# ---------------------------------------------------------------------------
# Resolve script directory (works under SLURM and direct invocation)
# ---------------------------------------------------------------------------
# Prefer LAUNCHER_SCRIPT_DIR (exported from launcher) over SLURM_SUBMIT_DIR
if [[ -n "${LAUNCHER_SCRIPT_DIR:-}" ]]; then
  SCRIPT_DIR="${LAUNCHER_SCRIPT_DIR}"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  SCRIPT_DIR="${SLURM_SUBMIT_DIR}"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

# ---------------------------------------------------------------------------
# Arguments (passed by aceScenicT.sh / aceScenicDJ.sh or manually)
# ---------------------------------------------------------------------------
INTEGRATION="${1:?ERROR: integration argument required (e.g. derived_batch)}"
PHENOTYPE="${2:?ERROR: phenotype argument required (e.g. tot_adverse_exp)}"
CELL_TYPE="${3:?ERROR: cell_type argument required (e.g. Mic)}"
SEX="${4:?ERROR: sex argument required (Male or Female)}"
# Cohort selects which cohort-specific DEG inputs / SCENIC output dir to use.
# It does NOT change the analysis method (micropool pool_size, GRNBoost2,
# cisTarget, AUCell, OLS are identical). Default keeps existing Tsai behaviour.
COHORT="${5:-tsai}"

# Validate sex
if [[ "${SEX}" != "Male" && "${SEX}" != "Female" ]]; then
  echo "ERROR: sex must be Male or Female, got: ${SEX}" >&2
  exit 1
fi

# Validate cohort
case "${COHORT}" in
  tsai)    COHORT_DIR="Tsai" ;;
  dejager) COHORT_DIR="DeJager" ;;
  *)
    echo "ERROR: cohort must be tsai or dejager, got: ${COHORT}" >&2
    exit 1
    ;;
esac

# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------
set +u
activate_env "${SCENIC_ANALYSIS_ENV}"
set -u

export HDF5_USE_FILE_LOCKING=FALSE

# ---------------------------------------------------------------------------
# Paths (cohort-specific DEG inputs and SCENIC outputs; shared method)
# ---------------------------------------------------------------------------
INPUT_DIR="${ACE_OUTPUT_ROOT}/DEG/${COHORT_DIR}/celltype_splits_${INTEGRATION}"
INPUT_H5AD="${INPUT_DIR}/${CELL_TYPE}.h5ad"
OUTPUT_DIR="${ACE_OUTPUT_ROOT}/SCENIC/${COHORT_DIR}/results_${INTEGRATION}/${PHENOTYPE}/${SEX}_${CELL_TYPE}"
mkdir -p "${OUTPUT_DIR}"

# ---------------------------------------------------------------------------
# Validate inputs
# ---------------------------------------------------------------------------
if [[ ! -f "${INPUT_H5AD}" ]]; then
  echo "ERROR: Input h5ad not found: ${INPUT_H5AD}" >&2
  exit 1
fi

if [[ ! -f "${ACE_SCORES_CSV}" ]]; then
  echo "ERROR: ACE scores CSV not found: ${ACE_SCORES_CSV}" >&2
  exit 1
fi

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
echo "=== ACE SCENIC: ${CELL_TYPE} / ${SEX} (cohort=${COHORT}) ==="
echo "Integration: ${INTEGRATION}"
echo "Phenotype:   ${PHENOTYPE}"
echo "Input h5ad:  ${INPUT_H5AD}"
echo "Output dir:  ${OUTPUT_DIR}"
echo "Ranking dir: ${SCENIC_RANKING_DIR}"
echo "TF list:     ${SCENIC_TF_LIST}"
echo "Pool size:   ${POOL_SIZE:-50}"
echo ""

"${SCENIC_ANALYSIS_ENV}/bin/python" "${SCRIPT_DIR}/scenic_analysis.py" \
  --cohort "${COHORT}" \
  --cell-type "${CELL_TYPE}" \
  --sex "${SEX}" \
  --phenotype "${PHENOTYPE}" \
  --integration "${INTEGRATION}" \
  --input-h5ad "${INPUT_H5AD}" \
  --pheno-csv "${ACE_SCORES_CSV}" \
  --output-dir "${OUTPUT_DIR}" \
  --ranking-dir "${SCENIC_RANKING_DIR}" \
  --tf-list "${SCENIC_TF_LIST}" \
  --pool-size "${POOL_SIZE:-50}" \
  --num-workers 4

echo "=== Done: ${CELL_TYPE} / ${SEX} (cohort=${COHORT}) ==="
