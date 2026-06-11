#!/bin/bash
#SBATCH -p pi_lhtsai,pi_manoli
#SBATCH -n 4
#SBATCH --mem=200G
#SBATCH -t 24:00:00
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
  # 'primary'/'canonical' resolve to the dataset's declared primary variant
  # (config/variants.yaml). 'batch' is a legacy alias for derived_batch. Any
  # other value is passed through and validated by rosmap_tx.config --resolve,
  # so adding a variant to variants.yaml makes it usable here automatically.
  case "$1" in
    batch) echo "derived_batch" ;;
    primary|canonical|main)
      PYTHONPATH="${REPO_ROOT}/src" python -m rosmap_tx.config --primary tsai ;;
    *)
      PYTHONPATH="${REPO_ROOT}/src" python -m rosmap_tx.config --resolve tsai "$1" 2>/dev/null \
        || { echo "ERROR: unknown integration/variant '$1' (see config/variants.yaml)" >&2; exit 1; }
      ;;
  esac
}

set +u
activate_env "${NEBULA_ENV}"
set -u

export HDF5_USE_FILE_LOCKING=FALSE

INTEGRATION_RAW="${1:?ERROR: integration argument required (e.g. primary, derived_batch, projid)}"
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
