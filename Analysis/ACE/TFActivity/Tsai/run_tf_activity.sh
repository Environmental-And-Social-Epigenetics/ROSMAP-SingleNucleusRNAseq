#!/bin/bash
#SBATCH -p ou_bcs_high,pi_lhtsai,pi_manoli,ou_bcs_low,mit_normal
#SBATCH -n 4
#SBATCH --mem=256G
#SBATCH -t 8:00:00
#SBATCH -o %j.out
#SBATCH -e %j.err

set -euo pipefail

# Under SLURM, BASH_SOURCE points to a temp copy in /var/spool/slurmd/;
# use SLURM_SUBMIT_DIR (the directory where sbatch was run) instead.
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

# -----------------------------------------------------------------------
# Validate arguments
# -----------------------------------------------------------------------

INTEGRATION="${1:?ERROR: integration argument required (e.g. derived_batch)}"
PHENOTYPE="${2:?ERROR: phenotype argument required (tot_adverse_exp, early_hh_ses, or ace_aggregate)}"

# Normalize integration label
case "${INTEGRATION}" in
  batch|derived_batch) INTEGRATION="derived_batch" ;;
  projid) INTEGRATION="projid" ;;
  *)
    echo "ERROR: integration must be one of: derived_batch, projid" >&2
    exit 1
    ;;
esac

# -----------------------------------------------------------------------
# Environment setup
# -----------------------------------------------------------------------

set +u
activate_env "${DECOUPLER_ENV}"
set -u

export HDF5_USE_FILE_LOCKING=FALSE

# -----------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------

INPUT_DIR="${ACE_OUTPUT_ROOT}/DEG/Tsai/celltype_splits_${INTEGRATION}"
OUTPUT_ROOT="${ACE_OUTPUT_ROOT}/TFActivity/Tsai"
RESULTS_DIR="${OUTPUT_ROOT}/results_${INTEGRATION}/${PHENOTYPE}"
FIGURES_DIR="${OUTPUT_ROOT}/figures_${INTEGRATION}/${PHENOTYPE}"
mkdir -p "${RESULTS_DIR}" "${FIGURES_DIR}"

echo "============================================="
echo "ACE TF Activity Analysis"
echo "============================================="
echo "  Integration: ${INTEGRATION}"
echo "  Phenotype:   ${PHENOTYPE}"
echo "  Input dir:   ${INPUT_DIR}"
echo "  Results dir: ${RESULTS_DIR}"
echo "  Figures dir: ${FIGURES_DIR}"
echo "  Pheno CSV:   ${ACE_SCORES_CSV}"
echo "  Conda env:   ${DECOUPLER_ENV}"
echo "============================================="

# -----------------------------------------------------------------------
# Run analysis
# -----------------------------------------------------------------------

"${DECOUPLER_ENV}/bin/python" "${SCRIPT_DIR}/tf_activity_analysis.py" \
  --integration "${INTEGRATION}" \
  --phenotype "${PHENOTYPE}" \
  --input-dir "${INPUT_DIR}" \
  --pheno-csv "${ACE_SCORES_CSV}" \
  --output-dir "${RESULTS_DIR}" \
  --run-progeny

# -----------------------------------------------------------------------
# Run visualization (only if tf_summary.csv was produced)
# -----------------------------------------------------------------------

if [[ -f "${RESULTS_DIR}/tf_summary.csv" ]]; then
  echo ""
  echo "Running visualization..."
  "${DECOUPLER_ENV}/bin/python" "${SCRIPT_DIR}/tf_activity_visualize.py" \
    --results-dir "${RESULTS_DIR}" \
    --output-dir "${FIGURES_DIR}" \
    --phenotype "${PHENOTYPE}" \
    --sex Both
  echo "Visualization complete."
else
  echo "WARNING: tf_summary.csv not found; skipping visualization."
fi

echo ""
echo "Done."
