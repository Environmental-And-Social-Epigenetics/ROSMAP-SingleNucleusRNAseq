#!/bin/bash
# TF / pathway activity for the non-ANCOVA male ACE AD-model arms.
# One SLURM job per arm; each infers TF (DoRothEA/CollecTRI) + PROGENy activity on
# the signal-bearing cell types (males) and runs the per-arm OLS (covariate set
# matches the arm). Output:
#   ${ACE_OUTPUT_ROOT}/TFActivity/Tsai/results_derived_batch_<ARM>/tot_adverse_exp/
#
# Activity inference itself is cheap (decoupler MLM on pseudobulk), so each arm
# recomputes activities then applies its own formula -- no shared cache needed.
#
# Usage:
#   bash aceTfActT_male_arms.sh
#   SMOKE_FLAG=--smoke bash aceTfActT_male_arms.sh
#   ARMS_OVERRIDE="MaleContAD" bash aceTfActT_male_arms.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"
source "${REPO_ROOT}/Analysis/ACE/_shared/arm_covariates.sh"

INTEGRATION="derived_batch"
PHENOTYPE="tot_adverse_exp"
INPUT_DIR="${ACE_OUTPUT_ROOT}/DEG/Tsai/celltype_splits_${INTEGRATION}"
OUTPUT_ROOT="${ACE_OUTPUT_ROOT}/TFActivity/Tsai"
LOG_DIR="${OUTPUT_ROOT}/logs"
mkdir -p "${LOG_DIR}"

PART="-p pi_lhtsai,pi_manoli"
SMOKE_FLAG="${SMOKE_FLAG:-}"
CELLTYPES_CSV="$(IFS=,; echo "${SIGNAL_CELLTYPES[*]}")"

read -r -a ARMS <<< "${ARMS_OVERRIDE:-${NON_ANCOVA_ARMS[*]}}"

echo "=== ACE TF activity: male AD-model arms ==="
echo "Arms:       ${ARMS[*]}"
echo "Cell types: ${CELLTYPES_CSV}"
echo "Smoke:      ${SMOKE_FLAG:-(none)}"
echo ""

for ARM in "${ARMS[@]}"; do
  RESULTS_DIR="${OUTPUT_ROOT}/results_${INTEGRATION}_${ARM}/${PHENOTYPE}"
  mkdir -p "${RESULTS_DIR}"
  jid=$(sbatch --parsable ${PART} -n 8 --mem=96G -t 8:00:00 \
    --export=ALL,LAUNCHER_SCRIPT_DIR="${SCRIPT_DIR}" \
    -o "${LOG_DIR}/%j_tfact_${ARM}.out" \
    -e "${LOG_DIR}/%j_tfact_${ARM}.err" \
    --job-name="tfact_${ARM}" \
    --wrap="set -e; source ${REPO_ROOT}/config/paths.sh; export HDF5_USE_FILE_LOCKING=FALSE; \
      \"\${DECOUPLER_ENV}/bin/python\" ${SCRIPT_DIR}/tf_activity_analysis.py \
        --integration ${INTEGRATION} --phenotype ${PHENOTYPE} \
        --input-dir ${INPUT_DIR} --pheno-csv \"\${ACE_SCORES_CSV}\" \
        --output-dir ${RESULTS_DIR} --run-progeny \
        --arm ${ARM} --sex-filter Male --celltypes ${CELLTYPES_CSV} ${SMOKE_FLAG}")
  echo "  ${ARM}: job ${jid}"
done

echo ""
echo "Monitor: squeue -u \$USER | grep tfact_"
