#!/bin/bash
# Submit Tsai ACE TF/Pathway Activity jobs using the shared config/output contract.
#
# Usage:
#   bash aceTfActT.sh [INTEGRATION]
#
# Arguments:
#   INTEGRATION  Integration method label (default: derived_batch)
#
# Submits one SLURM job per phenotype (tot_adverse_exp, early_hh_ses,
# ace_aggregate).  Each job processes ALL cell types internally, running
# DoRothEA, CollecTRI, and PROGENy.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

INTEGRATION="${1:-derived_batch}"
PHENOTYPES=("tot_adverse_exp" "early_hh_ses" "ace_aggregate")
INPUT_DIR="${ACE_OUTPUT_ROOT}/DEG/Tsai/celltype_splits_${INTEGRATION}"
OUTPUT_ROOT="${ACE_OUTPUT_ROOT}/TFActivity/Tsai"
LOG_DIR="${OUTPUT_ROOT}/logs"
mkdir -p "${LOG_DIR}"

echo "=== ACE TF Activity Pipeline: Tsai ==="
echo "Integration: ${INTEGRATION}"
echo "Phenotypes:  ${PHENOTYPES[*]}"
echo "Input dir:   ${INPUT_DIR}"
echo "Output root: ${OUTPUT_ROOT}"
echo ""

# Verify input directory exists
if [[ ! -d "${INPUT_DIR}" ]]; then
  echo "ERROR: Input directory does not exist: ${INPUT_DIR}"
  echo "       Run the DEG prep step first: bash Analysis/ACE/DEG/Tsai/aceDegT.sh"
  exit 1
fi

for PHENO in "${PHENOTYPES[@]}"; do
  RESULTS_DIR="${OUTPUT_ROOT}/results_${INTEGRATION}/${PHENO}"
  mkdir -p "${RESULTS_DIR}"

  JOB_ID=$(sbatch --parsable --export=ALL,REPO_ROOT="${REPO_ROOT}",LAUNCHER_SCRIPT_DIR="${SCRIPT_DIR}" \
    -o "${LOG_DIR}/%j_${PHENO}_${INTEGRATION}.out" \
    -e "${LOG_DIR}/%j_${PHENO}_${INTEGRATION}.err" \
    --job-name="ace_tf_${PHENO}" \
    "${SCRIPT_DIR}/run_tf_activity.sh" "${INTEGRATION}" "${PHENO}")

  echo "  ${PHENO}: job ${JOB_ID}"
done

echo ""
echo "All jobs submitted. Monitor with: squeue -u \$USER"
