#!/bin/bash
# Submit Tsai Resilient DEG jobs using the shared config/output contract.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

INTEGRATIONS=("derived_batch" "projid")
PHENOTYPES=("resilience_group")
OUTPUT_ROOT="${RESILIENT_OUTPUT_ROOT}/DEG/Tsai"
LOG_DIR="${OUTPUT_ROOT}/logs"
mkdir -p "${LOG_DIR}"

CONDA_SETUP="source ${CONDA_INIT_SCRIPT} && conda activate ${NEBULA_ENV} && export HDF5_USE_FILE_LOCKING=FALSE"

echo "=== Resilient DEG Pipeline: Tsai ==="
echo "Integrations: ${INTEGRATIONS[*]}"
echo "Phenotypes:   ${PHENOTYPES[*]}"
echo "Output root:  ${OUTPUT_ROOT}"
echo ""

for INTEGRATION in "${INTEGRATIONS[@]}"; do
  SPLIT_DIR="${OUTPUT_ROOT}/celltype_splits_${INTEGRATION}"
  PREP_JOB=$(sbatch --parsable \
    -p pi_lhtsai,pi_manoli \
    -n 16 \
    --mem=200G \
    -t 24:00:00 \
    -o "${LOG_DIR}/%j_prep_${INTEGRATION}.out" \
    -e "${LOG_DIR}/%j_prep_${INTEGRATION}.err" \
    --job-name="resil_prep_${INTEGRATION}" \
    --wrap="${CONDA_SETUP} && python ${SCRIPT_DIR}/prep_celltype_splits.py --integration ${INTEGRATION} --output-dir ${SPLIT_DIR}")

  echo "  Preprocessing ${INTEGRATION}: job ${PREP_JOB}"

  for PHENO in "${PHENOTYPES[@]}"; do
    DEG_JOB=$(sbatch --parsable \
      --dependency=afterok:${PREP_JOB} \
      -o "${LOG_DIR}/%j_${PHENO}_${INTEGRATION}.out" \
      -e "${LOG_DIR}/%j_${PHENO}_${INTEGRATION}.err" \
      --job-name="resil_${PHENO}_${INTEGRATION}" \
      "${SCRIPT_DIR}/run_deg.sh" "${INTEGRATION}" "${PHENO}")

    echo "  DEG ${PHENO}/${INTEGRATION}: job ${DEG_JOB}"
  done
done

echo ""
echo "All jobs submitted. Monitor with: squeue -u \$USER"
