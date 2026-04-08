#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

INTEGRATIONS=("library_id")
PHENOTYPES=("tot_adverse_exp" "early_hh_ses" "ace_aggregate")
OUTPUT_ROOT="${ANALYSIS_OUTPUT_ROOT}/ACE/DEG/DeJager"
LOG_DIR="${OUTPUT_ROOT}/logs"
mkdir -p "${LOG_DIR}"

CONDA_SETUP="source ${CONDA_INIT_SCRIPT} && conda activate ${NEBULA_ENV} && export HDF5_USE_FILE_LOCKING=FALSE"

echo "=== ACE DEG Pipeline: DeJager ==="
echo "Integrations: ${INTEGRATIONS[*]}"
echo "Phenotypes:   ${PHENOTYPES[*]}"

for INTEGRATION in "${INTEGRATIONS[@]}"; do
  SPLIT_DIR="${OUTPUT_ROOT}/celltype_splits_${INTEGRATION}"
  PREP_JOB=$(sbatch --parsable \
    -p pi_lhtsai,pi_tpoggio,pi_manoli \
    -n 16 \
    --mem=200G \
    -t 24:00:00 \
    -o "${LOG_DIR}/%j_prep_${INTEGRATION}.out" \
    -e "${LOG_DIR}/%j_prep_${INTEGRATION}.err" \
    --job-name="aceDJ_prep_${INTEGRATION}" \
    --wrap="${CONDA_SETUP} && python ${SCRIPT_DIR}/prep_celltype_splits.py --integration ${INTEGRATION} --output-dir ${SPLIT_DIR}")

  echo "  Preprocessing ${INTEGRATION}: job ${PREP_JOB}"

  for PHENO in "${PHENOTYPES[@]}"; do
    DEG_JOB=$(sbatch --parsable \
      --dependency=afterok:${PREP_JOB} \
      -o "${LOG_DIR}/%j_${PHENO}_${INTEGRATION}.out" \
      -e "${LOG_DIR}/%j_${PHENO}_${INTEGRATION}.err" \
      --job-name="aceDJ_${PHENO}_${INTEGRATION}" \
      "${SCRIPT_DIR}/run_deg.sh" "${INTEGRATION}" "${PHENO}")
    echo "  DEG ${PHENO}/${INTEGRATION}: job ${DEG_JOB}"
  done
done

echo "All jobs submitted. Monitor with: squeue -u \$USER"
