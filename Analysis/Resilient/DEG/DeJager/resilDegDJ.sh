#!/bin/bash
# Submit DeJager Resilient DEG jobs.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

OUTPUT_ROOT="${RESILIENT_OUTPUT_ROOT}/DEG/DeJager"
LOG_DIR="${OUTPUT_ROOT}/logs"
mkdir -p "${LOG_DIR}"

CONDA_SETUP="source ${CONDA_INIT_SCRIPT} && conda activate ${NEBULA_ENV} && export HDF5_USE_FILE_LOCKING=FALSE"

echo "=== Resilient DEG Pipeline: DeJager ==="

INTEGRATION="library_id"
SPLIT_DIR="${OUTPUT_ROOT}/celltype_splits_${INTEGRATION}"

PREP_JOB=$(sbatch --parsable \
  -p pi_lhtsai,pi_manoli \
  -n 16 \
  --mem=200G \
  -t 24:00:00 \
  -o "${LOG_DIR}/%j_prep_${INTEGRATION}.out" \
  -e "${LOG_DIR}/%j_prep_${INTEGRATION}.err" \
  --job-name="resil_dj_prep" \
  --wrap="${CONDA_SETUP} && python ${SCRIPT_DIR}/prep_celltype_splits.py --integration ${INTEGRATION} --output-dir ${SPLIT_DIR}")

echo "  Preprocessing: job ${PREP_JOB}"

for PHENO in resilience_group; do
  RESULTS_DIR="${OUTPUT_ROOT}/results_${INTEGRATION}/${PHENO}"
  mkdir -p "${RESULTS_DIR}"

  DEG_JOB=$(sbatch --parsable \
    --dependency=afterok:${PREP_JOB} \
    -p pi_lhtsai,pi_manoli \
    -n 4 \
    --mem=100G \
    -t 12:00:00 \
    -o "${LOG_DIR}/%j_${PHENO}_${INTEGRATION}.out" \
    -e "${LOG_DIR}/%j_${PHENO}_${INTEGRATION}.err" \
    --job-name="resil_dj_${PHENO}" \
    --wrap="${CONDA_SETUP} && Rscript ${SCRIPT_DIR}/resilDegDJ.Rscript --integration ${INTEGRATION} --phenotype ${PHENO} --input-dir ${SPLIT_DIR} --output-dir ${RESULTS_DIR} --pheno-csv ${RESILIENT_PHENOTYPE_CSV}")

  echo "  DEG ${PHENO}: job ${DEG_JOB}"
done

echo ""
echo "All jobs submitted."
