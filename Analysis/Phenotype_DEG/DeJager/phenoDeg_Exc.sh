#!/bin/bash
# phenoDeg_Exc.sh -- pseudobulk Exc in Python (avoids R dgCMatrix limit),
# then run DESeq2. Single (integration × phenotype) job; the Python
# pseudobulk only needs to run once per integration but we re-run it
# per-phenotype for simplicity (it's fast: ~20 min).
#
# Usage:  sbatch phenoDeg_Exc.sh <integration> <phenotype>

#SBATCH -p ou_bcs_normal,mit_normal
#SBATCH -n 8
#SBATCH --mem=128G
#SBATCH -t 06:00:00
#SBATCH -o %j.out
#SBATCH -e %j.err

set -euo pipefail

if [[ -z "${REPO_ROOT:-}" ]]; then
  if git -C "${PWD}" rev-parse --show-toplevel >/dev/null 2>&1; then
    REPO_ROOT="$(git -C "${PWD}" rev-parse --show-toplevel)"
  else
    echo "ERROR: REPO_ROOT is not set and could not be inferred from the current directory."
    exit 1
  fi
fi
SCRIPT_DIR="${REPO_ROOT}/Analysis/Phenotype_DEG/DeJager"
source "${REPO_ROOT}/config/paths.sh"

INTEGRATION="${1:?integration required}"
PHENOTYPE="${2:?phenotype required (msex|cogdx_binary)}"

OUTPUT_ROOT="${ANALYSIS_OUTPUT_ROOT:-${WORKSPACE_ROOT}/Analysis_Outputs}/Phenotype_DEG/DeJager"
SPLIT_DIR="${OUTPUT_ROOT}/celltype_splits_${INTEGRATION}"
PSEUDO_DIR="${OUTPUT_ROOT}/results_${INTEGRATION}/Exc_pseudobulk"
RESULTS_DIR="${OUTPUT_ROOT}/results_${INTEGRATION}/${PHENOTYPE}"
mkdir -p "${PSEUDO_DIR}" "${RESULTS_DIR}"

PHENO_CSV="${ACE_SCORES_CSV:-${REPO_ROOT}/Data/Phenotypes/TSAI_DEJAGER_all_patients_wACEscores.csv}"

set +u
source "${CONDA_INIT_SCRIPT}"
conda activate "${NEBULA_ENV}"
set -u
export HDF5_USE_FILE_LOCKING=FALSE
export PYTHONUNBUFFERED=1

echo "=== Exc pseudobulk + DEG: ${INTEGRATION} × ${PHENOTYPE} ==="

# Step 1: Python pseudobulk (idempotent — skip if csv exists)
if [[ -s "${PSEUDO_DIR}/pseudobulk_counts_Exc.csv" ]]; then
  echo "Pseudobulk CSVs already present in ${PSEUDO_DIR}, skipping Python step."
else
  python -u "${SCRIPT_DIR}/pseudobulk_Exc_python.py" \
    --integration "${INTEGRATION}" \
    --input-dir "${SPLIT_DIR}" \
    --output-dir "${PSEUDO_DIR}"
fi

# Step 2: R DESeq2
Rscript "${SCRIPT_DIR}/phenoDeg_Exc.Rscript" \
  --integration "${INTEGRATION}" \
  --phenotype "${PHENOTYPE}" \
  --input-dir "${PSEUDO_DIR}" \
  --output-dir "${RESULTS_DIR}" \
  --pheno-csv "${PHENO_CSV}"

echo "Done."
