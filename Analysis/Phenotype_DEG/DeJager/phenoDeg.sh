#!/bin/bash
# phenoDeg.sh -- Phase D2 wrapper: run a single (integration × phenotype × celltype)
# DESeq2 fit using phenoDeg.Rscript.
#
# Usage:
#   sbatch phenoDeg.sh <integration> <phenotype> <celltype>
#   e.g.  sbatch phenoDeg.sh library_id msex Mic

#SBATCH -p ou_bcs_normal,mit_normal
#SBATCH -n 4
#SBATCH --mem=64G
#SBATCH -t 04:00:00
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

INTEGRATION="${1:?integration required (library_id|patient_id|pool_batch|derived_batch|sequencing_date)}"
PHENOTYPE="${2:?phenotype required (msex|cogdx_binary)}"
CELLTYPE="${3:?celltype required (Exc|Inh|Ast|Mic|Oli|OPC)}"

OUTPUT_ROOT="${ANALYSIS_OUTPUT_ROOT:-${WORKSPACE_ROOT}/Analysis_Outputs}/Phenotype_DEG/DeJager"
INPUT_DIR="${OUTPUT_ROOT}/celltype_splits_${INTEGRATION}"
RESULTS_DIR="${OUTPUT_ROOT}/results_${INTEGRATION}/${PHENOTYPE}"
mkdir -p "${RESULTS_DIR}"

# Use the ACE phenotype CSV — it already contains msex, cogdx, age_death, pmi
PHENO_CSV="${ACE_SCORES_CSV:-${REPO_ROOT}/Data/Phenotypes/TSAI_DEJAGER_all_patients_wACEscores.csv}"

set +u
source "${CONDA_INIT_SCRIPT}"
conda activate "${NEBULA_ENV}"
set -u
export HDF5_USE_FILE_LOCKING=FALSE

echo "=== Phenotype DEG: ${INTEGRATION} × ${PHENOTYPE} × ${CELLTYPE} ==="
echo "Input dir:  ${INPUT_DIR}"
echo "Results:    ${RESULTS_DIR}"
echo "Pheno CSV:  ${PHENO_CSV}"

Rscript "${SCRIPT_DIR}/phenoDeg.Rscript" \
  --integration "${INTEGRATION}" \
  --phenotype "${PHENOTYPE}" \
  --celltype "${CELLTYPE}" \
  --input-dir "${INPUT_DIR}" \
  --output-dir "${RESULTS_DIR}" \
  --pheno-csv "${PHENO_CSV}"

echo "Done."
