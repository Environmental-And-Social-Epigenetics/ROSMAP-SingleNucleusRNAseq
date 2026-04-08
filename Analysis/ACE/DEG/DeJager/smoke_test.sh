#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

DEFAULT_ANALYSIS_ROOT="${ANALYSIS_OUTPUT_ROOT}"
SMOKE_ROOT="${DEFAULT_ANALYSIS_ROOT}/ACE/Smoke/DeJager/DEG"
FIXTURE_ROOT="${SMOKE_ROOT}/fixtures"
export ANALYSIS_OUTPUT_ROOT="${SMOKE_ROOT}/analysis_root"

activate_env "${NEBULA_ENV}"
export HDF5_USE_FILE_LOCKING=FALSE

python "${REPO_ROOT}/Analysis/ACE/_smoke/generate_fixtures.py" \
  --cohort dejager \
  --output-root "${FIXTURE_ROOT}"

export ACE_SCORES_CSV="${FIXTURE_ROOT}/ace_scores.csv"

python "${SCRIPT_DIR}/prep_celltype_splits.py" \
  --integration library_id \
  --annotated-h5ad "${FIXTURE_ROOT}/dejager/annotated/dejager_annotated.h5ad" \
  --doublet-dir "${FIXTURE_ROOT}/dejager/doublet_removed" \
  --celltype-filter broad_Exc \
  --overwrite

RESULTS_DIR="${ANALYSIS_OUTPUT_ROOT}/ACE/DEG/DeJager/results_library_id/tot_adverse_exp"
mkdir -p "${RESULTS_DIR}"

Rscript "${SCRIPT_DIR}/aceDegDJ.Rscript" \
  --integration library_id \
  --phenotype tot_adverse_exp \
  --input-dir "${ANALYSIS_OUTPUT_ROOT}/ACE/DEG/DeJager/celltype_splits_library_id" \
  --output-dir "${RESULTS_DIR}" \
  --pheno-csv "${ACE_SCORES_CSV}" \
  --celltype broad_Exc \
  --smoke

test -f "${RESULTS_DIR}/smoke_summary_tot_adverse_exp_Exc.csv"
echo "DeJager ACE DEG smoke test passed"
