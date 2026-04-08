#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

DEFAULT_ANALYSIS_ROOT="${ANALYSIS_OUTPUT_ROOT}"
SMOKE_ROOT="${DEFAULT_ANALYSIS_ROOT}/ACE/Smoke/Tsai/DEG"
FIXTURE_ROOT="${SMOKE_ROOT}/fixtures"
export ANALYSIS_OUTPUT_ROOT="${SMOKE_ROOT}/analysis_root"

activate_env "${NEBULA_ENV}"
export HDF5_USE_FILE_LOCKING=FALSE

python "${REPO_ROOT}/Analysis/ACE/_smoke/generate_fixtures.py" \
  --cohort tsai \
  --output-root "${FIXTURE_ROOT}"

export ACE_SCORES_CSV="${FIXTURE_ROOT}/ace_scores.csv"

python "${SCRIPT_DIR}/prep_celltype_splits.py" \
  --integration derived_batch \
  --annotated-h5ad "${FIXTURE_ROOT}/tsai/annotated/tsai_annotated.h5ad" \
  --doublet-dir "${FIXTURE_ROOT}/tsai/doublet_removed" \
  --celltype-filter broad_Exc \
  --overwrite

RESULTS_DIR="${ANALYSIS_OUTPUT_ROOT}/ACE/DEG/Tsai/results_derived_batch/tot_adverse_exp"
mkdir -p "${RESULTS_DIR}"

Rscript "${SCRIPT_DIR}/aceDegT.Rscript" \
  --integration derived_batch \
  --phenotype tot_adverse_exp \
  --input-dir "${ANALYSIS_OUTPUT_ROOT}/ACE/DEG/Tsai/celltype_splits_derived_batch" \
  --output-dir "${RESULTS_DIR}" \
  --pheno-csv "${ACE_SCORES_CSV}" \
  --celltype broad_Exc \
  --smoke

test -f "${RESULTS_DIR}/smoke_summary_tot_adverse_exp_Exc.csv"
echo "Tsai ACE DEG smoke test passed"
