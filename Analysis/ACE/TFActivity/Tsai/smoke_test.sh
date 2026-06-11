#!/bin/bash
# Smoke test for TF Activity analysis.
#
# Generates tiny fixtures, runs the analysis in --smoke mode (DoRothEA only,
# first cell type only), and verifies the output CSV was produced.
#
# Usage:
#   bash smoke_test.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

DEFAULT_ANALYSIS_ROOT="${ANALYSIS_OUTPUT_ROOT}"
SMOKE_ROOT="${DEFAULT_ANALYSIS_ROOT}/ACE/Smoke/Tsai/TFActivity"
FIXTURE_ROOT="${SMOKE_ROOT}/fixtures"
export ANALYSIS_OUTPUT_ROOT="${SMOKE_ROOT}/analysis_root"

activate_env "${DECOUPLER_ENV}"
export HDF5_USE_FILE_LOCKING=FALSE

# -----------------------------------------------------------------------
# Generate fixtures
# -----------------------------------------------------------------------

echo "Generating smoke-test fixtures..."
python "${REPO_ROOT}/Analysis/ACE/_smoke/generate_fixtures.py" \
  --cohort tsai \
  --output-root "${FIXTURE_ROOT}"

export ACE_SCORES_CSV="${FIXTURE_ROOT}/ace_scores.csv"

# -----------------------------------------------------------------------
# Prep cell-type splits (reuse DEG prep step)
# -----------------------------------------------------------------------

SPLIT_DIR="${ANALYSIS_OUTPUT_ROOT}/ACE/DEG/Tsai/celltype_splits_derived_batch"

# If DEG prep script is available, use it; otherwise create a minimal split
if [[ -f "${REPO_ROOT}/Analysis/ACE/DEG/Tsai/prep_celltype_splits.py" ]]; then
  echo "Running cell-type split preparation..."
  python "${REPO_ROOT}/Analysis/ACE/DEG/Tsai/prep_celltype_splits.py" \
    --integration derived_batch \
    --annotated-h5ad "${FIXTURE_ROOT}/tsai/annotated/tsai_annotated.h5ad" \
    --doublet-dir "${FIXTURE_ROOT}/tsai/doublet_removed" \
    --overwrite
else
  echo "DEG prep script not found; creating minimal fixture split..."
  mkdir -p "${SPLIT_DIR}"
  python -c "
import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

np.random.seed(42)
n_cells = 40
n_genes = 100

# Create fake gene names that overlap with DoRothEA targets
genes = [f'Gene{i}' for i in range(n_genes)]
patients = ['1001', '1002', '1003', '1004'] * 10

X = sp.csr_matrix(np.random.poisson(5, (n_cells, n_genes)).astype(np.float32))
obs = pd.DataFrame({
    'patient_id': patients,
    'cell_type': 'Ast',
}, index=[f'cell_{i}' for i in range(n_cells)])

adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes))
adata.write_h5ad('${SPLIT_DIR}/Ast.h5ad')
print(f'Wrote fixture: ${SPLIT_DIR}/Ast.h5ad  shape={adata.shape}')
"
fi

# -----------------------------------------------------------------------
# Run smoke test
# -----------------------------------------------------------------------

RESULTS_DIR="${ANALYSIS_OUTPUT_ROOT}/ACE/TFActivity/Tsai/results_derived_batch/tot_adverse_exp"
mkdir -p "${RESULTS_DIR}"

echo "Running TF activity analysis in smoke mode..."
python "${SCRIPT_DIR}/tf_activity_analysis.py" \
  --integration derived_batch \
  --phenotype tot_adverse_exp \
  --input-dir "${SPLIT_DIR}" \
  --pheno-csv "${ACE_SCORES_CSV}" \
  --output-dir "${RESULTS_DIR}" \
  --smoke

# -----------------------------------------------------------------------
# Verify outputs
# -----------------------------------------------------------------------

echo ""
echo "Checking outputs..."

# In smoke mode, the script writes smoke_summary_<phenotype>.csv
SMOKE_SUMMARY="${RESULTS_DIR}/smoke_summary_tot_adverse_exp.csv"
if [[ -f "${SMOKE_SUMMARY}" ]]; then
  echo "  PASS: ${SMOKE_SUMMARY} exists"
  ROWS=$(wc -l < "${SMOKE_SUMMARY}")
  echo "  INFO: ${ROWS} lines (including header)"
else
  echo "  FAIL: ${SMOKE_SUMMARY} not found"
  exit 1
fi

echo ""
echo "Tsai ACE TF Activity smoke test passed."
