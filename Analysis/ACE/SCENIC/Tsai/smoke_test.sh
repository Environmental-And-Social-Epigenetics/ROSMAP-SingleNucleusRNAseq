#!/bin/bash
# Smoke test for ACE SCENIC Tsai pipeline.
#
# Generates a tiny synthetic h5ad and phenotype CSV, then runs
# scenic_analysis.py in --smoke mode to verify the pipeline does not crash.
#
# Usage:
#   bash smoke_test.sh
#
# Requirements:
#   - SCENIC_ANALYSIS_ENV configured in config/paths.local.sh
#   - SCENIC_RANKING_DIR configured and databases present

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

# ---------------------------------------------------------------------------
# Validate prerequisites
# ---------------------------------------------------------------------------
if [[ "${SCENIC_RANKING_DIR}" == *"__UNCONFIGURED__"* ]]; then
  echo "SKIP: SCENIC_RANKING_DIR is not configured. Cannot run smoke test."
  exit 0
fi

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
DEFAULT_ANALYSIS_ROOT="${ANALYSIS_OUTPUT_ROOT}"
SMOKE_ROOT="${DEFAULT_ANALYSIS_ROOT}/ACE/Smoke/Tsai/SCENIC"
FIXTURE_DIR="${SMOKE_ROOT}/fixtures"
SMOKE_OUTPUT="${SMOKE_ROOT}/output"

rm -rf "${SMOKE_ROOT}"
mkdir -p "${FIXTURE_DIR}" "${SMOKE_OUTPUT}"
export FIXTURE_DIR

set +u
activate_env "${SCENIC_ANALYSIS_ENV}"
set -u
export HDF5_USE_FILE_LOCKING=FALSE

# ---------------------------------------------------------------------------
# Generate synthetic fixture data
# ---------------------------------------------------------------------------
echo "Generating synthetic test data ..."

"${SCENIC_ANALYSIS_ENV}/bin/python" - <<'PYEOF'
import anndata as ad
import numpy as np
import pandas as pd
import os
import scipy.sparse as sp

np.random.seed(42)

fixture_dir = os.environ["FIXTURE_DIR"]

n_patients = 15
n_cells_per_patient = 100
n_genes = 600

patient_ids = [str(1000 + i) for i in range(n_patients)]

obs_rows = []
for pid in patient_ids:
    for j in range(n_cells_per_patient):
        obs_rows.append({"patient_id": pid, "cell_type": "Ast"})
obs = pd.DataFrame(obs_rows)
obs.index = [f"cell_{i}" for i in range(len(obs))]

gene_names = [f"Gene{i}" for i in range(n_genes)]
X = sp.random(len(obs), n_genes, density=0.1, format="csr", dtype=np.float32)
X.data = np.round(X.data * 100).astype(np.float32)

adata = ad.AnnData(X=X, obs=obs)
adata.var_names = gene_names

h5ad_path = os.path.join(fixture_dir, "Ast.h5ad")
adata.write_h5ad(h5ad_path)
print(f"Wrote {h5ad_path}: {adata.shape}")

# Phenotype CSV
pheno = pd.DataFrame({
    "projid": patient_ids,
    "tot_adverse_exp": np.random.randint(0, 5, n_patients),
    "early_hh_ses": np.random.randint(0, 3, n_patients),
    "ace_aggregate": np.random.uniform(0, 4, n_patients),
    "msex": [1] * 8 + [0] * 7,
    "age_death": np.random.uniform(70, 95, n_patients),
    "pmi": np.random.uniform(2, 12, n_patients),
    "niareagansc": np.random.choice([1, 2, 3, 4], n_patients),
})
pheno_path = os.path.join(fixture_dir, "ace_scores.csv")
pheno.to_csv(pheno_path, index=False)
print(f"Wrote {pheno_path}")
PYEOF

# ---------------------------------------------------------------------------
# Run smoke test
# ---------------------------------------------------------------------------
echo ""
echo "Running scenic_analysis.py in smoke mode ..."

"${SCENIC_ANALYSIS_ENV}/bin/python" "${SCRIPT_DIR}/scenic_analysis.py" \
  --cell-type Ast \
  --sex Male \
  --phenotype tot_adverse_exp \
  --input-h5ad "${FIXTURE_DIR}/Ast.h5ad" \
  --pheno-csv "${FIXTURE_DIR}/ace_scores.csv" \
  --output-dir "${SMOKE_OUTPUT}/Male_Ast" \
  --ranking-dir "${SCENIC_RANKING_DIR}" \
  --tf-list "${SCENIC_TF_LIST}" \
  --pool-size 10 \
  --num-workers 2 \
  --smoke

# ---------------------------------------------------------------------------
# Validate outputs
# ---------------------------------------------------------------------------
echo ""
echo "Validating outputs ..."

EXPECTED_FILES=(
  "adjacencies.csv.gz"
  "regulons.pkl"
  "regulons_list.csv"
  "auc_matrix.csv"
  "regression_results.csv"
)

PASS=true
for f in "${EXPECTED_FILES[@]}"; do
  fpath="${SMOKE_OUTPUT}/Male_Ast/${f}"
  if [[ -f "${fpath}" ]]; then
    echo "  OK: ${f}"
  else
    echo "  FAIL: ${f} not found"
    PASS=false
  fi
done

echo ""
if ${PASS}; then
  echo "Tsai ACE SCENIC smoke test PASSED"
else
  echo "Tsai ACE SCENIC smoke test FAILED"
  exit 1
fi
