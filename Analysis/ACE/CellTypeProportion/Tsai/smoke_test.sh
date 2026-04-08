#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

DEFAULT_ANALYSIS_ROOT="${ANALYSIS_OUTPUT_ROOT}"
SMOKE_ROOT="${DEFAULT_ANALYSIS_ROOT}/ACE/Smoke/Tsai/CellTypeProportion"
FIXTURE_ROOT="${SMOKE_ROOT}/fixtures"
export ANALYSIS_OUTPUT_ROOT="${SMOKE_ROOT}/analysis_root"

activate_env "${NEBULA_ENV}"
export HDF5_USE_FILE_LOCKING=FALSE

python "${REPO_ROOT}/Analysis/ACE/_smoke/generate_fixtures.py" \
  --cohort tsai \
  --output-root "${FIXTURE_ROOT}"

export ACE_SCORES_CSV="${FIXTURE_ROOT}/ace_scores.csv"

python "${SCRIPT_DIR}/prep_counts.py" \
  --integration derived_batch \
  --h5ad "${FIXTURE_ROOT}/tsai/annotated/tsai_annotated.h5ad" \
  --pheno-csv "${ACE_SCORES_CSV}" \
  --output-root "${ANALYSIS_OUTPUT_ROOT}/ACE/CellTypeProportion/Tsai"

activate_env "${SCCOMP_ENV}"

Rscript "${SCRIPT_DIR}/sccomp_analysis.R" \
  --integration derived_batch \
  --sex all \
  --resolution fine \
  --output-root "${ANALYSIS_OUTPUT_ROOT}/ACE/CellTypeProportion/Tsai" \
  --smoke

test -f "${ANALYSIS_OUTPUT_ROOT}/ACE/CellTypeProportion/Tsai/results_derived_batch/smoke/sccomp_smoke_summary.csv"
echo "Tsai ACE cell-type proportion smoke test passed"
