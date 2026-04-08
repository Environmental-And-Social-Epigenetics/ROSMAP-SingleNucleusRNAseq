#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

DEFAULT_ANALYSIS_ROOT="${ANALYSIS_OUTPUT_ROOT}"
SMOKE_ROOT="${DEFAULT_ANALYSIS_ROOT}/ACE/Smoke/DeJager/CellTypeProportion"
FIXTURE_ROOT="${SMOKE_ROOT}/fixtures"
export ANALYSIS_OUTPUT_ROOT="${SMOKE_ROOT}/analysis_root"

activate_env "${NEBULA_ENV}"
export HDF5_USE_FILE_LOCKING=FALSE

python "${REPO_ROOT}/Analysis/ACE/_smoke/generate_fixtures.py" \
  --cohort dejager \
  --output-root "${FIXTURE_ROOT}"

export ACE_SCORES_CSV="${FIXTURE_ROOT}/ace_scores.csv"

python "${SCRIPT_DIR}/prep_counts.py" \
  --integration library_id \
  --h5ad "${FIXTURE_ROOT}/dejager/annotated/dejager_annotated.h5ad" \
  --pheno-csv "${ACE_SCORES_CSV}" \
  --output-root "${ANALYSIS_OUTPUT_ROOT}/ACE/CellTypeProportion/DeJager"

activate_env "${SCCOMP_ENV}"

Rscript "${SCRIPT_DIR}/sccomp_analysis.R" \
  --integration library_id \
  --sex all \
  --resolution fine \
  --output-root "${ANALYSIS_OUTPUT_ROOT}/ACE/CellTypeProportion/DeJager" \
  --smoke

test -f "${ANALYSIS_OUTPUT_ROOT}/ACE/CellTypeProportion/DeJager/results_library_id/smoke/sccomp_smoke_summary.csv"
echo "DeJager ACE cell-type proportion smoke test passed"
