#!/bin/bash
# =============================================================================
# Smoke test for ACE hdWGCNA analysis
#
# Validates that the hdWGCNA pipeline runs to completion on fixture data.
# Requires the WGCNA conda environment.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

SMOKE_ROOT="${ACE_OUTPUT_ROOT}/Smoke/Tsai/hdWGCNA"
FIXTURE_ROOT="${SMOKE_ROOT}/fixtures"
RESULTS_DIR="${SMOKE_ROOT}/analysis_root/results_derived_batch/Male_Exc"

mkdir -p "${FIXTURE_ROOT}" "${RESULTS_DIR}"

echo "=== hdWGCNA Smoke Test ==="

# Check for fixtures
FIXTURE_H5AD="${FIXTURE_ROOT}/Exc.h5ad"
FIXTURE_PHENO="${FIXTURE_ROOT}/ace_scores.csv"

if [[ ! -f "${FIXTURE_H5AD}" ]]; then
    echo "Generating fixtures..."
    set +u; activate_env "${WGCNA_ENV}"; set -u

    python3 "${REPO_ROOT}/Analysis/ACE/_smoke/generate_fixtures.py" \
        --cohort tsai --output-root "${FIXTURE_ROOT}" 2>/dev/null || {
        echo "SKIP: Could not generate fixtures"
        exit 0
    }
fi

if [[ ! -f "${FIXTURE_H5AD}" ]]; then
    echo "SKIP: No fixture h5ad available"
    exit 0
fi

set +u; activate_env "${WGCNA_ENV}"; set -u

echo "Running hdWGCNA analysis in smoke mode..."
Rscript "${SCRIPT_DIR}/wgcna_analysis.R" \
    --cell-type Exc \
    --sex Male \
    --phenotype tot_adverse_exp \
    --input-h5ad "${FIXTURE_H5AD}" \
    --pheno-csv "${FIXTURE_PHENO}" \
    --output-dir "${RESULTS_DIR}" \
    --smoke

PASS=true
for f in module_assignments.csv module_eigengenes.csv module_trait_correlations.csv; do
    if [[ -f "${RESULTS_DIR}/${f}" ]]; then
        echo "  PASS: ${f} exists"
    else
        echo "  FAIL: ${f} missing"
        PASS=false
    fi
done

echo ""
if ${PASS}; then
    echo "hdWGCNA smoke test PASSED"
else
    echo "hdWGCNA smoke test FAILED"
    exit 1
fi
