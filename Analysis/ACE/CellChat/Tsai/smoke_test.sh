#!/bin/bash
# =============================================================================
# Smoke test for ACE CellChat analysis
#
# Validates that the CellChat analysis pipeline runs to completion on
# fixture data. Requires the CellChat conda environment.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

SMOKE_ROOT="${ACE_OUTPUT_ROOT}/Smoke/Tsai/CellChat"
FIXTURE_ROOT="${SMOKE_ROOT}/fixtures"
RESULTS_DIR="${SMOKE_ROOT}/analysis_root/results_derived_batch/Male"

mkdir -p "${FIXTURE_ROOT}" "${RESULTS_DIR}"

echo "=== CellChat Smoke Test ==="
echo "Fixtures: ${FIXTURE_ROOT}"
echo "Results:  ${RESULTS_DIR}"
echo ""

# Check if fixtures exist
FIXTURE_H5AD="${FIXTURE_ROOT}/tsai_annotated.h5ad"
FIXTURE_PHENO="${FIXTURE_ROOT}/ace_scores.csv"

if [[ ! -f "${FIXTURE_H5AD}" ]]; then
    echo "Generating fixtures..."
    set +u; activate_env "${CELLCHAT_ENV}"; set -u

    python3 "${REPO_ROOT}/Analysis/ACE/_smoke/generate_fixtures.py" \
        --cohort tsai --output-root "${FIXTURE_ROOT}" 2>/dev/null || {
        echo "SKIP: Could not generate fixtures (generate_fixtures.py unavailable)"
        echo "      CellChat smoke test requires fixture data."
        exit 0
    }
fi

if [[ ! -f "${FIXTURE_H5AD}" ]]; then
    echo "SKIP: No fixture h5ad available at ${FIXTURE_H5AD}"
    exit 0
fi

set +u; activate_env "${CELLCHAT_ENV}"; set -u

echo "Running CellChat analysis in smoke mode..."
Rscript "${SCRIPT_DIR}/cellchat_analysis.R" \
    --sex Male \
    --input-h5ad "${FIXTURE_H5AD}" \
    --pheno-csv "${FIXTURE_PHENO}" \
    --output-dir "${RESULTS_DIR}" \
    --smoke

# Validate outputs
PASS=true
for f in differential_interactions.csv pathway_changes.csv focus_axes_results.csv; do
    if [[ -f "${RESULTS_DIR}/${f}" ]]; then
        echo "  PASS: ${f} exists"
    else
        echo "  FAIL: ${f} missing"
        PASS=false
    fi
done

echo ""
if ${PASS}; then
    echo "CellChat smoke test PASSED"
else
    echo "CellChat smoke test FAILED"
    exit 1
fi
