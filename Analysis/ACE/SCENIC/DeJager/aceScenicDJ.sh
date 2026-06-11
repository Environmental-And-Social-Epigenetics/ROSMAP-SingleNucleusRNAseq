#!/bin/bash
# Submit DeJager ACE SCENIC jobs using the shared config/output contract.
#
# Mirrors Analysis/ACE/SCENIC/Tsai/aceScenicT.sh, but points at the DeJager DEG
# celltype splits and the DeJager SCENIC output dir, and runs the SAME modular
# pipeline (scenic_analysis.py) with --cohort dejager and the SAME pool_size
# (default 50). Submits one independent SLURM job per cell_type x sex.
#
# Usage:
#   bash aceScenicDJ.sh [INTEGRATION] [PHENOTYPE]
#
# Defaults:
#   INTEGRATION = library_id    (DeJager canonical integration; see DEG/DeJager)
#   PHENOTYPE   = tot_adverse_exp

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
INTEGRATION="${1:-library_id}"
PHENOTYPE="${2:-tot_adverse_exp}"

# ---------------------------------------------------------------------------
# Validate SCENIC ranking databases
# ---------------------------------------------------------------------------
if [[ "${SCENIC_RANKING_DIR}" == *"__UNCONFIGURED__"* ]]; then
  echo "ERROR: SCENIC_RANKING_DIR is not configured."
  echo "       Set it in config/paths.local.sh. See config/paths.sh for details."
  exit 1
fi

if [[ ! -d "${SCENIC_RANKING_DIR}" ]]; then
  echo "ERROR: SCENIC_RANKING_DIR does not exist: ${SCENIC_RANKING_DIR}"
  exit 1
fi

# ---------------------------------------------------------------------------
# Paths (DeJager DEG splits in, DeJager SCENIC results out)
# ---------------------------------------------------------------------------
INPUT_DIR="${ACE_OUTPUT_ROOT}/DEG/DeJager/celltype_splits_${INTEGRATION}"
OUTPUT_ROOT="${ACE_OUTPUT_ROOT}/SCENIC/DeJager/results_${INTEGRATION}/${PHENOTYPE}"
LOG_DIR="${ACE_OUTPUT_ROOT}/SCENIC/DeJager/logs"
mkdir -p "${LOG_DIR}"

# ---------------------------------------------------------------------------
# Cell types per sex
# ---------------------------------------------------------------------------
# Use the same broad_Exc/broad_Inh split stems produced by the DEG prep step
# so jobs are not silently skipped when a split file is missing.
MALE_CELL_TYPES=("broad_Inh" "Mic" "Ast" "In-PV_Basket" "Oli" "broad_Exc")
FEMALE_CELL_TYPES=("Oli" "Ast" "Ex-L2_3" "broad_Exc" "OPC" "broad_Inh" "Mic")

# ---------------------------------------------------------------------------
# Submit
# ---------------------------------------------------------------------------
echo "=== ACE SCENIC Pipeline: DeJager ==="
echo "Integration: ${INTEGRATION}"
echo "Phenotype:   ${PHENOTYPE}"
echo "Input dir:   ${INPUT_DIR}"
echo "Output root: ${OUTPUT_ROOT}"
echo "Method:      shared modular scenic_analysis.py (--cohort dejager, pool_size=50)"
echo ""

submit_job() {
  local sex="$1"
  local cell_type="$2"

  local h5ad="${INPUT_DIR}/${cell_type}.h5ad"
  if [[ ! -f "${h5ad}" ]]; then
    echo "  WARNING: Input not found, skipping: ${h5ad}"
    return
  fi

  local out_dir="${OUTPUT_ROOT}/${sex}_${cell_type}"
  mkdir -p "${out_dir}"

  local job_name="scenicDJ_${sex}_${cell_type}"

  local job_id
  job_id=$(sbatch --parsable --export=ALL,REPO_ROOT="${REPO_ROOT}",LAUNCHER_SCRIPT_DIR="${SCRIPT_DIR}" \
    -o "${LOG_DIR}/%j_${sex}_${cell_type}.out" \
    -e "${LOG_DIR}/%j_${sex}_${cell_type}.err" \
    --job-name="${job_name}" \
    "${SCRIPT_DIR}/run_scenic.sh" \
    "${INTEGRATION}" "${PHENOTYPE}" "${cell_type}" "${sex}")

  echo "  Submitted ${sex}/${cell_type}: job ${job_id}"
}

echo "--- Male cell types ---"
for ct in "${MALE_CELL_TYPES[@]}"; do
  submit_job "Male" "${ct}"
done

echo ""
echo "--- Female cell types ---"
for ct in "${FEMALE_CELL_TYPES[@]}"; do
  submit_job "Female" "${ct}"
done

echo ""
echo "All jobs submitted. Monitor with: squeue -u \$USER"
echo "Results will appear in: ${OUTPUT_ROOT}/"
