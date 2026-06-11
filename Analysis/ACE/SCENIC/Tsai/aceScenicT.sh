#!/bin/bash
# Submit Tsai ACE SCENIC jobs using the shared config/output contract.
#
# Submits one SLURM job per cell_type x sex combination.
# Each job is independent (no dependencies).
#
# Usage:
#   bash aceScenicT.sh [INTEGRATION] [PHENOTYPE]
#
# Defaults:
#   INTEGRATION = derived_batch
#   PHENOTYPE   = tot_adverse_exp

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
INTEGRATION="${1:-derived_batch}"
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
# Paths
# ---------------------------------------------------------------------------
INPUT_DIR="${ACE_OUTPUT_ROOT}/DEG/Tsai/celltype_splits_${INTEGRATION}"
OUTPUT_ROOT="${ACE_OUTPUT_ROOT}/SCENIC/Tsai/results_${INTEGRATION}/${PHENOTYPE}"
LOG_DIR="${ACE_OUTPUT_ROOT}/SCENIC/Tsai/logs"
mkdir -p "${LOG_DIR}"

# ---------------------------------------------------------------------------
# Cell types per sex
# ---------------------------------------------------------------------------
# Priority cell types (order determines submission order, not priority).
# Broad excitatory/inhibitory split files are named broad_Exc/broad_Inh by
# the DEG prep step; use those exact stems so jobs are not silently skipped.
MALE_CELL_TYPES=("broad_Inh" "Mic" "Ast" "In-PV_Basket" "Oli" "broad_Exc")
FEMALE_CELL_TYPES=("Oli" "Ast" "Ex-L2_3" "broad_Exc" "OPC" "broad_Inh" "Mic")

# ---------------------------------------------------------------------------
# Submit
# ---------------------------------------------------------------------------
echo "=== ACE SCENIC Pipeline: Tsai ==="
echo "Integration: ${INTEGRATION}"
echo "Phenotype:   ${PHENOTYPE}"
echo "Input dir:   ${INPUT_DIR}"
echo "Output root: ${OUTPUT_ROOT}"
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

  local job_name="scenic_${sex}_${cell_type}"

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
