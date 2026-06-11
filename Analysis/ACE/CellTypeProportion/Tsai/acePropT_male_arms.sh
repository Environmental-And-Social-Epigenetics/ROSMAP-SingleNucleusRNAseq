#!/bin/bash
# Cell-type proportion (sccomp) for the non-ANCOVA male ACE AD-model arms.
# One SLURM job per arm: males-only, fine resolution, single primary
# tot_adverse_exp composition model whose covariates match the DEG arm:
#   ~ tot_adverse_exp + age_death + pmi [+ arm AD covars] [+ AD_binary:tot_adverse_exp]
# Reuses the already-prepped count + metadata CSVs (no prep re-run). Output:
#   ${ACE_PROP_OUTPUT_ROOT}/results_derived_batch_<ARM>/primary/
#
# Usage:
#   bash acePropT_male_arms.sh
#   SMOKE_FLAG=--smoke bash acePropT_male_arms.sh
#   ARMS_OVERRIDE="MaleContAD" bash acePropT_male_arms.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"
source "${REPO_ROOT}/Analysis/ACE/_shared/arm_covariates.sh"

INTEGRATION="derived_batch"
RESOLUTION="fine"
OUTPUT_ROOT="${ACE_PROP_OUTPUT_ROOT:-${ACE_OUTPUT_ROOT}/CellTypeProportion/Tsai}"
LOG_DIR="${OUTPUT_ROOT}/logs"
mkdir -p "${LOG_DIR}"

PART="-p pi_lhtsai,pi_manoli"
SMOKE_FLAG="${SMOKE_FLAG:-}"

read -r -a ARMS <<< "${ARMS_OVERRIDE:-${NON_ANCOVA_ARMS[*]}}"

# Verify prepped inputs exist (prep is phenotype-independent and already done).
COUNTS="${OUTPUT_ROOT}/data/cell_counts_${RESOLUTION}_${INTEGRATION}.csv"
META="${OUTPUT_ROOT}/data/metadata_${INTEGRATION}.csv"
if [[ ! -s "${COUNTS}" || ! -s "${META}" ]]; then
  echo "ERROR: prepped inputs missing (${COUNTS} / ${META}). Run run_prep.sh first." >&2
  exit 1
fi

echo "=== ACE cell-type proportion (sccomp): male AD-model arms ==="
echo "Arms:       ${ARMS[*]}"
echo "Resolution: ${RESOLUTION} (males, primary tot_adverse_exp composition)"
echo "Smoke:      ${SMOKE_FLAG:-(none)}"
echo ""

for ARM in "${ARMS[@]}"; do
  jid=$(sbatch --parsable ${PART} -n 8 --mem=64G -t 08:00:00 \
    --export=ALL,ACE_PROP_OUTPUT_ROOT="${OUTPUT_ROOT}" \
    -o "${LOG_DIR}/%j_sccomp_${ARM}.out" \
    -e "${LOG_DIR}/%j_sccomp_${ARM}.err" \
    --job-name="sccomp_${ARM}" \
    "${SCRIPT_DIR}/run_sccomp.sh" "${INTEGRATION}" "male" "${RESOLUTION}" --arm "${ARM}" ${SMOKE_FLAG})
  echo "  ${ARM}: job ${jid}"
done

echo ""
echo "Monitor: squeue -u \$USER | grep sccomp_"
