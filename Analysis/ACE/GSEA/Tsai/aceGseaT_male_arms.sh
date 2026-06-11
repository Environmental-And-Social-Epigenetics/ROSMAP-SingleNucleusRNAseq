#!/bin/bash
# GSEA for the non-ANCOVA male ACE AD-model arms.
# One SLURM job per arm; each runs WebGestaltR GSEA over the signal-bearing cell
# types, reading that arm's DEG .rda files. Output:
#   ${ACE_OUTPUT_ROOT}/GSEA/Tsai/results_derived_batch_<ARM>/tot_adverse_exp/
#
# Usage:
#   bash aceGseaT_male_arms.sh
#   SMOKE_FLAG=--smoke bash aceGseaT_male_arms.sh
#   ARMS_OVERRIDE="MaleContAD" bash aceGseaT_male_arms.sh   # single arm

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"
source "${REPO_ROOT}/Analysis/ACE/_shared/arm_covariates.sh"

INTEGRATION="derived_batch"
PHENOTYPE="tot_adverse_exp"
OUTPUT_ROOT="${ACE_OUTPUT_ROOT}/GSEA/Tsai"
DEG_BASE="${ACE_OUTPUT_ROOT}/DEG/Tsai"
LOG_DIR="${OUTPUT_ROOT}/logs"
mkdir -p "${LOG_DIR}"

PART="-p pi_lhtsai,pi_manoli"
SMOKE_FLAG="${SMOKE_FLAG:-}"
CELLTYPES_CSV="$(IFS=,; echo "${SIGNAL_CELLTYPES[*]}")"

read -r -a ARMS <<< "${ARMS_OVERRIDE:-${NON_ANCOVA_ARMS[*]}}"

echo "=== ACE GSEA: male AD-model arms ==="
echo "Arms:       ${ARMS[*]}"
echo "Cell types: ${CELLTYPES_CSV}"
echo "Smoke:      ${SMOKE_FLAG:-(none)}"
echo ""

# Fan out ONE job per (arm x cell type): WebGestaltR over 8 DBs for a single cell
# type is ~1-2h, so 6 cell types in one job overran the 12h wall. Per-CT jobs run
# concurrently and finish well within the limit. ACE_GSEA_SKIP_EXISTING=1 makes
# reruns skip already-written per-DB .rds (so a resubmit resumes cheaply).
for ARM in "${ARMS[@]}"; do
  DEG_ROOT="${DEG_BASE}/results_${INTEGRATION}_${ARM}_AllCellTypes"
  if [[ ! -d "${DEG_ROOT}/${PHENOTYPE}" ]]; then
    echo "  WARN: DEG dir missing for ${ARM}: ${DEG_ROOT}/${PHENOTYPE} -- skipping"
    continue
  fi
  for CT in "${SIGNAL_CELLTYPES[@]}"; do
    jid=$(sbatch --parsable ${PART} -n 8 --mem=64G -t 06:00:00 \
      --export=ALL,LAUNCHER_SCRIPT_DIR="${SCRIPT_DIR}",GSEA_CELLTYPES="${CT}",SMOKE_FLAG="${SMOKE_FLAG}",ACE_GSEA_SKIP_EXISTING=1 \
      -o "${LOG_DIR}/%j_gsea_${ARM}_${CT}.out" \
      -e "${LOG_DIR}/%j_gsea_${ARM}_${CT}.err" \
      --job-name="gsea_${ARM}_${CT}" \
      "${SCRIPT_DIR}/run_enrichment.sh" "${INTEGRATION}" "${PHENOTYPE}" "Male" "${DEG_ROOT}" "${OUTPUT_ROOT}" "${ARM}")
    echo "  ${ARM} / ${CT}: job ${jid}"
  done
done

echo ""
echo "Monitor: squeue -u \$USER | grep gsea_"
