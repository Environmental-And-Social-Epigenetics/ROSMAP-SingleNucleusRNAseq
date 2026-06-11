#!/bin/bash
# Efficient orchestrator for the male ACE AD-confounding DEG arms.
#
# Fans out ONE SLURM job per cell type (all concurrent), each running the new
# arms via --celltype so the (large) h5ad is read once per cell type. Then
# submits the cross-model summarizer with an afterok dependency on the DEG jobs.
#
# Arms launched (NEW):
#   MaleNiaReagan : ~ age_death + pmi + niareagansc + tot_adverse_exp        (Model 0)
#   MaleAncovaAD  : per AD var v: ~ age_death + pmi + v + tot_adverse_exp     (Model 3)
#                   v in {niareagansc, tangsqrt, braaksc, amylsqrt}
#
# The existing arms (NoADadj, ContAD, BinaryAD, AceByAD) are already complete on
# disk and are NOT re-run unless RERUN_EXISTING=1.
#
# Env knobs:
#   SMOKE_FLAG=--smoke   dry run (load + pseudobulk dims, no DESeq2)
#   RERUN_EXISTING=1     also re-fire the four existing arm wrappers per cell type
#
# Usage:
#   bash run_male_ad_models.sh
#   SMOKE_FLAG=--smoke bash run_male_ad_models.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

OUTPUT_ROOT="${ANALYSIS_OUTPUT_ROOT}/ACE/DEG/Tsai"
LOG_DIR="${OUTPUT_ROOT}/logs"
mkdir -p "${LOG_DIR}"

PART="-p pi_lhtsai,pi_manoli"   # NEVER pi_tpoggio
SMOKE_FLAG="${SMOKE_FLAG:-}"
RERUN_EXISTING="${RERUN_EXISTING:-0}"

# New arms (run per cell type). Existing arms appended only if RERUN_EXISTING=1.
NEW_ARMS=("MaleNiaReagan" "MaleAncovaAD")
EXISTING_ARMS=("MaleNoADadj" "MaleContAD" "MaleBinaryAD" "MaleAceByAD")
ARMS=("${NEW_ARMS[@]}")
if [[ "${RERUN_EXISTING}" == "1" ]]; then
  ARMS=("${NEW_ARMS[@]}" "${EXISTING_ARMS[@]}")
fi

# Cell types by size tier. NOTE: pooled inhibitory uses file broad_Inh.h5ad
# (label "Inh"); broad_Exc is intentionally excluded (>2^31 nnz / 101 GB).
L_CT=(Ex-L2_3 broad_Inh Oli)
M_CT=(Ex-L4 Ex-L4_5 In-VIP Ast Ex-L5_6-CC In-PV_Basket OPC)
S_CT=(Ex-L5 In-SST Mic In-Rosehip Ex-L5_6 In-PV_Chandelier Ex-NRGN Endo)

tier_mem()  { case "$1" in L) echo 200G;; M) echo 96G;; S) echo 48G;; esac; }
tier_time() { case "$1" in L) echo 08:00:00;; M) echo 04:00:00;; S) echo 02:00:00;; esac; }
tier_ntask(){ case "$1" in L) echo 4;; M) echo 4;; S) echo 2;; esac; }

echo "=== Male ACE AD-model orchestrator ==="
echo "Output root:   ${OUTPUT_ROOT}"
echo "Arms:          ${ARMS[*]}"
echo "Smoke flag:    ${SMOKE_FLAG:-(none)}"
echo "Rerun existing:${RERUN_EXISTING}"
echo ""

DEG_JOBS=()

submit_celltype() {
  local ct="$1" tier="$2"
  local mem time ntask
  mem="$(tier_mem "$tier")"; time="$(tier_time "$tier")"; ntask="$(tier_ntask "$tier")"

  # Build the per-cell-type command: run each arm's wrapper with CELLTYPE set.
  local inner=""
  for arm in "${ARMS[@]}"; do
    inner+="CELLTYPE=${ct} SMOKE_FLAG=${SMOKE_FLAG} bash ${SCRIPT_DIR}/aceDegT_${arm}_AllCellTypes.sh; "
  done

  local jid
  jid=$(sbatch --parsable ${PART} -n "${ntask}" --mem="${mem}" -t "${time}" \
    -o "${LOG_DIR}/%j_maleAD_${ct}.out" \
    -e "${LOG_DIR}/%j_maleAD_${ct}.err" \
    --job-name="maleAD_${ct}" \
    --export=ALL,SLURM_SUBMIT_DIR="${SCRIPT_DIR}" \
    --wrap="set -e; ${inner}")
  DEG_JOBS+=("${jid}")
  printf "  %-18s [%s mem=%s t=%s n=%s] job %s\n" "${ct}" "${tier}" "${mem}" "${time}" "${ntask}" "${jid}"
}

for ct in "${L_CT[@]}"; do submit_celltype "${ct}" L; done
for ct in "${M_CT[@]}"; do submit_celltype "${ct}" M; done
for ct in "${S_CT[@]}"; do submit_celltype "${ct}" S; done

echo ""
echo "Submitted ${#DEG_JOBS[@]} per-cell-type DEG jobs."

# Cross-model summary, gated on all DEG jobs succeeding.
DEP="$(IFS=:; echo "${DEG_JOBS[*]}")"
SUM_JOB=$(sbatch --parsable ${PART} -n 2 --mem=32G -t 01:00:00 \
  --dependency=afterok:"${DEP}" \
  -o "${LOG_DIR}/%j_summarize_male_ad.out" \
  -e "${LOG_DIR}/%j_summarize_male_ad.err" \
  --job-name="maleAD_summary" \
  --wrap="bash ${SCRIPT_DIR}/report/summarize_ACE_AD_models.sh")

echo "Summary job ${SUM_JOB} (dependency afterok:${DEP})"
echo ""
echo "Monitor with: squeue -u \$USER"
