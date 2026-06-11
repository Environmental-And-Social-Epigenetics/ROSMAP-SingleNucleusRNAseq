#!/bin/bash
# submit_all_rerun.sh — Resubmit all DEG jobs with tiered resources
#
# Tiered by h5ad file size:
#   Small  (<7GB):  8h, 64G   — most subtypes
#   Medium (7-22GB): 16h, 100G — In-VIP, Ex-L4_5, Ex-L4, Oli, broad_Inh
#   Large  (>22GB): 36h, 250G — Ex-L2_3, broad_Exc (pi_lhtsai only)
#
# Usage:
#   bash submit_all_rerun.sh [batch_prep_job_id]
#
#   If batch_prep_job_id is given, broad_Exc batch jobs depend on it
#   (because broad_Exc.h5ad is missing from celltype_splits_batch).

set -e

SCRIPT_DIR="/orcd/data/lhtsai/001/mabdel03/ROSMAP_Code/Transcriptomics/Analysis/ACE/DEG/Tsai"
BATCH_PREP_JOB="${1:-}"

PHENOTYPES=("tot_adverse_exp" "early_hh_ses" "ace_aggregate")

SMALL=("Endo" "In-PV_Chandelier" "Ex-NRGN" "Ex-L5_6" "Mic" "In-Rosehip" "Ex-L5" "In-SST" "OPC" "In-PV_Basket" "Ex-L5_6-CC" "Ast")
MEDIUM=("In-VIP" "Ex-L4_5" "Ex-L4" "Oli" "broad_Inh")
LARGE=("Ex-L2_3" "broad_Exc")

COUNT=0

submit_job() {
  local CT="$1" TIME="$2" MEM="$3" PARTS="$4" INTEGRATION="$5" PHENO="$6" DEP="$7"
  local DEP_FLAG=""
  [ -n "$DEP" ] && DEP_FLAG="--dependency=afterok:${DEP}"

  sbatch --parsable \
    -p "$PARTS" \
    -t "$TIME" --mem="$MEM" \
    $DEP_FLAG \
    -o "${SCRIPT_DIR}/logs/%j_${CT}_${PHENO}_${INTEGRATION}.out" \
    -e "${SCRIPT_DIR}/logs/%j_${CT}_${PHENO}_${INTEGRATION}.err" \
    --job-name="deg_${CT}_${PHENO}_${INTEGRATION}" \
    "${SCRIPT_DIR}/run_deg.sh" "$INTEGRATION" "$PHENO" "$CT" &

  COUNT=$((COUNT + 1))
}

for INTEGRATION in projid batch; do
  for PHENO in "${PHENOTYPES[@]}"; do

    for CT in "${SMALL[@]}"; do
      submit_job "$CT" "8:00:00" "64G" "pi_lhtsai,pi_manoli" "$INTEGRATION" "$PHENO" ""
    done

    for CT in "${MEDIUM[@]}"; do
      submit_job "$CT" "16:00:00" "100G" "pi_lhtsai,pi_manoli" "$INTEGRATION" "$PHENO" ""
    done

    for CT in "${LARGE[@]}"; do
      DEP=""
      if [ "$INTEGRATION" = "batch" ] && [ "$CT" = "broad_Exc" ] && [ -n "$BATCH_PREP_JOB" ]; then
        DEP="$BATCH_PREP_JOB"
      fi
      submit_job "$CT" "36:00:00" "250G" "pi_lhtsai" "$INTEGRATION" "$PHENO" "$DEP"
    done

  done
done

wait
echo "Submitted $COUNT jobs total (expected: 114)"
