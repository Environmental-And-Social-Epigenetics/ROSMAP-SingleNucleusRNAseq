#!/bin/bash
# submit_remaining_rerun.sh — Submit DEG jobs that hit the QOS limit
#
# Usage:
#   bash submit_remaining_rerun.sh [batch_prep_job_id]

set -e

SCRIPT_DIR="/orcd/data/lhtsai/001/mabdel03/ROSMAP_Code/Transcriptomics/Analysis/ACE/DEG/Tsai"
BATCH_PREP_JOB="${1:-}"

COUNT=0

submit_job() {
  local CT="$1" TIME="$2" MEM="$3" PARTS="$4" INTEGRATION="$5" PHENO="$6" DEP="$7"
  local DEP_FLAG=""
  [ -n "$DEP" ] && DEP_FLAG="--dependency=afterok:${DEP}"

  local JOBNAME="deg_${CT}_${PHENO}_${INTEGRATION}"

  # Skip if already queued
  if squeue -u mabdel03 -o "%j" 2>/dev/null | grep -qx "$JOBNAME"; then
    echo "  SKIP (already queued): $JOBNAME"
    return
  fi

  sbatch --parsable \
    -p "$PARTS" \
    -t "$TIME" --mem="$MEM" \
    $DEP_FLAG \
    -o "${SCRIPT_DIR}/logs/%j_${CT}_${PHENO}_${INTEGRATION}.out" \
    -e "${SCRIPT_DIR}/logs/%j_${CT}_${PHENO}_${INTEGRATION}.err" \
    --job-name="$JOBNAME" \
    "${SCRIPT_DIR}/run_deg.sh" "$INTEGRATION" "$PHENO" "$CT"

  COUNT=$((COUNT + 1))
}

# Tier definitions (same as submit_all_rerun.sh)
SMALL_TIME="8:00:00";  SMALL_MEM="64G";  SMALL_PARTS="pi_lhtsai,pi_manoli"
MED_TIME="16:00:00";   MED_MEM="100G";   MED_PARTS="pi_lhtsai,pi_manoli"
LARGE_TIME="36:00:00"; LARGE_MEM="250G"; LARGE_PARTS="pi_lhtsai"

# --- Missing jobs (from QOS limit during initial submission) ---

# Small tier
for CT in Ast Endo Ex-L5_6 Ex-L5_6-CC Ex-NRGN Ex-L5 In-PV_Basket In-PV_Chandelier In-Rosehip In-SST Mic OPC; do
  for combo in \
    "ace_aggregate:projid" "ace_aggregate:batch" \
    "early_hh_ses:projid" "early_hh_ses:batch" \
    "tot_adverse_exp:projid" "tot_adverse_exp:batch"; do
    PHENO="${combo%%:*}"
    INT="${combo##*:}"
    submit_job "$CT" "$SMALL_TIME" "$SMALL_MEM" "$SMALL_PARTS" "$INT" "$PHENO" ""
  done
done

# Medium tier
for CT in In-VIP Ex-L4_5 Ex-L4 Oli broad_Inh; do
  for combo in \
    "ace_aggregate:projid" "ace_aggregate:batch" \
    "early_hh_ses:projid" "early_hh_ses:batch" \
    "tot_adverse_exp:projid" "tot_adverse_exp:batch"; do
    PHENO="${combo%%:*}"
    INT="${combo##*:}"
    submit_job "$CT" "$MED_TIME" "$MED_MEM" "$MED_PARTS" "$INT" "$PHENO" ""
  done
done

# Large tier
for CT in Ex-L2_3 broad_Exc; do
  for combo in \
    "ace_aggregate:projid" "ace_aggregate:batch" \
    "early_hh_ses:projid" "early_hh_ses:batch" \
    "tot_adverse_exp:projid" "tot_adverse_exp:batch"; do
    PHENO="${combo%%:*}"
    INT="${combo##*:}"
    DEP=""
    if [ "$INT" = "batch" ] && [ "$CT" = "broad_Exc" ] && [ -n "$BATCH_PREP_JOB" ]; then
      DEP="$BATCH_PREP_JOB"
    fi
    submit_job "$CT" "$LARGE_TIME" "$LARGE_MEM" "$LARGE_PARTS" "$INT" "$PHENO" "$DEP"
  done
done

echo ""
echo "Submitted $COUNT new jobs"
