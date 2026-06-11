#!/bin/bash
# Submit all DEG jobs: per cell type × per phenotype × per integration
# Each job runs a single (integration, phenotype, celltype) combination
#
# Usage: bash submit_parallel.sh [batch_prep_job_id]
#   If batch_prep_job_id is given, batch DEG jobs depend on it.
#   Projid jobs are submitted without dependency (prep already done).

set -e

SCRIPT_DIR="/orcd/data/lhtsai/001/mabdel03/ROSMAP_Code/Transcriptomics/Analysis/ACE/DEG/Tsai"
PARTS="pi_lhtsai,pi_manoli"
BATCH_PREP_JOB="${1:-}"

PHENOTYPES=("tot_adverse_exp" "early_hh_ses" "ace_aggregate")
CELLTYPES=(
  "Ast" "Endo" "Ex-L2_3" "Ex-L4" "Ex-L4_5" "Ex-L5" "Ex-L5_6" "Ex-L5_6-CC" "Ex-NRGN"
  "In-PV_Basket" "In-PV_Chandelier" "In-Rosehip" "In-SST" "In-VIP"
  "Mic" "Oli" "OPC"
  "broad_Exc" "broad_Inh"
)

COUNT=0

for INTEGRATION in projid batch; do
  DEP_FLAG=""
  if [ "$INTEGRATION" = "batch" ] && [ -n "$BATCH_PREP_JOB" ]; then
    DEP_FLAG="--dependency=afterok:${BATCH_PREP_JOB}"
  fi

  for PHENO in "${PHENOTYPES[@]}"; do
    for CT in "${CELLTYPES[@]}"; do
      # Skip broad types for batch if their h5ad doesn't exist yet
      # (they'll be created by the prep job, dependency handles ordering)

      JOBNAME="deg_${CT}_${PHENO}_${INTEGRATION}"
      sbatch --parsable \
        -p $PARTS \
        $DEP_FLAG \
        -o "${SCRIPT_DIR}/logs/%j_${CT}_${PHENO}_${INTEGRATION}.out" \
        -e "${SCRIPT_DIR}/logs/%j_${CT}_${PHENO}_${INTEGRATION}.err" \
        --job-name="$JOBNAME" \
        "${SCRIPT_DIR}/run_deg.sh" "$INTEGRATION" "$PHENO" "$CT" &

      COUNT=$((COUNT + 1))
    done
  done
done

wait
echo "Submitted $COUNT jobs total"
