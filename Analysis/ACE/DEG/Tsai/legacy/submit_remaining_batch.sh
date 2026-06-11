#!/bin/bash
# Submit remaining batch DEG jobs (early_hh_ses and ace_aggregate)
# Run this after some jobs complete to free up QOS slots
# Also submits broad_Exc jobs if broad_Exc.h5ad now exists

set -e

SCRIPT_DIR="/orcd/data/lhtsai/001/mabdel03/ROSMAP_Code/Transcriptomics/Analysis/ACE/DEG/Tsai"
PARTS="pi_lhtsai,pi_tpoggio,pi_manoli"

CELLTYPES=("Ast" "Endo" "Ex-L2_3" "Ex-L4" "Ex-L4_5" "Ex-L5" "Ex-L5_6" "Ex-L5_6-CC" "Ex-NRGN" "In-PV_Basket" "In-PV_Chandelier" "In-Rosehip" "In-SST" "In-VIP" "Mic" "Oli" "OPC" "broad_Inh")

# Add broad_Exc if it exists
if [ -f "${SCRIPT_DIR}/celltype_splits_batch/broad_Exc.h5ad" ]; then
  CELLTYPES+=("broad_Exc")
fi

COUNT=0
for PHENO in early_hh_ses ace_aggregate; do
  for CT in "${CELLTYPES[@]}"; do
    # Check if result already exists (skip if so)
    sbatch --parsable \
      -p $PARTS \
      -o "${SCRIPT_DIR}/logs/%j_${CT}_${PHENO}_batch.out" \
      -e "${SCRIPT_DIR}/logs/%j_${CT}_${PHENO}_batch.err" \
      --job-name="deg_${CT}_${PHENO}_batch" \
      "${SCRIPT_DIR}/run_deg.sh" "batch" "$PHENO" "$CT" &
    COUNT=$((COUNT + 1))
  done
done

wait
echo "Submitted $COUNT batch jobs"
