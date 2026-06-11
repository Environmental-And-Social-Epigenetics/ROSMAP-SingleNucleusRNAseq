#!/bin/bash
# Submit broad_Exc DEG jobs (Python pseudobulk → R DESeq2)
# Usage: bash submit_broad_Exc.sh

set -e

SCRIPT_DIR="/orcd/data/lhtsai/001/mabdel03/ROSMAP_Code/Transcriptomics/Analysis/ACE/DEG/Tsai"

for INTEGRATION in projid batch; do
  # Step 1: Python pseudobulk
  PB_JOB=$(sbatch --parsable \
    -t 12:00:00 --mem=64G \
    -o "${SCRIPT_DIR}/logs/%j_broad_Exc_pseudobulk_${INTEGRATION}.out" \
    -e "${SCRIPT_DIR}/logs/%j_broad_Exc_pseudobulk_${INTEGRATION}.err" \
    --job-name="pb_Exc_${INTEGRATION}" \
    "${SCRIPT_DIR}/run_pseudobulk_broad_Exc.sh" "${INTEGRATION}")
  echo "Submitted pseudobulk ${INTEGRATION}: ${PB_JOB}"

  # Step 2: R DESeq2 (depends on pseudobulk completing)
  DEG_JOB=$(sbatch --parsable \
    -t 4:00:00 --mem=32G \
    --dependency=afterok:${PB_JOB} \
    -o "${SCRIPT_DIR}/logs/%j_broad_Exc_deseq_${INTEGRATION}.out" \
    -e "${SCRIPT_DIR}/logs/%j_broad_Exc_deseq_${INTEGRATION}.err" \
    --job-name="deseq_Exc_${INTEGRATION}" \
    "${SCRIPT_DIR}/run_deseq_broad_Exc.sh" "${INTEGRATION}")
  echo "Submitted DESeq2 ${INTEGRATION}: ${DEG_JOB} (depends on ${PB_JOB})"
done

echo "Done. 4 jobs submitted (2 pseudobulk + 2 DESeq2)."
