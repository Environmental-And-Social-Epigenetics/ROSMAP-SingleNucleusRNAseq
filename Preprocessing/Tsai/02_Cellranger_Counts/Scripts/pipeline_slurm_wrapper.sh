#!/bin/bash
#SBATCH --job-name=pipeline_master
#SBATCH --partition=mit_preemptable
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/02_Cellranger_Counts/Logs/Outs/pipeline_master_%j.out
#SBATCH --error=/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/02_Cellranger_Counts/Logs/Errs/pipeline_master_%j.err
#SBATCH --mail-user=mabdel03@mit.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --requeue

# Master orchestrator that runs batches sequentially
# Requeues automatically if preempted, will continue from where it left off
# Uses tracking files to know what's already done

# Use absolute paths (SLURM copies scripts, breaking BASH_SOURCE)
PIPELINE_DIR="/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/02_Cellranger_Counts"
SCRIPT_DIR="${PIPELINE_DIR}/Scripts"
source "${PIPELINE_DIR}/Config/cellranger_config.sh"

# Find next incomplete batch
find_next_batch() {
    for batch_num in $(seq 1 16); do
        # Check if all patients in this batch have completed CellBender
        local batch_patients=$(tail -n +2 "${TRACKING_DIR}/batch_assignments.csv" | awk -F',' -v b="$batch_num" '$3 == b {print $1}')
        local all_done=true
        
        for projid in $batch_patients; do
            if ! grep -q "^${projid}$" "${TRACKING_DIR}/cellbender_completed.txt" 2>/dev/null; then
                all_done=false
                break
            fi
        done
        
        if [[ "$all_done" == "false" ]]; then
            echo "$batch_num"
            return
        fi
    done
    echo "0"  # All done
}

# Check if Cell Ranger is done for a batch
is_cellranger_done() {
    local batch_num=$1
    local batch_patients=$(tail -n +2 "${TRACKING_DIR}/batch_assignments.csv" | awk -F',' -v b="$batch_num" '$3 == b {print $1}')
    
    for projid in $batch_patients; do
        if ! grep -q "^${projid}$" "${TRACKING_DIR}/cellranger_completed.txt" 2>/dev/null; then
            return 1
        fi
    done
    return 0
}

echo "=============================================="
echo "Pipeline Master - Starting/Resuming"
echo "Time: $(date)"
echo "=============================================="

NEXT_BATCH=$(find_next_batch)

if [[ "$NEXT_BATCH" == "0" ]]; then
    echo "All batches complete!"
    exit 0
fi

echo "Next batch to process: $NEXT_BATCH"

# Run batches until time runs out or all done
for BATCH_NUM in $(seq "$NEXT_BATCH" 16); do
    echo ""
    echo "=== Processing Batch $BATCH_NUM ==="
    
    # Determine what to skip
    SKIP_ARG=""
    if is_cellranger_done "$BATCH_NUM"; then
        SKIP_ARG="--skip-cellranger-batch $BATCH_NUM"
        echo "Cell Ranger already complete for batch $BATCH_NUM"
    fi
    
    # Run single batch
    "${SCRIPT_DIR}/run_all_batches.sh" --start-batch "$BATCH_NUM" --end-batch "$BATCH_NUM" $SKIP_ARG
done

echo ""
echo "Pipeline complete!"
