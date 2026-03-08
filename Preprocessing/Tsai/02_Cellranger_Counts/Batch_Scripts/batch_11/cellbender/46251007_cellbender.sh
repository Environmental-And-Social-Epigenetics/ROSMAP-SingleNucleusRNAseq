#!/bin/bash
#SBATCH --job-name=cb_46251007
#SBATCH --gres=gpu:1
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=/om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/02_Cellranger_Counts/Logs/Outs/cellbender_46251007_%j.out
#SBATCH --error=/om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/02_Cellranger_Counts/Logs/Errs/cellbender_46251007_%j.err
#SBATCH --mail-user=mabdel03@mit.edu
#SBATCH --mail-type=FAIL

# CellBender for patient 46251007 (Batch 11)

set -e

# Initialize conda
source /orcd/data/lhtsai/001/om2/mabdel03/miniforge3/etc/profile.d/conda.sh
conda activate /orcd/data/lhtsai/001/om2/mabdel03/conda_envs/Cellbender_env

# Create output directory
mkdir -p /om/scratch/Mon/mabdel03/Tsai/Cellbender_Output/46251007

# Check input exists
if [[ ! -f "/om/scratch/Mon/mabdel03/Tsai_Data/Cellranger_Outputs/46251007/outs/raw_feature_bc_matrix.h5" ]]; then
    echo "ERROR: Input file not found: /om/scratch/Mon/mabdel03/Tsai_Data/Cellranger_Outputs/46251007/outs/raw_feature_bc_matrix.h5"
    echo "46251007" >> /om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/02_Cellranger_Counts/Tracking/cellbender_failed.txt
    exit 1
fi

# Run CellBender
cellbender remove-background \
    --input /om/scratch/Mon/mabdel03/Tsai_Data/Cellranger_Outputs/46251007/outs/raw_feature_bc_matrix.h5 \
    --output /om/scratch/Mon/mabdel03/Tsai/Cellbender_Output/46251007/cellbender_output.h5 \
    --expected-cells 5000 \
    --total-droplets-included 20000 \
    --fpr 0.01 \
    --epochs 150 \
    --cuda

# Copy to permanent storage
mkdir -p /orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/Tsai/Preprocessed_Counts/46251007
cp /om/scratch/Mon/mabdel03/Tsai/Cellbender_Output/46251007/cellbender_output.h5 /orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/Tsai/Preprocessed_Counts/46251007/
cp /om/scratch/Mon/mabdel03/Tsai/Cellbender_Output/46251007/cellbender_output_filtered.h5 /orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/Tsai/Preprocessed_Counts/46251007/ 2>/dev/null || true

# Mark completion
echo "46251007" >> /om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/02_Cellranger_Counts/Tracking/cellbender_completed.txt

echo "CellBender completed for 46251007"
