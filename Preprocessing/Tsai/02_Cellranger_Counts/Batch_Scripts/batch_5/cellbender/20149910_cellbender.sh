#!/bin/bash
#SBATCH --job-name=cb_20149910
#SBATCH --gres=gpu:1
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output=/om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/02_Cellranger_Counts/Logs/Outs/cellbender_20149910_%j.out
#SBATCH --error=/om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/02_Cellranger_Counts/Logs/Errs/cellbender_20149910_%j.err
#SBATCH --mail-user=mabdel03@mit.edu
#SBATCH --mail-type=FAIL

# CellBender for patient 20149910 (Batch 5)

set -e

# Initialize conda
source /orcd/data/lhtsai/001/om2/mabdel03/miniforge3/etc/profile.d/conda.sh
conda activate /orcd/data/lhtsai/001/om2/mabdel03/conda_envs/Cellbender_env

# Create output directory
mkdir -p /om/scratch/Mon/mabdel03/Tsai/Cellbender_Output/20149910

# Check input exists
if [[ ! -f "/om/scratch/Mon/mabdel03/Tsai_Data/Cellranger_Outputs/20149910/outs/raw_feature_bc_matrix.h5" ]]; then
    echo "ERROR: Input file not found: /om/scratch/Mon/mabdel03/Tsai_Data/Cellranger_Outputs/20149910/outs/raw_feature_bc_matrix.h5"
    echo "20149910" >> /om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/02_Cellranger_Counts/Tracking/cellbender_failed.txt
    exit 1
fi

# Run CellBender
cellbender remove-background \
    --input /om/scratch/Mon/mabdel03/Tsai_Data/Cellranger_Outputs/20149910/outs/raw_feature_bc_matrix.h5 \
    --output /om/scratch/Mon/mabdel03/Tsai/Cellbender_Output/20149910/cellbender_output.h5 \
    --expected-cells 5000 \
    --total-droplets-included 20000 \
    --fpr 0.01 \
    --epochs 150 \
    --cuda

# Copy to permanent storage
mkdir -p /orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/Tsai/Preprocessed_Counts/20149910
cp /om/scratch/Mon/mabdel03/Tsai/Cellbender_Output/20149910/cellbender_output.h5 /orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/Tsai/Preprocessed_Counts/20149910/
cp /om/scratch/Mon/mabdel03/Tsai/Cellbender_Output/20149910/cellbender_output_filtered.h5 /orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/Tsai/Preprocessed_Counts/20149910/ 2>/dev/null || true

# Mark completion
echo "20149910" >> /om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/02_Cellranger_Counts/Tracking/cellbender_completed.txt

echo "CellBender completed for 20149910"
