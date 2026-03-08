#!/bin/bash
#SBATCH --job-name=cr_20173942
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/02_Cellranger_Counts/Logs/Outs/cellranger_20173942_%j.out
#SBATCH --error=/om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/02_Cellranger_Counts/Logs/Errs/cellranger_20173942_%j.err
#SBATCH --mail-user=mabdel03@mit.edu
#SBATCH --mail-type=FAIL

# Cell Ranger count for patient 20173942 (Batch 5)
# Library IDs: D17-8796
# FASTQ directories: 1

set -e

# Add Cell Ranger to PATH
export PATH=/om2/user/mabdel03/apps/yard/cellranger-8.0.0:$PATH

# Create output directory
mkdir -p /om/scratch/Mon/mabdel03/Tsai_Data/Cellranger_Outputs/20173942

# Run Cell Ranger count
cellranger count \
    --create-bam=false \
    --include-introns=true \
    --nosecondary \
    --r1-length=26 \
    --id=20173942 \
    --transcriptome=/om2/user/mabdel03/yard/references/human/refdata-gex-GRCh38-2020-A \
    --sample=D17-8796 \
    --fastqs=/om/scratch/Mon/mabdel03/Tsai_Data/FASTQs/20173942/D17-8796 \
    --output-dir=/om/scratch/Mon/mabdel03/Tsai_Data/Cellranger_Outputs/20173942

# Mark completion
echo "20173942" >> /om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/02_Cellranger_Counts/Tracking/cellranger_completed.txt

echo "Cell Ranger completed for 20173942"
