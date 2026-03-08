#!/bin/bash
#SBATCH --job-name=cr_50107583
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/02_Cellranger_Counts/Logs/Outs/cellranger_50107583_%j.out
#SBATCH --error=/om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/02_Cellranger_Counts/Logs/Errs/cellranger_50107583_%j.err
#SBATCH --mail-user=mabdel03@mit.edu
#SBATCH --mail-type=FAIL

# Cell Ranger count for patient 50107583 (Batch 12)
# Library IDs: D19-5044,D19-5044
# FASTQ directories: 2

set -e

# Add Cell Ranger to PATH
export PATH=/om2/user/mabdel03/apps/yard/cellranger-8.0.0:$PATH

# Create output directory
mkdir -p /om/scratch/Mon/mabdel03/Tsai_Data/Cellranger_Outputs/50107583

# Run Cell Ranger count
cellranger count \
    --create-bam=false \
    --include-introns=true \
    --nosecondary \
    --r1-length=26 \
    --id=50107583 \
    --transcriptome=/om2/user/mabdel03/yard/references/human/refdata-gex-GRCh38-2020-A \
    --sample=D19-5044,D19-5044 \
    --fastqs=/om/scratch/Mon/mabdel03/Tsai_Data/FASTQs/50107583/D19-5044/HL3JJDMXX,/om/scratch/Mon/mabdel03/Tsai_Data/FASTQs/50107583/D19-5044/10x-3914H \
    --output-dir=/om/scratch/Mon/mabdel03/Tsai_Data/Cellranger_Outputs/50107583

# Mark completion
echo "50107583" >> /om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/02_Cellranger_Counts/Tracking/cellranger_completed.txt

echo "Cell Ranger completed for 50107583"
