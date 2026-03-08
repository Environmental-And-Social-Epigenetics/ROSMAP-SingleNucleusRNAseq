#!/bin/bash
#SBATCH --job-name=cr_50105301
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/02_Cellranger_Counts/Logs/Outs/cellranger_50105301_%j.out
#SBATCH --error=/om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/02_Cellranger_Counts/Logs/Errs/cellranger_50105301_%j.err
#SBATCH --mail-user=mabdel03@mit.edu
#SBATCH --mail-type=FAIL

# Cell Ranger count for patient 50105301 (Batch 11)
# Library IDs: D19-4131,D19-4131,D20-7458,D20-7458
# FASTQ directories: 4

set -e

# Add Cell Ranger to PATH
export PATH=/om2/user/mabdel03/apps/yard/cellranger-8.0.0:$PATH

# Create output directory
mkdir -p /om/scratch/Mon/mabdel03/Tsai_Data/Cellranger_Outputs/50105301

# Run Cell Ranger count
cellranger count \
    --create-bam=false \
    --include-introns=true \
    --nosecondary \
    --r1-length=26 \
    --id=50105301 \
    --transcriptome=/om2/user/mabdel03/yard/references/human/refdata-gex-GRCh38-2020-A \
    --sample=D19-4131,D19-4131,D20-7458,D20-7458 \
    --fastqs=/om/scratch/Mon/mabdel03/Tsai_Data/FASTQs/50105301/D19-4131/HKYYMDMXX,/om/scratch/Mon/mabdel03/Tsai_Data/FASTQs/50105301/D19-4131/10x-3855H,/om/scratch/Mon/mabdel03/Tsai_Data/FASTQs/50105301/D20-7458/10x-4819F,/om/scratch/Mon/mabdel03/Tsai_Data/FASTQs/50105301/D20-7458/10x-4826F \
    --output-dir=/om/scratch/Mon/mabdel03/Tsai_Data/Cellranger_Outputs/50105301

# Mark completion
echo "50105301" >> /om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/02_Cellranger_Counts/Tracking/cellranger_completed.txt

echo "Cell Ranger completed for 50105301"
