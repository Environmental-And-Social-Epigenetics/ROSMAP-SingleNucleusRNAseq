#!/bin/bash
#SBATCH --job-name=cr_10490993
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/02_Cellranger_Counts/Logs/Outs/cellranger_10490993_%j.out
#SBATCH --error=/om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/02_Cellranger_Counts/Logs/Errs/cellranger_10490993_%j.err
#SBATCH --mail-user=mabdel03@mit.edu
#SBATCH --mail-type=FAIL

# Cell Ranger count for patient 10490993 (Batch 2)
# Library IDs: D19-5951,D19-5951,D20-7452,D20-7452
# FASTQ directories: 4

set -e

# Add Cell Ranger to PATH
export PATH=/om2/user/mabdel03/apps/yard/cellranger-8.0.0:$PATH

# Create output directory
mkdir -p /om/scratch/Mon/mabdel03/Tsai_Data/Cellranger_Outputs/10490993

# Run Cell Ranger count
cellranger count \
    --create-bam=false \
    --include-introns=true \
    --nosecondary \
    --r1-length=26 \
    --id=10490993 \
    --transcriptome=/om2/user/mabdel03/yard/references/human/refdata-gex-GRCh38-2020-A \
    --sample=D19-5951,D19-5951,D20-7452,D20-7452 \
    --fastqs=/om/scratch/Mon/mabdel03/Tsai_Data/FASTQs/10490993/D19-5951/HKMG5DMXX,/om/scratch/Mon/mabdel03/Tsai_Data/FASTQs/10490993/D19-5951/HL2FCDMXX,/om/scratch/Mon/mabdel03/Tsai_Data/FASTQs/10490993/D20-7452/10x-4819F,/om/scratch/Mon/mabdel03/Tsai_Data/FASTQs/10490993/D20-7452/10x-4826F \
    --output-dir=/om/scratch/Mon/mabdel03/Tsai_Data/Cellranger_Outputs/10490993

# Mark completion
echo "10490993" >> /om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/02_Cellranger_Counts/Tracking/cellranger_completed.txt

echo "Cell Ranger completed for 10490993"
