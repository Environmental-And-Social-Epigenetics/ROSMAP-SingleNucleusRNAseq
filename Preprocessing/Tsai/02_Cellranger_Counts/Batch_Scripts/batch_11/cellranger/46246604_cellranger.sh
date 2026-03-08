#!/bin/bash
#SBATCH --job-name=cr_46246604
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/02_Cellranger_Counts/Logs/Outs/cellranger_46246604_%j.out
#SBATCH --error=/om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/02_Cellranger_Counts/Logs/Errs/cellranger_46246604_%j.err
#SBATCH --mail-user=mabdel03@mit.edu
#SBATCH --mail-type=FAIL

# Cell Ranger count for patient 46246604 (Batch 11)
# Library IDs: D19-4783,D19-4783,D20-7450,D20-7450
# FASTQ directories: 4

set -e

# Add Cell Ranger to PATH
export PATH=/om2/user/mabdel03/apps/yard/cellranger-8.0.0:$PATH

# Create output directory
mkdir -p /om/scratch/Mon/mabdel03/Tsai_Data/Cellranger_Outputs/46246604

# Run Cell Ranger count
cellranger count \
    --create-bam=false \
    --include-introns=true \
    --nosecondary \
    --r1-length=26 \
    --id=46246604 \
    --transcriptome=/om2/user/mabdel03/yard/references/human/refdata-gex-GRCh38-2020-A \
    --sample=D19-4783,D19-4783,D20-7450,D20-7450 \
    --fastqs=/om/scratch/Mon/mabdel03/Tsai_Data/FASTQs/46246604/D19-4783/10x-3877H,/om/scratch/Mon/mabdel03/Tsai_Data/FASTQs/46246604/D19-4783/HL2TKDMXX,/om/scratch/Mon/mabdel03/Tsai_Data/FASTQs/46246604/D20-7450/10x-4819F,/om/scratch/Mon/mabdel03/Tsai_Data/FASTQs/46246604/D20-7450/10x-4826F \
    --output-dir=/om/scratch/Mon/mabdel03/Tsai_Data/Cellranger_Outputs/46246604

# Mark completion
echo "46246604" >> /om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/02_Cellranger_Counts/Tracking/cellranger_completed.txt

echo "Cell Ranger completed for 46246604"
