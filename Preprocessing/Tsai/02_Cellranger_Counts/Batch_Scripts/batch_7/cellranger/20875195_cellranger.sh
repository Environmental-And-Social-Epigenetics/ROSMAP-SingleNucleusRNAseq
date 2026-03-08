#!/bin/bash
#SBATCH --job-name=cr_20875195
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/02_Cellranger_Counts/Logs/Outs/cellranger_20875195_%j.out
#SBATCH --error=/om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/02_Cellranger_Counts/Logs/Errs/cellranger_20875195_%j.err
#SBATCH --mail-user=mabdel03@mit.edu
#SBATCH --mail-type=FAIL

# Cell Ranger count for patient 20875195 (Batch 7)
# Library IDs: D19-4108,D19-4108,D20-7436,D20-7436
# FASTQ directories: 4

set -e

# Add Cell Ranger to PATH
export PATH=/om2/user/mabdel03/apps/yard/cellranger-8.0.0:$PATH

# Create output directory
mkdir -p /om/scratch/Mon/mabdel03/Tsai_Data/Cellranger_Outputs/20875195

# Run Cell Ranger count
cellranger count \
    --create-bam=false \
    --include-introns=true \
    --nosecondary \
    --r1-length=26 \
    --id=20875195 \
    --transcriptome=/om2/user/mabdel03/yard/references/human/refdata-gex-GRCh38-2020-A \
    --sample=D19-4108,D19-4108,D20-7436,D20-7436 \
    --fastqs=/om/scratch/Mon/mabdel03/Tsai_Data/FASTQs/20875195/D19-4108/10x-3855H,/om/scratch/Mon/mabdel03/Tsai_Data/FASTQs/20875195/D19-4108/HKYYMDMXX,/om/scratch/Mon/mabdel03/Tsai_Data/FASTQs/20875195/D20-7436/10x-4819F,/om/scratch/Mon/mabdel03/Tsai_Data/FASTQs/20875195/D20-7436/10x-4826F \
    --output-dir=/om/scratch/Mon/mabdel03/Tsai_Data/Cellranger_Outputs/20875195

# Mark completion
echo "20875195" >> /om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/02_Cellranger_Counts/Tracking/cellranger_completed.txt

echo "Cell Ranger completed for 20875195"
