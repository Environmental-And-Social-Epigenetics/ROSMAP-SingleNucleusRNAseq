#!/bin/bash
#SBATCH --job-name=cr_66406040
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/02_Cellranger_Counts/Logs/Outs/cellranger_66406040_%j.out
#SBATCH --error=/om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/02_Cellranger_Counts/Logs/Errs/cellranger_66406040_%j.err
#SBATCH --mail-user=mabdel03@mit.edu
#SBATCH --mail-type=FAIL

# Cell Ranger count for patient 66406040 (Batch 14)
# Library IDs: D19-7295,D19-7295
# FASTQ directories: 2

set -e

# Add Cell Ranger to PATH
export PATH=/om2/user/mabdel03/apps/yard/cellranger-8.0.0:$PATH

# Create output directory
mkdir -p /om/scratch/Mon/mabdel03/Tsai_Data/Cellranger_Outputs/66406040

# Run Cell Ranger count
cellranger count \
    --create-bam=false \
    --include-introns=true \
    --nosecondary \
    --r1-length=26 \
    --id=66406040 \
    --transcriptome=/om2/user/mabdel03/yard/references/human/refdata-gex-GRCh38-2020-A \
    --sample=D19-7295,D19-7295 \
    --fastqs=/om/scratch/Mon/mabdel03/Tsai_Data/FASTQs/66406040/D19-7295/HL2CHDMXX,/om/scratch/Mon/mabdel03/Tsai_Data/FASTQs/66406040/D19-7295/HLFG3DMXX \
    --output-dir=/om/scratch/Mon/mabdel03/Tsai_Data/Cellranger_Outputs/66406040

# Mark completion
echo "66406040" >> /om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/02_Cellranger_Counts/Tracking/cellranger_completed.txt

echo "Cell Ranger completed for 66406040"
