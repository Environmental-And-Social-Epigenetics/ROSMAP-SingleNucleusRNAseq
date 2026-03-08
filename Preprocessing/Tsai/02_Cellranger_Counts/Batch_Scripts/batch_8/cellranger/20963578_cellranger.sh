#!/bin/bash
#SBATCH --job-name=cr_20963578
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=/om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/02_Cellranger_Counts/Logs/Outs/cellranger_20963578_%j.out
#SBATCH --error=/om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/02_Cellranger_Counts/Logs/Errs/cellranger_20963578_%j.err
#SBATCH --mail-user=mabdel03@mit.edu
#SBATCH --mail-type=FAIL

# Cell Ranger count for patient 20963578 (Batch 8)
# Library IDs: D19-4136,D19-4136,D19-5029,D19-5029
# FASTQ directories: 4

set -e

# Add Cell Ranger to PATH
export PATH=/om2/user/mabdel03/apps/yard/cellranger-8.0.0:$PATH

# Create output directory
mkdir -p /om/scratch/Mon/mabdel03/Tsai_Data/Cellranger_Outputs/20963578

# Run Cell Ranger count
cellranger count \
    --create-bam=false \
    --include-introns=true \
    --nosecondary \
    --r1-length=26 \
    --id=20963578 \
    --transcriptome=/om2/user/mabdel03/yard/references/human/refdata-gex-GRCh38-2020-A \
    --sample=D19-4136,D19-4136,D19-5029,D19-5029 \
    --fastqs=/om/scratch/Mon/mabdel03/Tsai_Data/FASTQs/20963578/D19-4136/HL32WDMXX,/om/scratch/Mon/mabdel03/Tsai_Data/FASTQs/20963578/D19-4136/10x-3862H,/om/scratch/Mon/mabdel03/Tsai_Data/FASTQs/20963578/D19-5029/10x-3914H,/om/scratch/Mon/mabdel03/Tsai_Data/FASTQs/20963578/D19-5029/HL3JJDMXX \
    --output-dir=/om/scratch/Mon/mabdel03/Tsai_Data/Cellranger_Outputs/20963578

# Mark completion
echo "20963578" >> /om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/02_Cellranger_Counts/Tracking/cellranger_completed.txt

echo "Cell Ranger completed for 20963578"
