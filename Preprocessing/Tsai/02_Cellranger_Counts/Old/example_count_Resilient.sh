#!/bin/bash
#SBATCH -t 47:00:00
#SBATCH -n 32
#SBATCH --mem=128G
#SBATCH --output=/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/Tsai/Preprocessing/Counts/Resilient/Batch_Scripts_Outs/slurm-%j.out
#SBATCH --error=/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/Tsai/Preprocessing/Counts/Resilient/Batch_Scripts_Err/slurm-%j.err
#SBATCH --mail-user=mabdel03@mit.edu
#SBATCH --mail-type=FAIL
export PATH=/om2/user/$USER/apps/yard/cellranger-8.0.0:$PATH
cellranger count --create-bam false --include-introns true --nosecondary --r1-length 26 --id D19-7499 --transcriptome=/orcd/data/lhtsai/001/om2/mabdel03/yard/references/human/refdata-gex-GRCh38-2020-A --sample D19-7499 --fastqs /om/scratch/Mon/mabdel03/Tsai/Resilient/FASTQs/10101589/D19-7499_1,/om/scratch/Mon/mabdel03/Tsai/Resilient/FASTQs/10101589/D19-7499_2 --output-dir=/om/scratch/Mon/mabdel03/Tsai/Resilient/Counts/10101589
