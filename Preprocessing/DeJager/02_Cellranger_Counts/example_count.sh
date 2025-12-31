#!/bin/bash
#SBATCH -t 47:00:00
#SBATCH -n 32
#SBATCH --mem=128G
#SBATCH --output=/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/DeJager/Original_Counts/Counts_outs/slurm-%j.out
#SBATCH --error=/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/DeJager/Original_Counts/Counts_errors/slurm-%j.err
#SBATCH --mail-user=mabdel03@mit.edu
#SBATCH --mail-type=FAIL
export PATH=/orcd/data/lhtsai/001/om2/mabdel03/apps/yard/cellranger-8.0.0:$PATH
cellranger count --create-bam true --include-introns true --nosecondary --r1-length 26 --id 190403-B4-A --transcriptome=/orcd/data/lhtsai/001/om2/mabdel03/yard/references/human/refdata-gex-GRCh38-2020-A --sample 190403-B4-A_Broad --fastqs /om/scratch/Mon/mabdel03/FASTQs/190403-B4-A --output-dir=/om/scratch/Mon/mabdel03/Counts/190403-B4-A
