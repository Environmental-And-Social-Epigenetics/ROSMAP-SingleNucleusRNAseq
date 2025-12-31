#!/bin/bash
#SBATCH -t 47:00:00
#SBATCH -n 32
#SBATCH --mem=250G
#SBATCH --output=/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/Tsai/Preprocessing/Counts/SocIsl/Batch_Scripts_Outs/slurm-%j.out
#SBATCH --error=/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/Tsai/Preprocessing/Counts/SocIsl/Batch_Scripts_Err/slurm-%j.err
#SBATCH --mail-user=nkhera@college.harvard.edu
#SBATCH --mail-type=FAIL
export PATH=/om2/user/$USER/apps/yard/cellranger-8.0.0:$PATH
cellranger count --create-bam false --include-introns true --nosecondary --r1-length 26 --id D19-8345 --transcriptome=/orcd/data/lhtsai/001/om2/mabdel03/yard/references/human/refdata-gex-GRCh38-2020-A --sample D19-8345 --fastqs /om/scratch/Mon/mabdel03/SocIsl/FASTQs/20907534/D19-8345_1 --output-dir=/om/scratch/Mon/mabdel03/SocIsl/Counts/20907534
