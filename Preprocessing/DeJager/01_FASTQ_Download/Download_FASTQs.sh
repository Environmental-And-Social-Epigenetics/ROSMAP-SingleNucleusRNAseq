#!/bin/bash
#SBATCH -t 47:00:00
#SBATCH -n 64
#SBATCH --mem=512G
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-user=mabdel03@mit.edu
#SBATCH --mail-type=BEGIN,END,FAIL

source /orcd/data/lhtsai/001/om2/mabdel03/miniforge3/etc/profile.d/conda.sh

conda activate /orcd/data/lhtsai/001/om2/mabdel03/conda_envs/synapse_env

python /orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/DeJager/FASTQs_Download/FASTQ_Downloads_Scripts/Download_FASTQs.py



