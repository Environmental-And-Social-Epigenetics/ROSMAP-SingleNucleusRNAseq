#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -n 16
#SBATCH --mem=64G
#SBATCH --output=/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/DeJager/Original_Counts/Counts_outs/slurm-%j.out
#SBATCH --error=/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/DeJager/Original_Counts/Counts_errors/slurm-%j.err
#SBATCH --mail-user=mabdel03@mit.edu
#SBATCH --mail-type=BEGIN,END,FAIL

source /orcd/data/lhtsai/001/om2/mabdel03/miniforge3/etc/profile.d/conda.sh

conda activate /orcd/data/lhtsai/001/om2/mabdel03/conda_envs/python_env

python /orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/DeJager/Original_Counts/Counts_Master_Scripts/Count_DeJager.py
