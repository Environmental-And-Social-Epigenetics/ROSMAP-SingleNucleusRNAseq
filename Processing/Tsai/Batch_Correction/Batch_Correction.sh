#!/bin/bash  
#SBATCH -t 48:00:00
#SBATCH -n 128
#SBATCH --mem=512G
#SBATCH --gres=gpu:a100:1
#SBATCH --mail-user=mabdel03@mit.edu
#SBATCH --mail-type=BEGIN,END,FAIL

source /orcd/data/lhtsai/001/om2/mabdel03/miniforge3/etc/profile.d/conda.sh

conda activate /orcd/data/lhtsai/001/om2/mabdel03/conda_envs/BatchCorrection_SingleCell

python /orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Scripts/pipeline_scripts/Batch_Correction/Batch_Correction.py
