#!/bin/bash
#SBATCH -n 40                  	 # Number of cores requested
#SBATCH -t 20:00:00                # Runtime in hours
#SBATCH --mem=700G              # GB memory needed (memory PER CORE)
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-user=nkhera@college.harvard.edu
#SBATCH --mail-type=ALL

source /orcd/data/lhtsai/001/om2/mabdel03/miniforge3/etc/profile.d/conda.sh

cd /net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Scripts/pipeline_scripts/Batch_Correction

# conda activate /orcd/data/lhtsai/001/om2/mabdel03/conda_envs/BatchCorrection_SingleCell
conda activate /orcd/data/lhtsai/001/om2/mabdel03/conda_envs/NormFeatDimRedClust_SingleCell2

python3 plotBeyond.py
