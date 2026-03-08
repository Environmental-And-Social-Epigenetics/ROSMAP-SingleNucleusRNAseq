#!/bin/bash
#SBATCH -n 40                  	 # Number of cores requested
#SBATCH -t 3:00:00                # Runtime in hours
#SBATCH --mem=500G              # GB memory needed (memory PER CORE)
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-user=nkhera@college.harvard.edu
#SBATCH --mail-type=ALL

source /orcd/data/lhtsai/001/om2/mabdel03/miniforge3/etc/profile.d/conda.sh

#REPLACE WITH NAME OF DIRECTORY YOUR FILE IS LOCATED IN

cd /net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Scripts/pipeline_scripts/Batch_Correction

conda activate /orcd/data/lhtsai/001/om2/mabdel03/conda_envs/BatchCorrection_SingleCell

# wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz

# gunzip gencode.v43.annotation.gtf.gz

python3 batchCluster.py

# conda activate /orcd/data/lhtsai/001/om2/mabdel03/conda_envs/NormFeatDimRedClust_SingleCell2

