# DEPRECATED: Legacy data preparation script from Openmind era.
# The current Processing pipeline (Processing/*/Pipeline/) supersedes this.
# Paths may still reference decommissioned Openmind directories.
#
#!/bin/bash
#SBATCH -n 40                  	 # Number of cores requested
#SBATCH -t 10:00:00                # Runtime in hours
#SBATCH --mem=500G              # GB memory needed (memory PER CORE)
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-user=nkhera@college.harvard.edu
#SBATCH --mail-type=ALL

source /om2/user/mabdel03/anaconda/etc/profile.d/conda.sh

cd /om/scratch/Mon/mabdel03/SocialIsolation/Tsai

# conda activate /om2/user/mabdel03/conda_envs/nebulaAnalysis7

# export HDF5_USE_FILE_LOCKING=FALSE

# Rscript h5SeuratWork.Rscript

conda activate /om2/user/mabdel03/conda_envs/BatchCorrection_SingleCell

export HDF5_USE_FILE_LOCKING=FALSE

python3 vars.py

# sbatch micropool.sh

# sbatch newScenic.sh