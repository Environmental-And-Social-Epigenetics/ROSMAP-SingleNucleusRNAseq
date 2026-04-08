# DEPRECATED: Legacy data preparation script from Openmind era.
# The current Processing pipeline (Processing/*/Pipeline/) supersedes this.
# Paths may still reference decommissioned Openmind directories.
#
#!/bin/bash
#SBATCH -n 16                    # Reduce cores to reduce thread/memory pressure
#SBATCH -t 15:00:00
#SBATCH --mem=400G              # Total memory (not per core!)
#SBATCH -o tsaiAdataScenic.out
#SBATCH -e %j.err
#SBATCH --mail-user=nkhera@college.harvard.edu
#SBATCH --mail-type=ALL

source /om2/user/mabdel03/anaconda/etc/profile.d/conda.sh

conda activate /om2/user/mabdel03/conda_envs/scenicAnalysis

cd /om/scratch/Mon/mabdel03/SocialIsolation/Tsai

export HDF5_USE_FILE_LOCKING=FALSE

python3 getMatrices.py