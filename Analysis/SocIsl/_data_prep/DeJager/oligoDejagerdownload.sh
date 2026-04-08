# DEPRECATED: Legacy data preparation script from Openmind era.
# The current Processing pipeline (Processing/*/Pipeline/) supersedes this.
# Paths may still reference decommissioned Openmind directories.
#
#!/bin/bash
#SBATCH -n 50                    # Number of cores requested
#SBATCH -t 5:00:00               # Runtime in hours
#SBATCH --mem=400G               # GB memory needed (memory PER CORE)
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-user=nkhera@college.harvard.edu
#SBATCH --mail-type=ALL

source /om2/user/mabdel03/anaconda/etc/profile.d/conda.sh

# conda create -n myJupyterEnv python pip jupyterlab notebook
cd /om/scratch/Mon/mabdel03/SocialIsolation

conda activate /net/vast-storage/scratch/vast/lhtsai/mabdel03/conda_envs/bcftools_env

python3 oligoDownload.py