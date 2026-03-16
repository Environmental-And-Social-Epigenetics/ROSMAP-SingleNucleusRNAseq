#!/bin/bash

#SBATCH -n 45                    # Number of cores requested
#SBATCH -t 2:00:00                # Runtime in hours
#SBATCH --mem=700G              # GB memory needed (memory PER CORE)
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-user=nkhera@college.harvard.edu
#SBATCH --mail-type=ALL

source /om2/user/mabdel03/anaconda/etc/profile.d/conda.sh

conda activate /om2/user/mabdel03/conda_envs/nebulaAnalysis7

cd /om/scratch/Mon/mabdel03/SocialIsolation/Tsai/

Rscript tsaiGseaResults.Rscript

# cd /om/scratch/Mon/mabdel03/SocialIsolation/

# Rscript gseaResults.Rscript