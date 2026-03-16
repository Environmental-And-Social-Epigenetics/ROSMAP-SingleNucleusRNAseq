#!/bin/bash
#SBATCH -n 40                  	 # Number of cores requested
#SBATCH -t 24:00:00                # Runtime in hours
#SBATCH --mem=600G              # GB memory needed (memory PER CORE)
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-user=nkhera@college.harvard.edu
#SBATCH --mail-type=ALL

source /om2/user/mabdel03/anaconda/etc/profile.d/conda.sh

cd /om/scratch/Mon/mabdel03/SocialIsolation/Tsai

#python ver = 3.9

#python3 -m pip install git+https://github.com/yoseflab/Compass.git --upgrade

conda activate /om2/user/mabdel03/conda_envs/cplex_env2

export CPLEX_STUDIO_DIR=/om/scratch/Mon/mabdel03/SocialIsolation/opt/ibm/ILOG/CPLEX_Studio2211
export PATH=$CPLEX_STUDIO_DIR/cplex/python/3.9/x86-64_linux:$PATH
export PYTHONPATH=$CPLEX_STUDIO_DIR/cplex/python/3.9/x86-64_linux:$PYTHONPATH

echo "y" | compass --data matrix500_male_Inh.tsv --num-processes 40 --species homo_sapiens --output-dir CompassPMInhNew

echo "y" | compass --data matrix500_female_Inh.tsv --num-processes 40 --species homo_sapiens --output-dir CompassPFInhNew