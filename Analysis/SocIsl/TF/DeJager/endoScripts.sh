#!/bin/bash

#SBATCH -n 45                    # Number of cores requested
#SBATCH -t 48:00:00                # Runtime in hours
#SBATCH --mem=600G              # GB memory needed (memory PER CORE)
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-user=nkhera@college.harvard.edu
#SBATCH --mail-type=ALL

source /om2/user/mabdel03/anaconda/etc/profile.d/conda.sh

conda activate /om2/user/mabdel03/conda_envs/nebulaAnalysis7

cd cd /net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/Tsai/Preprocessing/Preprocessed_Counts/Resilient/
Rscript endoKV.Rscript
Rscript endoRV.Rscript
Rscript endoNV.Rscript


cd /om/scratch/Wed/mabdel03/Subfolder

Rscript endoKVD.Rscript
Rscript endoRVD.Rscript
Rscript endoNVD.Rscript

