#!/bin/bash
#SBATCH -n 40                  	 # Number of cores requested
#SBATCH -t 5:00:00                # Runtime in hours
#SBATCH --mem=500G              # GB memory needed (memory PER CORE)
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-user=__SET_YOUR_EMAIL__
#SBATCH --mail-type=ALL

# Source central configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/paths.sh"
init_conda

cd "${SOCISL_OUTPUT_ROOT}/Tsai"

activate_env "${NEBULA_ENV}"

export HDF5_USE_FILE_LOCKING=FALSE

Rscript socIslDegT.Rscript

cd "${SOCISL_OUTPUT_ROOT}"

sbatch gsea.sh

# cd /net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/Tsai/Preprocessing/Preprocessed_Counts/Resilient/

# sbatch plotPatients.sh


# conda activate /om2/user/mabdel03/conda_envs/NormFeatDimRedClust_SingleCell2

