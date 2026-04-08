#!/bin/bash
#SBATCH -n 40                  	 # Number of cores requested
#SBATCH -t 4:00:00                # Runtime in hours
#SBATCH --mem=800G              # GB memory needed (memory PER CORE)
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-user=__SET_YOUR_EMAIL__
#SBATCH --mail-type=ALL

# Source central configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/paths.sh"
init_conda

cd "${SOCISL_OUTPUT_ROOT}"

activate_env "${BATCHCORR_ENV}"

export HDF5_USE_FILE_LOCKING=FALSE

python3 dejagerOverlapTry.py

