#!/bin/bash
#SBATCH -n 40                    # Reduce cores to reduce thread/memory pressure
#SBATCH -t 48:00:00
#SBATCH --mem=600G              # Total memory (not per core!)
#SBATCH -o newAdataScenic.out
#SBATCH -e %j.err
#SBATCH --mail-user=__SET_YOUR_EMAIL__
#SBATCH --mail-type=ALL

# Source central configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/paths.sh"
init_conda

activate_env "${SCENIC_ANALYSIS_ENV}"

cd "${SOCISL_OUTPUT_ROOT}/Tsai"

export HDF5_USE_FILE_LOCKING=FALSE

python3 newScenic.py