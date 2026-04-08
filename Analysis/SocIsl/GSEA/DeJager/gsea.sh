#!/bin/bash

#SBATCH -n 45                    # Number of cores requested
#SBATCH -t 2:00:00                # Runtime in hours
#SBATCH --mem=700G              # GB memory needed (memory PER CORE)
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-user=__SET_YOUR_EMAIL__
#SBATCH --mail-type=ALL

# Source central configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/paths.sh"
init_conda

activate_env "${NEBULA_ENV}"

cd "${SOCISL_OUTPUT_ROOT}/Tsai"

Rscript tsaiGseaResults.Rscript

# cd /om/scratch/Mon/mabdel03/SocialIsolation/

# Rscript gseaResults.Rscript