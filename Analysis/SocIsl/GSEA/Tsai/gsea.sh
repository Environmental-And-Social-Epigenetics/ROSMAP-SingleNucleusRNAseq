#!/bin/bash

#SBATCH -n 45                    # Number of cores requested
#SBATCH -t 5:00:00                # Runtime in hours
#SBATCH --mem=100G              # GB memory needed (memory PER CORE)
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

cd "${SOCISL_OUTPUT_ROOT}"

Rscript gseaResults.Rscript