#!/bin/bash
#SBATCH -n 40                  	 # Number of cores requested
#SBATCH -t 24:00:00                # Runtime in hours
#SBATCH --mem=600G              # GB memory needed (memory PER CORE)
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-user=__SET_YOUR_EMAIL__
#SBATCH --mail-type=ALL

# Source central configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/paths.sh"
init_conda

cd "${SOCISL_OUTPUT_ROOT}/Tsai"

#python ver = 3.9

#python3 -m pip install git+https://github.com/yoseflab/Compass.git --upgrade

activate_env "${COMPASS_ANALYSIS_ENV}"

export CPLEX_STUDIO_DIR="${CPLEX_DIR}"
export PATH=$CPLEX_STUDIO_DIR/cplex/python/3.9/x86-64_linux:$PATH
export PYTHONPATH=$CPLEX_STUDIO_DIR/cplex/python/3.9/x86-64_linux:$PYTHONPATH

echo "y" | compass --data female_matrixMic.tsv --num-processes 10 --species homo_sapiens --output-dir FemaleMicCompass
