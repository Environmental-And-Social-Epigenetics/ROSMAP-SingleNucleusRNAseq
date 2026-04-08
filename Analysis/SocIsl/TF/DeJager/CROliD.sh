# Source central configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../../../config/paths.sh"
init_conda

cd "${SOCISL_OUTPUT_ROOT}"

#python ver = 3.9

#python3 -m pip install git+https://github.com/yoseflab/Compass.git --upgrade

activate_env "${COMPASS_ANALYSIS_ENV}"

export CPLEX_STUDIO_DIR="${CPLEX_DIR}"
export PATH=$CPLEX_STUDIO_DIR/cplex/python/3.9/x86-64_linux:$PATH
export PYTHONPATH=$CPLEX_STUDIO_DIR/cplex/python/3.9/x86-64_linux:$PYTHONPATH

echo "y" | compass --data dejpseudobulk_male_Oli.tsv --num-processes 10 --species homo_sapiens --output-dir CompassPMOli

echo "y" | compass --data dejpseudobulk_female_Oli.tsv --num-processes 10 --species homo_sapiens --output-dir CompassPFOli

./CROPCD.sh