#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -n 16
#SBATCH --mem=64G
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

source "${REPO_ROOT}/config/paths.sh"

source "${CONDA_INIT_SCRIPT}"
conda activate "${PYTHON_ENV}"

python "${SCRIPT_DIR}/Count_DeJager.py"
