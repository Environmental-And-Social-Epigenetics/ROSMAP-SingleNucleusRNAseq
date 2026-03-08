#!/bin/bash
#SBATCH -t 47:00:00
#SBATCH -n 64
#SBATCH --mem=512G
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

source "${REPO_ROOT}/config/paths.sh"

source "${CONDA_INIT_SCRIPT}"
conda activate "${SYNAPSE_ENV}"

python "${SCRIPT_DIR}/Download_FASTQs.py"
