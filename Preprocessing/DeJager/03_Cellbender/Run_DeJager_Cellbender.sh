#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

source "${REPO_ROOT}/config/paths.sh"

set +u
source "${CONDA_INIT_SCRIPT}"
conda activate "${PYTHON_ENV}"
set -u

python "${SCRIPT_DIR}/Run_DeJager_Cellbender.py" "$@"
