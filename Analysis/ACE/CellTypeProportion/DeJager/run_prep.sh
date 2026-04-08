#!/bin/bash
#SBATCH -p pi_lhtsai,mit_normal
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH --mem=16G
#SBATCH -o %j.out
#SBATCH -e %j.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

INTEGRATION="${1:?ERROR: integration argument required (library_id, patient_id, pool_batch, or derived_batch)}"
shift || true

OUTPUT_ROOT="${ACE_PROP_OUTPUT_ROOT:-${ANALYSIS_OUTPUT_ROOT}/ACE/CellTypeProportion/DeJager}"

activate_env "${NEBULA_ENV}"
export HDF5_USE_FILE_LOCKING=FALSE

python "${SCRIPT_DIR}/prep_counts.py" \
    --integration "${INTEGRATION}" \
    --output-root "${OUTPUT_ROOT}" \
    "$@"
