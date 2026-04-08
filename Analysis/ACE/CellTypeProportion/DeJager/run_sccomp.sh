#!/bin/bash
#SBATCH -p pi_lhtsai,mit_normal
#SBATCH -n 8
#SBATCH -t 08:00:00
#SBATCH --mem=64G
#SBATCH -o %j.out
#SBATCH -e %j.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

INTEGRATION="${1:?ERROR: integration argument required (library_id, patient_id, pool_batch, or derived_batch)}"
SEX="${2:?ERROR: sex argument required (all, male, or female)}"
RESOLUTION="${3:?ERROR: resolution argument required (fine or broad)}"
shift 3 || true

OUTPUT_ROOT="${ACE_PROP_OUTPUT_ROOT:-${ANALYSIS_OUTPUT_ROOT}/ACE/CellTypeProportion/DeJager}"

activate_env "${SCCOMP_ENV}"
export HDF5_USE_FILE_LOCKING=FALSE

Rscript "${SCRIPT_DIR}/sccomp_analysis.R" \
    --integration "${INTEGRATION}" \
    --sex "${SEX}" \
    --resolution "${RESOLUTION}" \
    --output-root "${OUTPUT_ROOT}" \
    "$@"
