#!/bin/bash
#SBATCH -p pi_lhtsai,mit_normal
#SBATCH -n 4
#SBATCH -t 01:00:00
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

activate_env "${SCCOMP_ENV}"

Rscript "${SCRIPT_DIR}/sccomp_visualize.R" \
    --integration "${INTEGRATION}" \
    --output-root "${OUTPUT_ROOT}" \
    "$@"
