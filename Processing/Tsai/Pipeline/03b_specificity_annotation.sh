#!/bin/bash
#SBATCH -J tsai_annot3b
#SBATCH -p pi_lhtsai,pi_manoli
#SBATCH -t 12:00:00
#SBATCH -n 16
#SBATCH --mem=500G
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

# REPO_ROOT may be pre-exported by a submit wrapper. Otherwise derive it.
# Under SLURM, BASH_SOURCE points at the spooled copy under /var, so prefer
# SLURM_SUBMIT_DIR; fall back to BASH_SOURCE for standalone (non-SLURM) use.
if [[ -z "${REPO_ROOT:-}" ]]; then
    if [[ -n "${SLURM_SUBMIT_DIR:-}" && -f "${SLURM_SUBMIT_DIR}/03b_specificity_annotation.py" ]]; then
        _SCRIPT_DIR="${SLURM_SUBMIT_DIR}"
    else
        _SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    fi
    REPO_ROOT="$(cd "${_SCRIPT_DIR}/../../.." && pwd)"
fi
if [[ ! -f "${REPO_ROOT}/config/paths.sh" ]]; then
    echo "ERROR: could not locate repo root (REPO_ROOT=${REPO_ROOT}). Export REPO_ROOT or submit from Processing/Tsai/Pipeline." >&2
    exit 1
fi
SCRIPT_DIR="${REPO_ROOT}/Processing/Tsai/Pipeline"

source "${REPO_ROOT}/config/paths.sh"

mkdir -p "${TSAI_PROCESSING_LOGS}"

export HDF5_USE_FILE_LOCKING=FALSE
export PYTHONPATH="${REPO_ROOT}/src:${PYTHONPATH:-}"

# Relax nounset for conda activation (activate.d scripts reference unset vars).
set +u
source "${CONDA_INIT_SCRIPT}"
conda activate "${BATCHCORR_ENV}"
set -u
if [[ -z "${CONDA_PREFIX:-}" ]]; then
    echo "ERROR: Failed to activate conda environment: ${BATCHCORR_ENV}"
    exit 1
fi

# Pass-through knobs (all optional; sensible defaults in the script):
#   TOP_N, TAU_MIN, MIN_DETECT, PRIMARY_SCORE, SUBSAMPLE, EXTRA_ARGS
EXTRA=()
[[ -n "${TOP_N:-}" ]]         && EXTRA+=(--top-n "${TOP_N}")
[[ -n "${TAU_MIN:-}" ]]       && EXTRA+=(--tau-min "${TAU_MIN}")
[[ -n "${MIN_DETECT:-}" ]]    && EXTRA+=(--min-detect "${MIN_DETECT}")
[[ -n "${PRIMARY_SCORE:-}" ]] && EXTRA+=(--primary-score "${PRIMARY_SCORE}")
[[ -n "${SUBSAMPLE:-}" ]]     && EXTRA+=(--subsample "${SUBSAMPLE}")
[[ "${TRIM_ONLY:-}" == "1" ]] && EXTRA+=(--trim-only)
[[ "${ALLOW_TAU_FALLBACK:-}" == "1" ]] && EXTRA+=(--allow-tau-fallback)
# shellcheck disable=SC2206
[[ -n "${EXTRA_ARGS:-}" ]]    && EXTRA+=(${EXTRA_ARGS})

echo "REPO_ROOT:    ${REPO_ROOT}"
echo "BATCHCORR_ENV:${BATCHCORR_ENV}"
echo "TSAI_INTEGRATED (configured): ${TSAI_INTEGRATED}"
echo "Extra args:   ${EXTRA[*]:-(none)}"

python "${SCRIPT_DIR}/03b_specificity_annotation.py" "${EXTRA[@]}"
