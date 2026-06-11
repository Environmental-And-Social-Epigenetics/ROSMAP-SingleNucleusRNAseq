#!/bin/bash
#SBATCH -p pi_lhtsai,pi_manoli,ou_bcs_low,mit_preemptable
#SBATCH -n 20
#SBATCH --mem=200G
#SBATCH -t 24:00:00
#SBATCH -o %j.out
#SBATCH -e %j.err
#
# DeJager ACE SCENIC SLURM batch script for ONE cell-type x sex.
#
# Thin wrapper: the analysis method lives ENTIRELY in the Tsai modular pipeline
# (Analysis/ACE/SCENIC/Tsai/run_scenic.sh + scenic_analysis.py). This wrapper
# just forwards to it with cohort=dejager so the SAME code, SAME pool_size, and
# SAME stats run on the DeJager cohort. The Tsai run_scenic.sh derives the
# cohort-specific DEG input dir and SCENIC output dir from the cohort flag:
#   inputs  -> ${ACE_OUTPUT_ROOT}/DEG/DeJager/celltype_splits_${INTEGRATION}/${CELL_TYPE}.h5ad
#   outputs -> ${ACE_OUTPUT_ROOT}/SCENIC/DeJager/results_${INTEGRATION}/${PHENOTYPE}/${SEX}_${CELL_TYPE}
#
# Usage (normally invoked by aceScenicDJ.sh):
#   bash run_scenic.sh INTEGRATION PHENOTYPE CELL_TYPE SEX

set -euo pipefail

# ---------------------------------------------------------------------------
# Resolve this script's directory (works under SLURM and direct invocation)
# ---------------------------------------------------------------------------
if [[ -n "${LAUNCHER_SCRIPT_DIR:-}" ]]; then
  SCRIPT_DIR="${LAUNCHER_SCRIPT_DIR}"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  SCRIPT_DIR="${SLURM_SUBMIT_DIR}"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

# Tsai modular run_scenic.sh is the single source of truth.
TSAI_RUN_SCENIC="$(cd "${SCRIPT_DIR}/../Tsai" && pwd)/run_scenic.sh"
if [[ ! -f "${TSAI_RUN_SCENIC}" ]]; then
  echo "ERROR: Tsai run_scenic.sh not found at ${TSAI_RUN_SCENIC}" >&2
  exit 1
fi

# ---------------------------------------------------------------------------
# Arguments (passed by aceScenicDJ.sh or manually)
# ---------------------------------------------------------------------------
INTEGRATION="${1:?ERROR: integration argument required (e.g. library_id)}"
PHENOTYPE="${2:?ERROR: phenotype argument required (e.g. tot_adverse_exp)}"
CELL_TYPE="${3:?ERROR: cell_type argument required (e.g. Mic)}"
SEX="${4:?ERROR: sex argument required (Male or Female)}"

# Pin the Tsai run_scenic.sh self-resolution to the Tsai directory so its
# REPO_ROOT and scenic_analysis.py lookups are correct regardless of how this
# wrapper was launched.
export LAUNCHER_SCRIPT_DIR="$(cd "$(dirname "${TSAI_RUN_SCENIC}")" && pwd)"

# Delegate to the shared modular pipeline with cohort=dejager.
exec bash "${TSAI_RUN_SCENIC}" "${INTEGRATION}" "${PHENOTYPE}" "${CELL_TYPE}" "${SEX}" "dejager"
