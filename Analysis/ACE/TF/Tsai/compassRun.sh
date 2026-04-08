#!/bin/bash
#SBATCH -n 40
#SBATCH -t 24:00:00
#SBATCH --mem=600G
#SBATCH -o compass_%j.out
#SBATCH -e compass_%j.err
#SBATCH -p pi_lhtsai,pi_manoli

# ACE COMPASS metabolic flux analysis — Tsai cohort
# Usage: sbatch compassRun.sh <celltype> <sex>
# Example: sbatch compassRun.sh Exc male

set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  SCRIPT_DIR="${SLURM_SUBMIT_DIR}"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

CELLTYPE="${1:?ERROR: cell type argument required (e.g. Exc, Inh, Ast, Mic, Oli, OPC)}"
SEX="${2:?ERROR: sex argument required (male or female)}"

activate_env "${COMPASS_ANALYSIS_ENV}"

export CPLEX_STUDIO_DIR="${CPLEX_DIR}"
export PATH="${CPLEX_STUDIO_DIR}/cplex/python/3.9/x86-64_linux:${PATH}"
export PYTHONPATH="${CPLEX_STUDIO_DIR}/cplex/python/3.9/x86-64_linux:${PYTHONPATH:-}"

OUTPUT_DIR="${ACE_OUTPUT_ROOT}/TF/Tsai"
mkdir -p "${OUTPUT_DIR}"
cd "${OUTPUT_DIR}"

INPUT_TSV="${SEX}_matrix${CELLTYPE}.tsv"
OUTPUT_SUBDIR="${SEX^}${CELLTYPE}Compass"

if [[ ! -f "${INPUT_TSV}" ]]; then
  echo "ERROR: Input TSV not found: ${OUTPUT_DIR}/${INPUT_TSV}"
  echo "Run the SCENIC micropooling step first to generate expression matrices."
  exit 1
fi

echo "Running COMPASS for ${CELLTYPE} (${SEX})"
echo "Input: ${INPUT_TSV}"
echo "Output: ${OUTPUT_SUBDIR}/"

echo "y" | compass --data "${INPUT_TSV}" \
    --num-processes 10 \
    --species homo_sapiens \
    --output-dir "${OUTPUT_SUBDIR}"
