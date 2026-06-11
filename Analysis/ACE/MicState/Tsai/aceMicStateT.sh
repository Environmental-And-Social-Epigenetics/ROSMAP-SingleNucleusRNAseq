#!/bin/bash
# =============================================================================
# ACE Microglial State Analysis Launcher — Tsai Cohort
#
# Submits SLURM jobs to analyze microglial activation states in
# ACE-exposed individuals, stratified by sex.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

INTEGRATION="${1:-derived_batch}"
PHENOTYPE="${2:-all}"

INPUT_H5AD="${ACE_OUTPUT_ROOT}/DEG/Tsai/celltype_splits_${INTEGRATION}/Mic.h5ad"
OUTPUT_ROOT="${ACE_OUTPUT_ROOT}/MicState/Tsai"
LOG_DIR="${OUTPUT_ROOT}/logs"
mkdir -p "${LOG_DIR}"

if [[ ! -f "${INPUT_H5AD}" ]]; then
    echo "ERROR: Microglia h5ad not found: ${INPUT_H5AD}" >&2
    echo "Run the DEG prep step first: bash Analysis/ACE/DEG/Tsai/aceDegT.sh" >&2
    exit 1
fi

if [[ "${PHENOTYPE}" == "all" ]]; then
    PHENOTYPES=("tot_adverse_exp" "early_hh_ses" "ace_aggregate")
else
    PHENOTYPES=("${PHENOTYPE}")
fi
SEXES=("Male" "Female")

echo "=== ACE Microglial State Pipeline ==="
echo "Integration: ${INTEGRATION}"
echo "Phenotypes:  ${PHENOTYPES[*]}"
echo "Input:       ${INPUT_H5AD}"
echo "Output:      ${OUTPUT_ROOT}"
echo ""

n=0
for PHENO in "${PHENOTYPES[@]}"; do
    for SEX in "${SEXES[@]}"; do
        JOB=$(sbatch --parsable --export=ALL,REPO_ROOT="${REPO_ROOT}",LAUNCHER_SCRIPT_DIR="${SCRIPT_DIR}" \
            -p pi_lhtsai,pi_manoli,ou_bcs_low,mit_normal \
            -n 8 --mem=128G -t 6:00:00 \
            -o "${LOG_DIR}/%j_micstate_${PHENO}_${SEX}.out" \
            -e "${LOG_DIR}/%j_micstate_${PHENO}_${SEX}.err" \
            --job-name="ace_micstate_${PHENO}_${SEX}" \
            "${SCRIPT_DIR}/run_mic_state.sh" \
                "${INTEGRATION}" "${PHENO}" "${SEX}" "${INPUT_H5AD}" \
                "${OUTPUT_ROOT}")
        echo "MicState ${PHENO} ${SEX}: job ${JOB}"
        n=$((n + 1))
    done
done

echo ""
echo "Submitted ${n} microglial state analysis jobs."
