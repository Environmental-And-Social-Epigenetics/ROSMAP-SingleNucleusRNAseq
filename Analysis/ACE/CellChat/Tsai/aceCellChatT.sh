#!/bin/bash
# =============================================================================
# ACE Cell-Cell Communication Launcher — Tsai Cohort
#
# Submits SLURM jobs to run CellChat differential interaction analysis
# between ACE-high and ACE-low groups, stratified by sex.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

INTEGRATION="${1:-derived_batch}"
PHENOTYPE="${2:-all}"

# Resolve input h5ad based on integration
if [[ "${INTEGRATION}" == "derived_batch" ]]; then
    INPUT_H5AD="${TSAI_INTEGRATED}/tsai_annotated.h5ad"
elif [[ "${INTEGRATION}" == "projid" ]]; then
    INPUT_H5AD="${TSAI_INTEGRATED_PROJID}/tsai_annotated.h5ad"
else
    echo "ERROR: Unknown integration '${INTEGRATION}'. Use 'derived_batch' or 'projid'." >&2
    exit 1
fi

OUTPUT_ROOT="${ACE_OUTPUT_ROOT}/CellChat/Tsai"
LOG_DIR="${OUTPUT_ROOT}/logs"
mkdir -p "${LOG_DIR}"

if [[ "${PHENOTYPE}" == "all" ]]; then
    PHENOTYPES=("tot_adverse_exp" "early_hh_ses" "ace_aggregate")
else
    PHENOTYPES=("${PHENOTYPE}")
fi
SEXES=("Male" "Female")

echo "=== ACE CellChat Pipeline ==="
echo "Integration: ${INTEGRATION}"
echo "Phenotypes:  ${PHENOTYPES[*]}"
echo "Input h5ad:  ${INPUT_H5AD}"
echo "Output:      ${OUTPUT_ROOT}"
echo ""

n=0
for PHENO in "${PHENOTYPES[@]}"; do
    for SEX in "${SEXES[@]}"; do
        JOB=$(sbatch --parsable --export=ALL,REPO_ROOT="${REPO_ROOT}",LAUNCHER_SCRIPT_DIR="${SCRIPT_DIR}" \
            -p pi_lhtsai,pi_manoli \
            -n 16 --mem=200G -t 12:00:00 \
            -o "${LOG_DIR}/%j_cellchat_${PHENO}_${SEX}.out" \
            -e "${LOG_DIR}/%j_cellchat_${PHENO}_${SEX}.err" \
            --job-name="ace_cellchat_${PHENO}_${SEX}" \
            "${SCRIPT_DIR}/run_cellchat.sh" \
                "${INTEGRATION}" "${PHENO}" "${SEX}" "${INPUT_H5AD}" "${OUTPUT_ROOT}")
        echo "CellChat ${PHENO} ${SEX}: job ${JOB}"
        n=$((n + 1))
    done
done

echo ""
echo "Submitted ${n} CellChat jobs."
