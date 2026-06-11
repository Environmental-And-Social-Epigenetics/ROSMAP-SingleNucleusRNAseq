#!/bin/bash
# =============================================================================
# ACE hdWGCNA Co-Expression Network Launcher — Tsai Cohort
#
# Submits SLURM jobs to run hdWGCNA on priority cell types (those with
# the most DEGs). Each cell type × sex is an independent job.
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

INTEGRATION="${1:-derived_batch}"
PHENOTYPE="${2:-tot_adverse_exp}"

INPUT_DIR="${ACE_OUTPUT_ROOT}/DEG/Tsai/celltype_splits_${INTEGRATION}"
OUTPUT_ROOT="${ACE_OUTPUT_ROOT}/hdWGCNA/Tsai"
LOG_DIR="${OUTPUT_ROOT}/logs"
mkdir -p "${LOG_DIR}"

# Priority cell types by DEG count. Broad excitatory/inhibitory split files
# are named broad_Exc/broad_Inh by the DEG prep step.
MALE_CTS=("broad_Inh" "Mic" "Ast" "In-PV_Basket" "Oli" "broad_Exc")
FEMALE_CTS=("Oli" "Ex-L2_3" "Ast" "broad_Exc" "OPC")

if [[ "${PHENOTYPE}" == "all" ]]; then
    PHENOTYPES=("tot_adverse_exp" "early_hh_ses" "ace_aggregate")
else
    PHENOTYPES=("${PHENOTYPE}")
fi

echo "=== ACE hdWGCNA Pipeline ==="
echo "Integration: ${INTEGRATION}"
echo "Phenotypes:  ${PHENOTYPES[*]}"
echo "Input:       ${INPUT_DIR}"
echo "Output:      ${OUTPUT_ROOT}"
echo ""

submit_job() {
    local PHENO="$1"
    local SEX="$2"
    local CT="$3"

    local H5AD="${INPUT_DIR}/${CT}.h5ad"
    if [[ ! -f "${H5AD}" ]]; then
        echo "  SKIP: ${H5AD} not found"
        return
    fi

    local JOB=$(sbatch --parsable --export=ALL,REPO_ROOT="${REPO_ROOT}",LAUNCHER_SCRIPT_DIR="${SCRIPT_DIR}" \
        -p pi_lhtsai,pi_manoli \
        -n 16 --mem=300G -t 24:00:00 \
        -o "${LOG_DIR}/%j_${PHENO}_${SEX}_${CT}.out" \
        -e "${LOG_DIR}/%j_${PHENO}_${SEX}_${CT}.err" \
        --job-name="ace_wgcna_${PHENO}_${SEX}_${CT}" \
        "${SCRIPT_DIR}/run_wgcna.sh" \
            "${INTEGRATION}" "${PHENO}" "${SEX}" "${CT}" \
            "${H5AD}" "${OUTPUT_ROOT}")
    echo "  hdWGCNA ${PHENO} ${SEX} ${CT}: job ${JOB}"
}

for PHENO in "${PHENOTYPES[@]}"; do
    echo "--- Phenotype: ${PHENO} ---"
    echo "Male cell types:"
    for CT in "${MALE_CTS[@]}"; do
        submit_job "${PHENO}" "Male" "${CT}"
    done

    echo "Female cell types:"
    for CT in "${FEMALE_CTS[@]}"; do
        submit_job "${PHENO}" "Female" "${CT}"
    done
done

echo ""
echo "All hdWGCNA jobs submitted."
