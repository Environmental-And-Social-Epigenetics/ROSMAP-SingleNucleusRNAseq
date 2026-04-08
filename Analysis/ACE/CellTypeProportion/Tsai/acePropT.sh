#!/bin/bash
# Submit the ACE sccomp workflow for the Tsai cohort.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

OUTPUT_ROOT="${ACE_PROP_OUTPUT_ROOT:-${ANALYSIS_OUTPUT_ROOT}/ACE/CellTypeProportion/Tsai}"
LOG_DIR="${OUTPUT_ROOT}/logs"
mkdir -p "${LOG_DIR}"

if [[ $# -gt 0 ]]; then
    INTEGRATIONS=("$@")
else
    INTEGRATIONS=("derived_batch" "projid")
fi

SEX_STRATA=("all" "male" "female")
RESOLUTIONS=("fine" "broad")

echo "=== ACE Cell Type Proportion Pipeline: Tsai ==="
echo "Output root: ${OUTPUT_ROOT}"
echo "Integrations: ${INTEGRATIONS[*]}"
echo ""

for INTEGRATION in "${INTEGRATIONS[@]}"; do
    echo "--- Integration: ${INTEGRATION} ---"

    PREP_JOB=$(sbatch --parsable \
        --export=ALL,ACE_PROP_OUTPUT_ROOT="${OUTPUT_ROOT}" \
        -o "${LOG_DIR}/%j_prep_${INTEGRATION}.out" \
        -e "${LOG_DIR}/%j_prep_${INTEGRATION}.err" \
        --job-name="acepropT_prep_${INTEGRATION}" \
        "${SCRIPT_DIR}/run_prep.sh" "${INTEGRATION}")
    echo "  Preprocessing: job ${PREP_JOB}"

    SCCOMP_JOBS=()
    for SEX in "${SEX_STRATA[@]}"; do
        for RESOLUTION in "${RESOLUTIONS[@]}"; do
            SCCOMP_JOB=$(sbatch --parsable \
                --dependency=afterok:${PREP_JOB} \
                --export=ALL,ACE_PROP_OUTPUT_ROOT="${OUTPUT_ROOT}" \
                -o "${LOG_DIR}/%j_sccomp_${INTEGRATION}_${SEX}_${RESOLUTION}.out" \
                -e "${LOG_DIR}/%j_sccomp_${INTEGRATION}_${SEX}_${RESOLUTION}.err" \
                --job-name="acepropT_${INTEGRATION}_${SEX}_${RESOLUTION}" \
                "${SCRIPT_DIR}/run_sccomp.sh" "${INTEGRATION}" "${SEX}" "${RESOLUTION}")
            SCCOMP_JOBS+=("${SCCOMP_JOB}")
            echo "  sccomp ${SEX}/${RESOLUTION}: job ${SCCOMP_JOB}"
        done
    done

    DEP_STRING=$(IFS=:; echo "${SCCOMP_JOBS[*]}")
    VIZ_JOB=$(sbatch --parsable \
        --dependency=afterok:${DEP_STRING} \
        --export=ALL,ACE_PROP_OUTPUT_ROOT="${OUTPUT_ROOT}" \
        -o "${LOG_DIR}/%j_viz_${INTEGRATION}.out" \
        -e "${LOG_DIR}/%j_viz_${INTEGRATION}.err" \
        --job-name="acepropT_viz_${INTEGRATION}" \
        "${SCRIPT_DIR}/run_visualize.sh" "${INTEGRATION}")
    echo "  Visualization: job ${VIZ_JOB}"
    echo ""
done

echo "All jobs submitted. Monitor with: squeue -u \$USER"
