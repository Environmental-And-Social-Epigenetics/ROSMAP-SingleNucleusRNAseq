#!/bin/bash
# =============================================================================
# ACE Epigenomic Integration Launcher — Tsai Cohort
#
# Submits SLURM jobs that:
#   1. (optional) call download_epigenomic.sh if data is not yet on disk
#   2. call match_individuals.py to produce individual_overlap.csv
#   3. call integrate.py per phenotype (DEG x epigenetic mark joins)
#
# Usage:
#   bash aceEpiIntT.sh [INTEGRATION] [PHENOTYPE]
#
# Defaults:
#   INTEGRATION = derived_batch
#   PHENOTYPE   = all (loops over tot_adverse_exp, early_hh_ses, ace_aggregate)
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

INTEGRATION="${1:-derived_batch}"
PHENOTYPE="${2:-all}"

OUTPUT_ROOT="${ACE_OUTPUT_ROOT}/EpigenomicIntegration/Tsai"
LOG_DIR="${OUTPUT_ROOT}/logs"
EPI_DATA_DIR="${OUTPUT_ROOT}/epigenomic_data"
mkdir -p "${LOG_DIR}"

if [[ "${PHENOTYPE}" == "all" ]]; then
    PHENOTYPES=("tot_adverse_exp" "early_hh_ses" "ace_aggregate")
else
    PHENOTYPES=("${PHENOTYPE}")
fi

CONSORTIUM_PY="${CONSORTIUM_ENV:-${CONDA_ENV_BASE}/consortium}/bin/python"
NEBULA_RSCRIPT="${NEBULA_ENV:-${CONDA_ENV_BASE}/nebulaAnalysis7}/bin/Rscript"

# ----------------------------------------------------------------------------
# Step 1: download (skipped if already present)
# ----------------------------------------------------------------------------
NEED_DOWNLOAD=1
if [[ -d "${EPI_DATA_DIR}/h3k9ac" || -d "${EPI_DATA_DIR}/methylation" ]]; then
    NEED_DOWNLOAD=0
fi

DOWNLOAD_JOB=""
if [[ "${NEED_DOWNLOAD}" == "1" ]]; then
    if [[ -z "${SYNAPSE_AUTH_TOKEN:-}" ]]; then
        echo "ERROR: epigenomic data not on disk and SYNAPSE_AUTH_TOKEN is not set."
        echo "       export SYNAPSE_AUTH_TOKEN=<your_pat> and rerun."
        exit 1
    fi
    DOWNLOAD_JOB=$(sbatch --parsable --export=ALL,SYNAPSE_AUTH_TOKEN \
        -p pi_lhtsai,pi_manoli \
        -n 4 --mem=16G -t 12:00:00 \
        -o "${LOG_DIR}/%j_download.out" \
        -e "${LOG_DIR}/%j_download.err" \
        --job-name="ace_epi_download" \
        --wrap "bash ${SCRIPT_DIR}/download_epigenomic.sh")
    echo "Download job: ${DOWNLOAD_JOB}"
fi

# ----------------------------------------------------------------------------
# Step 2: match individuals
# ----------------------------------------------------------------------------
OVERLAP_CSV="${OUTPUT_ROOT}/individual_overlap.csv"
DEPENDS_FOR_MATCH=""
if [[ -n "${DOWNLOAD_JOB}" ]]; then
    DEPENDS_FOR_MATCH="--dependency=afterok:${DOWNLOAD_JOB}"
fi
MATCH_JOB=$(sbatch --parsable --export=ALL ${DEPENDS_FOR_MATCH} \
    -p pi_lhtsai,pi_manoli \
    -n 2 --mem=16G -t 1:00:00 \
    -o "${LOG_DIR}/%j_match.out" \
    -e "${LOG_DIR}/%j_match.err" \
    --job-name="ace_epi_match" \
    --wrap "set -euo pipefail; \
        '${CONSORTIUM_PY}' ${SCRIPT_DIR}/match_individuals.py \
            --tsai-h5ad '${TSAI_INTEGRATED}/tsai_annotated.h5ad' \
            --pheno-csv '${ACE_SCORES_CSV}' \
            --epi-data-dir '${EPI_DATA_DIR}' \
            --output-csv '${OVERLAP_CSV}'")
echo "Match job:    ${MATCH_JOB}"

# ----------------------------------------------------------------------------
# Step 3: integrate per phenotype
# ----------------------------------------------------------------------------
for PHENO in "${PHENOTYPES[@]}"; do
    DEG_PHENO_DIR="${ACE_OUTPUT_ROOT}/DEG/Tsai/results_${INTEGRATION}/${PHENO}"
    OUT_PHENO_DIR="${OUTPUT_ROOT}/${PHENO}"
    mkdir -p "${OUT_PHENO_DIR}"

    JOB=$(sbatch --parsable --export=ALL --dependency=afterok:${MATCH_JOB} \
        -p pi_lhtsai,pi_manoli \
        -n 4 --mem=32G -t 4:00:00 \
        -o "${LOG_DIR}/%j_integrate_${PHENO}.out" \
        -e "${LOG_DIR}/%j_integrate_${PHENO}.err" \
        --job-name="ace_epi_${PHENO}" \
        --wrap "set -euo pipefail; \
            '${CONSORTIUM_PY}' ${SCRIPT_DIR}/integrate.py \
                --phenotype '${PHENO}' \
                --deg-dir '${DEG_PHENO_DIR}' \
                --epi-data-dir '${EPI_DATA_DIR}' \
                --overlap-csv '${OVERLAP_CSV}' \
                --pheno-csv '${ACE_SCORES_CSV}' \
                --output-dir '${OUT_PHENO_DIR}' \
                --rscript '${NEBULA_RSCRIPT}'")
    echo "Integrate ${PHENO}: ${JOB}"
done

echo ""
echo "Submitted Epigenomic Integration pipeline."
echo "Steps chain via SLURM dependencies; final outputs land under:"
echo "  ${OUTPUT_ROOT}/<phenotype>/deg_x_epigenetic.tsv"
echo "  ${OUTPUT_ROOT}/<phenotype>/summary_by_celltype.tsv"
