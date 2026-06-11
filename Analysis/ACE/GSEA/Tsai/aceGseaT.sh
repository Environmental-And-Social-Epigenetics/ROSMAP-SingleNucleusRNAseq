#!/bin/bash
# ACE GSEA launcher for Tsai cohort
# Submits SLURM jobs for GSEA enrichment across all phenotypes and sexes

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

INTEGRATION="${1:-derived_batch}"
OUTPUT_ROOT="${ACE_OUTPUT_ROOT}/GSEA/Tsai"
DEG_ROOT="${ACE_OUTPUT_ROOT}/DEG/Tsai/results_${INTEGRATION}"
LOG_DIR="${OUTPUT_ROOT}/logs"
mkdir -p "${LOG_DIR}"

PHENOTYPES=("tot_adverse_exp" "early_hh_ses" "ace_aggregate")
SEXES=("Female" "Male")

echo "=== ACE GSEA Pipeline ==="
echo "Integration: ${INTEGRATION}"
echo "DEG results: ${DEG_ROOT}"
echo "Output:      ${OUTPUT_ROOT}"
echo ""

for PHENO in "${PHENOTYPES[@]}"; do
  for SEX in "${SEXES[@]}"; do
    JOB=$(sbatch --parsable --export=ALL,REPO_ROOT="${REPO_ROOT}",LAUNCHER_SCRIPT_DIR="${SCRIPT_DIR}",ACE_GSEA_SKIP_EXISTING="${ACE_GSEA_SKIP_EXISTING:-0}" \
      -p pi_lhtsai,pi_manoli,ou_bcs_low,mit_normal \
      -n 8 --mem=64G -t 12:00:00 \
      -o "${LOG_DIR}/%j_${PHENO}_${SEX}.out" \
      -e "${LOG_DIR}/%j_${PHENO}_${SEX}.err" \
      --job-name="ace_gsea_${PHENO}_${SEX}" \
      "${SCRIPT_DIR}/run_enrichment.sh" \
        "${INTEGRATION}" "${PHENO}" "${SEX}" "${DEG_ROOT}" "${OUTPUT_ROOT}")
    echo "GSEA ${PHENO} ${SEX}: job ${JOB}"
  done
done

echo ""
echo "Submitted $(( ${#PHENOTYPES[@]} * ${#SEXES[@]} )) GSEA jobs."
echo "Monitor with: squeue -u \$USER"
