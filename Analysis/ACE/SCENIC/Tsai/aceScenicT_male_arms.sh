#!/bin/bash
# SCENIC for the non-ANCOVA male ACE AD-model arms, build-once / associate-per-arm.
#
# Stage 1 (heavy, phenotype-independent): GRNBoost -> cisTarget -> AUCell, once per
#   signal-bearing cell type (males). Reuses existing auc_matrix.csv where present
#   (only missing cell types are (re)built). Output dir: the canonical
#   results_derived_batch/tot_adverse_exp/Male_<CT>/.
# Stage 2 (cheap): per-arm OLS of regulon activity over the cached AUCell, matching
#   each arm's covariate set. Output: results_derived_batch_<ARM>/tot_adverse_exp/<CT>/.
#
# Usage:
#   bash aceScenicT_male_arms.sh
#   SMOKE_FLAG=--smoke bash aceScenicT_male_arms.sh
#   SKIP_BUILD=1 bash aceScenicT_male_arms.sh    # only re-run per-arm association

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"
source "${REPO_ROOT}/Analysis/ACE/_shared/arm_covariates.sh"

INTEGRATION="derived_batch"
PHENOTYPE="tot_adverse_exp"
OUTPUT_ROOT="${ACE_OUTPUT_ROOT}/SCENIC/Tsai"
BUILD_ROOT="${OUTPUT_ROOT}/results_${INTEGRATION}/${PHENOTYPE}"   # canonical Male_<CT> AUCell dirs
LOG_DIR="${OUTPUT_ROOT}/logs"
mkdir -p "${LOG_DIR}"

PART="-p pi_lhtsai,pi_manoli"
SMOKE_FLAG="${SMOKE_FLAG:-}"
SKIP_BUILD="${SKIP_BUILD:-0}"

read -r -a ARMS <<< "${ARMS_OVERRIDE:-${NON_ANCOVA_ARMS[*]}}"

echo "=== ACE SCENIC: male AD-model arms (build-once / associate-per-arm) ==="
echo "Arms:       ${ARMS[*]}"
echo "Cell types: ${SIGNAL_CELLTYPES[*]}"
echo "Smoke:      ${SMOKE_FLAG:-(none)}   Skip build: ${SKIP_BUILD}"
echo ""

# ---------------------------------------------------------------------------
# Stage 1: build AUCell once per cell type (males), only where missing.
# ---------------------------------------------------------------------------
declare -A BUILD_JOB   # ct -> jobid (for afterok dependency)
for CT in "${SIGNAL_CELLTYPES[@]}"; do
  H5CT="$(ct_to_h5ad "${CT}")"
  AUC="${BUILD_ROOT}/Male_${H5CT}/auc_matrix.csv"
  if [[ "${SKIP_BUILD}" != "1" && ! -s "${AUC}" ]]; then
    echo "  [build] ${CT} (file ${H5CT}): AUCell missing -> submitting GRN/AUCell"
    jid=$(sbatch --parsable -p pi_lhtsai,pi_manoli,ou_bcs_low,mit_preemptable \
      -n 20 --mem=200G -t 24:00:00 \
      --export=ALL,LAUNCHER_SCRIPT_DIR="${SCRIPT_DIR}" \
      -o "${LOG_DIR}/%j_scenic_build_${H5CT}.out" \
      -e "${LOG_DIR}/%j_scenic_build_${H5CT}.err" \
      --job-name="scenic_build_${H5CT}" \
      "${SCRIPT_DIR}/run_scenic.sh" "${INTEGRATION}" "${PHENOTYPE}" "${H5CT}" "Male")
    BUILD_JOB[$CT]="${jid}"
    echo "    job ${jid}"
  else
    echo "  [build] ${CT}: AUCell present (or skip) -> reuse ${AUC}"
  fi
done

# ---------------------------------------------------------------------------
# Stage 2: per-arm association over cached AUCell (cheap), afterok on build.
# ---------------------------------------------------------------------------
echo ""
for ARM in "${ARMS[@]}"; do
  for CT in "${SIGNAL_CELLTYPES[@]}"; do
    H5CT="$(ct_to_h5ad "${CT}")"
    AUC="${BUILD_ROOT}/Male_${H5CT}/auc_matrix.csv"
    OUTDIR="${OUTPUT_ROOT}/results_${INTEGRATION}_${ARM}/${PHENOTYPE}/${CT}"
    DEP=""
    if [[ -n "${BUILD_JOB[$CT]:-}" ]]; then
      DEP="--dependency=afterok:${BUILD_JOB[$CT]}"
    fi
    jid=$(sbatch --parsable ${PART} -n 2 --mem=16G -t 01:00:00 ${DEP} \
      --export=ALL,LAUNCHER_SCRIPT_DIR="${SCRIPT_DIR}" \
      -o "${LOG_DIR}/%j_scenic_assoc_${ARM}_${CT}.out" \
      -e "${LOG_DIR}/%j_scenic_assoc_${ARM}_${CT}.err" \
      --job-name="scenic_assoc_${ARM}_${CT}" \
      --wrap="set -e; source ${REPO_ROOT}/config/paths.sh; export HDF5_USE_FILE_LOCKING=FALSE; \
        if [[ ! -s '${AUC}' ]]; then echo 'ERROR: AUCell missing: ${AUC}'; exit 1; fi; \
        \"\${SCENIC_ANALYSIS_ENV}/bin/python\" ${SCRIPT_DIR}/scenic_associate.py \
          --auc-matrix '${AUC}' --pheno-csv \"\${ACE_SCORES_CSV}\" \
          --arm ${ARM} --phenotype ${PHENOTYPE} --cell-type ${CT} \
          --output-dir '${OUTDIR}'")
    echo "  assoc ${ARM} / ${CT}: job ${jid} ${DEP}"
  done
done

echo ""
echo "Monitor: squeue -u \$USER | grep scenic_"
