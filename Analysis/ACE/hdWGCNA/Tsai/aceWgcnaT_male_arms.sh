#!/bin/bash
# hdWGCNA module enrichment for the non-ANCOVA male ACE AD-model arms.
# Build-once / associate-per-arm:
#   Stage 1 (heavy): detect co-expression modules once per signal-bearing cell
#     type (males) via wgcna_analysis.R. Reuses existing module_eigengenes.csv +
#     module_assignments.csv where present. Output: canonical
#     results_derived_batch/tot_adverse_exp/Male_<CT>/.
#   Stage 2 (cheap): per-arm module-trait (covariate-adjusted) + module-DEG
#     overlap via wgcna_associate.R. Output:
#     results_derived_batch_<ARM>/tot_adverse_exp/Male_<CT>/.
#
# Usage:
#   bash aceWgcnaT_male_arms.sh
#   SKIP_BUILD=1 bash aceWgcnaT_male_arms.sh     # only per-arm association
#   SMOKE_FLAG=--smoke bash aceWgcnaT_male_arms.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"
source "${REPO_ROOT}/Analysis/ACE/_shared/arm_covariates.sh"

INTEGRATION="derived_batch"
PHENOTYPE="tot_adverse_exp"
SPLIT_DIR="${ACE_OUTPUT_ROOT}/DEG/Tsai/celltype_splits_${INTEGRATION}"
OUTPUT_ROOT="${ACE_OUTPUT_ROOT}/hdWGCNA/Tsai"
BUILD_ROOT="${OUTPUT_ROOT}/results_${INTEGRATION}/${PHENOTYPE}"
DEG_BASE="${ACE_OUTPUT_ROOT}/DEG/Tsai"
LOG_DIR="${OUTPUT_ROOT}/logs"
mkdir -p "${LOG_DIR}"

PART="-p pi_lhtsai,pi_manoli"
SMOKE_FLAG="${SMOKE_FLAG:-}"
SKIP_BUILD="${SKIP_BUILD:-0}"

read -r -a ARMS <<< "${ARMS_OVERRIDE:-${NON_ANCOVA_ARMS[*]}}"

echo "=== ACE hdWGCNA: male AD-model arms (build-once / associate-per-arm) ==="
echo "Arms:       ${ARMS[*]}"
echo "Cell types: ${SIGNAL_CELLTYPES[*]}"
echo "Smoke:      ${SMOKE_FLAG:-(none)}   Skip build: ${SKIP_BUILD}"
echo ""

# ---------------------------------------------------------------------------
# Stage 1: detect modules once per cell type (males), only where missing.
# ---------------------------------------------------------------------------
declare -A BUILD_JOB
for CT in "${SIGNAL_CELLTYPES[@]}"; do
  H5CT="$(ct_to_h5ad "${CT}")"
  ME="${BUILD_ROOT}/Male_${CT}/module_eigengenes.csv"
  if [[ "${SKIP_BUILD}" != "1" && ! -s "${ME}" ]]; then
    INPUT_H5AD="${SPLIT_DIR}/${H5CT}.h5ad"
    echo "  [build] ${CT}: modules missing -> detect (input ${H5CT}.h5ad)"
    jid=$(sbatch --parsable ${PART} -n 16 --mem=300G -t 24:00:00 \
      --export=ALL,LAUNCHER_SCRIPT_DIR="${SCRIPT_DIR}" \
      -o "${LOG_DIR}/%j_wgcna_build_${CT}.out" \
      -e "${LOG_DIR}/%j_wgcna_build_${CT}.err" \
      --job-name="wgcna_build_${CT}" \
      "${SCRIPT_DIR}/run_wgcna.sh" "${INTEGRATION}" "${PHENOTYPE}" "Male" "${CT}" "${INPUT_H5AD}" "${OUTPUT_ROOT}")
    BUILD_JOB[$CT]="${jid}"
    echo "    job ${jid}"
  else
    echo "  [build] ${CT}: modules present (or skip) -> reuse ${ME}"
  fi
done

# ---------------------------------------------------------------------------
# Stage 2: per-arm module-trait + module-DEG overlap (cheap), afterok on build.
# ---------------------------------------------------------------------------
echo ""
WGCNA_RS="${WGCNA_ENV}/bin/Rscript"
for ARM in "${ARMS[@]}"; do
  DEG_DIR="${DEG_BASE}/results_${INTEGRATION}_${ARM}_AllCellTypes/${PHENOTYPE}"
  for CT in "${SIGNAL_CELLTYPES[@]}"; do
    MODULES_DIR="${BUILD_ROOT}/Male_${CT}"
    OUTDIR="${OUTPUT_ROOT}/results_${INTEGRATION}_${ARM}/${PHENOTYPE}/Male_${CT}"
    DEP=""
    if [[ -n "${BUILD_JOB[$CT]:-}" ]]; then DEP="--dependency=afterok:${BUILD_JOB[$CT]}"; fi
    jid=$(sbatch --parsable ${PART} -n 2 --mem=16G -t 01:00:00 ${DEP} \
      --export=ALL,LAUNCHER_SCRIPT_DIR="${SCRIPT_DIR}" \
      -o "${LOG_DIR}/%j_wgcna_assoc_${ARM}_${CT}.out" \
      -e "${LOG_DIR}/%j_wgcna_assoc_${ARM}_${CT}.err" \
      --job-name="wgcna_assoc_${ARM}_${CT}" \
      --wrap="set -e; source ${REPO_ROOT}/config/paths.sh; export HDF5_USE_FILE_LOCKING=FALSE; \
        if [[ ! -s '${MODULES_DIR}/module_eigengenes.csv' ]]; then echo 'ERROR: modules missing: ${MODULES_DIR}'; exit 1; fi; \
        \"${WGCNA_RS}\" ${SCRIPT_DIR}/wgcna_associate.R \
          --modules-dir '${MODULES_DIR}' --pheno-csv \"\${ACE_SCORES_CSV}\" \
          --arm ${ARM} --phenotype ${PHENOTYPE} --cell-type ${CT} \
          --deg-results-dir '${DEG_DIR}' --output-dir '${OUTDIR}'")
    echo "  assoc ${ARM} / ${CT}: job ${jid} ${DEP}"
  done
done

echo ""
echo "Monitor: squeue -u \$USER | grep wgcna_"
