#!/bin/bash
# Resubmit WGCNA cleanly after the merge/subset fix: 6 module builds (males) +
# 30 per-arm assoc jobs gated afterok on their build. Skips already-built CTs.
set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config/paths.sh"
source "${SCRIPT_DIR}/_shared/arm_covariates.sh"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
INTEGRATION="derived_batch"; PHENOTYPE="tot_adverse_exp"
WG="${ACE_OUTPUT_ROOT}/hdWGCNA/Tsai"; DEG_BASE="${ACE_OUTPUT_ROOT}/DEG/Tsai"
SPLIT="${DEG_BASE}/celltype_splits_${INTEGRATION}"; SD="${SCRIPT_DIR}/hdWGCNA/Tsai"
PART="-p pi_lhtsai,pi_manoli"
queued(){ squeue -u "$USER" -h -o "%j" 2>/dev/null; }

declare -A BJ
echo "=== WGCNA builds ==="
for CT in "${SIGNAL_CELLTYPES[@]}"; do
  ME="${WG}/results_${INTEGRATION}/${PHENOTYPE}/Male_${CT}/module_eigengenes.csv"
  [[ -s "$ME" ]] && { echo "  ${CT}: built, skip"; continue; }
  queued | grep -q "^wgcna_build_${CT}$" && { BJ[$CT]=$(squeue -u "$USER" -h -o "%i %j" | awk -v n=wgcna_build_${CT} '$2==n{print $1;exit}'); echo "  ${CT}: already queued (${BJ[$CT]})"; continue; }
  H5CT="$(ct_to_h5ad "$CT")"
  jid=$(sbatch --parsable ${PART} -n 16 --mem=300G -t 24:00:00 \
    --export=ALL,LAUNCHER_SCRIPT_DIR="$SD" \
    -o "${WG}/logs/%j_wgcna_build_${CT}.out" -e "${WG}/logs/%j_wgcna_build_${CT}.err" \
    --job-name="wgcna_build_${CT}" \
    "${SD}/run_wgcna.sh" "${INTEGRATION}" "${PHENOTYPE}" "Male" "${CT}" "${SPLIT}/${H5CT}.h5ad" "${WG}")
  BJ[$CT]="$jid"; echo "  ${CT}: build job ${jid}"
done

echo "=== WGCNA assoc (afterok build) ==="
for ARM in "${NON_ANCOVA_ARMS[@]}"; do
  DEG_DIR="${DEG_BASE}/results_${INTEGRATION}_${ARM}_AllCellTypes/${PHENOTYPE}"
  for CT in "${SIGNAL_CELLTYPES[@]}"; do
    queued | grep -q "^wgcna_assoc_${ARM}_${CT}$" && continue
    out="${WG}/results_${INTEGRATION}_${ARM}/${PHENOTYPE}/Male_${CT}/module_trait_correlations.csv"
    [[ -s "$out" ]] && continue
    MOD="${WG}/results_${INTEGRATION}/${PHENOTYPE}/Male_${CT}"
    DEP=""; [[ -n "${BJ[$CT]:-}" ]] && DEP="--dependency=afterok:${BJ[$CT]}"
    [[ -z "$DEP" && ! -s "${MOD}/module_eigengenes.csv" ]] && { echo "  ${ARM}/${CT}: no build, skip"; continue; }
    jid=$(sbatch --parsable ${PART} -n 2 --mem=16G -t 01:00:00 ${DEP} \
      -o "${WG}/logs/%j_wgcna_assoc_${ARM}_${CT}.out" -e "${WG}/logs/%j_wgcna_assoc_${ARM}_${CT}.err" \
      --job-name="wgcna_assoc_${ARM}_${CT}" \
      --wrap="set -e; source ${REPO_ROOT}/config/paths.sh; export HDF5_USE_FILE_LOCKING=FALSE; \
        \"\${WGCNA_ENV}/bin/Rscript\" ${SD}/wgcna_associate.R \
          --modules-dir '${MOD}' --pheno-csv \"\${ACE_SCORES_CSV}\" --arm ${ARM} \
          --phenotype ${PHENOTYPE} --cell-type ${CT} --deg-results-dir '${DEG_DIR}' \
          --output-dir '${WG}/results_${INTEGRATION}_${ARM}/${PHENOTYPE}/Male_${CT}'")
    echo "  ${ARM}/${CT}: assoc ${jid} ${DEP}"
  done
done
echo "Done."
