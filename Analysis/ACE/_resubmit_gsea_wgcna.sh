#!/bin/bash
# Resubmit the GSEA (per-CT fan-out, fixed) and WGCNA (reader='R' fixed) jobs once
# the SLURM queue has headroom. Idempotent: skips arms/CTs already queued or done.
# Exits 0 only when nothing remains to submit ("ALL GSEA+WGCNA SUBMITTED").
set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config/paths.sh"
source "${SCRIPT_DIR}/_shared/arm_covariates.sh"

CAP=140; NEED=10            # require >=NEED free slots before submitting a batch
depth=$(squeue -u "$USER" -h 2>/dev/null | wc -l); free=$(( CAP - depth ))
echo "Queue depth=${depth} free=${free}"
if (( free < NEED )); then echo "No headroom (need >=${NEED}); waiting."; exit 3; fi

INTEGRATION="derived_batch"; PHENOTYPE="tot_adverse_exp"
GS="${ACE_OUTPUT_ROOT}/GSEA/Tsai"; WG="${ACE_OUTPUT_ROOT}/hdWGCNA/Tsai"
DEG_BASE="${ACE_OUTPUT_ROOT}/DEG/Tsai"
SPLIT_DIR="${DEG_BASE}/celltype_splits_${INTEGRATION}"
PART="-p pi_lhtsai,pi_manoli"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
queued(){ squeue -u "$USER" -h -o "%j" 2>/dev/null; }
SUB=0

# ---- GSEA: per (arm x cell type); skip if ranked_genes already written or queued ----
for ARM in "${NON_ANCOVA_ARMS[@]}"; do
  DEG_ROOT="${DEG_BASE}/results_${INTEGRATION}_${ARM}_AllCellTypes"
  [[ -d "${DEG_ROOT}/${PHENOTYPE}" ]] || continue
  for CT in "${SIGNAL_CELLTYPES[@]}"; do
    (( free - SUB < NEED )) && { echo "headroom exhausted mid-GSEA"; break 2; }
    queued | grep -q "^gsea_${ARM}_${CT}$" && { echo "GSEA ${ARM}/${CT}: queued"; continue; }
    # done = summary for this CT exists (per-CT summary) OR ranked + all 8 rds
    sm="${GS}/results_${INTEGRATION}_${ARM}/${PHENOTYPE}/gsea_summary_${CT}.csv"
    [[ -s "$sm" ]] && { echo "GSEA ${ARM}/${CT}: done"; continue; }
    jid=$(sbatch --parsable ${PART} -n 8 --mem=64G -t 06:00:00 \
      --export=ALL,LAUNCHER_SCRIPT_DIR="${SCRIPT_DIR}/GSEA/Tsai",GSEA_CELLTYPES="${CT}",ACE_GSEA_SKIP_EXISTING=1 \
      -o "${GS}/logs/%j_gsea_${ARM}_${CT}.out" -e "${GS}/logs/%j_gsea_${ARM}_${CT}.err" \
      --job-name="gsea_${ARM}_${CT}" \
      "${SCRIPT_DIR}/GSEA/Tsai/run_enrichment.sh" "${INTEGRATION}" "${PHENOTYPE}" "Male" "${DEG_ROOT}" "${GS}" "${ARM}" 2>/dev/null) \
      && { echo "GSEA ${ARM}/${CT}: job ${jid}"; SUB=$((SUB+1)); } || { echo "GSEA ${ARM}/${CT}: submit failed (cap)"; break 2; }
  done
done

# ---- WGCNA: 6 module builds (reader=R fixed); skip if eigengenes exist or queued ----
declare -A WB
for CT in "${SIGNAL_CELLTYPES[@]}"; do
  (( free - SUB < NEED )) && { echo "headroom exhausted before WGCNA builds"; break; }
  ME="${WG}/results_${INTEGRATION}/${PHENOTYPE}/Male_${CT}/module_eigengenes.csv"
  [[ -s "$ME" ]] && { echo "WGCNA build ${CT}: done"; continue; }
  queued | grep -q "^wgcna_build_${CT}$" && { echo "WGCNA build ${CT}: queued"; continue; }
  H5CT="$(ct_to_h5ad "${CT}")"
  jid=$(sbatch --parsable ${PART} -n 16 --mem=300G -t 24:00:00 \
    --export=ALL,LAUNCHER_SCRIPT_DIR="${SCRIPT_DIR}/hdWGCNA/Tsai" \
    -o "${WG}/logs/%j_wgcna_build_${CT}.out" -e "${WG}/logs/%j_wgcna_build_${CT}.err" \
    --job-name="wgcna_build_${CT}" \
    "${SCRIPT_DIR}/hdWGCNA/Tsai/run_wgcna.sh" "${INTEGRATION}" "${PHENOTYPE}" "Male" "${CT}" "${SPLIT_DIR}/${H5CT}.h5ad" "${WG}" 2>/dev/null) \
    && { echo "WGCNA build ${CT}: job ${jid}"; WB[$CT]="$jid"; SUB=$((SUB+1)); } || { echo "WGCNA build ${CT}: submit failed (cap)"; break; }
done

# ---- WGCNA assoc per arm x CT (afterok build if just submitted, else if modules exist) ----
for ARM in "${NON_ANCOVA_ARMS[@]}"; do
  DEG_DIR="${DEG_BASE}/results_${INTEGRATION}_${ARM}_AllCellTypes/${PHENOTYPE}"
  for CT in "${SIGNAL_CELLTYPES[@]}"; do
    (( free - SUB < NEED )) && { echo "headroom exhausted in WGCNA assoc"; break 2; }
    queued | grep -q "^wgcna_assoc_${ARM}_${CT}$" && continue
    out="${WG}/results_${INTEGRATION}_${ARM}/${PHENOTYPE}/Male_${CT}/module_trait_correlations.csv"
    [[ -s "$out" ]] && continue
    MOD="${WG}/results_${INTEGRATION}/${PHENOTYPE}/Male_${CT}"
    # Gate on the build: prefer a just-submitted build (WB), else look up the
    # already-queued build job id, else (build done) no dependency needed.
    DEP=""
    if [[ -n "${WB[$CT]:-}" ]]; then
      DEP="--dependency=afterok:${WB[$CT]}"
    else
      qbuild=$(squeue -u "$USER" -h -o "%i %j" 2>/dev/null | awk -v n="wgcna_build_${CT}" '$2==n{print $1; exit}')
      [[ -n "$qbuild" ]] && DEP="--dependency=afterok:${qbuild}"
    fi
    # if no build (queued or just-submitted) and modules don't exist yet, skip
    [[ -z "$DEP" && ! -s "${MOD}/module_eigengenes.csv" ]] && continue
    jid=$(sbatch --parsable ${PART} -n 2 --mem=16G -t 01:00:00 ${DEP} \
      -o "${WG}/logs/%j_wgcna_assoc_${ARM}_${CT}.out" -e "${WG}/logs/%j_wgcna_assoc_${ARM}_${CT}.err" \
      --job-name="wgcna_assoc_${ARM}_${CT}" \
      --wrap="set -e; source ${REPO_ROOT}/config/paths.sh; export HDF5_USE_FILE_LOCKING=FALSE; \
        \"\${WGCNA_ENV}/bin/Rscript\" ${SCRIPT_DIR}/hdWGCNA/Tsai/wgcna_associate.R \
          --modules-dir '${MOD}' --pheno-csv \"\${ACE_SCORES_CSV}\" --arm ${ARM} \
          --phenotype ${PHENOTYPE} --cell-type ${CT} --deg-results-dir '${DEG_DIR}' \
          --output-dir '${WG}/results_${INTEGRATION}_${ARM}/${PHENOTYPE}/Male_${CT}'" 2>/dev/null) \
      && { echo "WGCNA assoc ${ARM}/${CT}: job ${jid} ${DEP}"; SUB=$((SUB+1)); } || { echo "WGCNA assoc ${ARM}/${CT}: cap"; break 2; }
  done
done

echo "Submitted ${SUB} job(s) this pass."
# remaining check
rem=0
for ARM in "${NON_ANCOVA_ARMS[@]}"; do for CT in "${SIGNAL_CELLTYPES[@]}"; do
  [[ -s "${GS}/results_${INTEGRATION}_${ARM}/${PHENOTYPE}/gsea_summary_${CT}.csv" ]] || queued | grep -q "^gsea_${ARM}_${CT}$" || rem=$((rem+1))
done; done
for CT in "${SIGNAL_CELLTYPES[@]}"; do
  [[ -s "${WG}/results_${INTEGRATION}/${PHENOTYPE}/Male_${CT}/module_eigengenes.csv" ]] || queued | grep -q "^wgcna_build_${CT}$" || rem=$((rem+1))
done
# WGCNA assoc jobs (done = module_trait_correlations.csv exists, or queued)
for ARM in "${NON_ANCOVA_ARMS[@]}"; do for CT in "${SIGNAL_CELLTYPES[@]}"; do
  [[ -s "${WG}/results_${INTEGRATION}_${ARM}/${PHENOTYPE}/Male_${CT}/module_trait_correlations.csv" ]] || \
    queued | grep -q "^wgcna_assoc_${ARM}_${CT}$" || rem=$((rem+1))
done; done
echo "Remaining unsubmitted (GSEA + WGCNA build/assoc): ${rem}"
[[ ${rem} -eq 0 ]] && echo "ALL GSEA+WGCNA SUBMITTED" && exit 0
exit 5
