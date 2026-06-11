#!/bin/bash
# One-shot resubmission of the downstream jobs that hit the QOS submit cap during
# the initial run_downstream_male_arms.sh launch. Submits only the missing pieces,
# and only if there's queue headroom. Safe to re-run: skips arms/CTs whose jobs
# are already queued or whose outputs already exist.
#
#   Missing after first launch:
#     TF:     MaleBinaryAD, MaleContAD, MaleAceByAD   (3 arm-jobs)
#     SCENIC: MaleContAD/{Oli,OPC,Inh} + MaleAceByAD/{all 6}  (9 assoc jobs)
#
# Usage: bash _resubmit_downstream_gaps.sh   (re-run until it reports "all gaps submitted")

set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../../config/paths.sh"
source "${SCRIPT_DIR}/_shared/arm_covariates.sh"

CAP=140                       # observed QOS submit cap
HEADROOM_MIN=15               # need at least this many free slots before submitting
INTEGRATION="derived_batch"; PHENOTYPE="tot_adverse_exp"
PART="-p pi_lhtsai,pi_manoli"

queued() { squeue -u "$USER" -h -o "%j" 2>/dev/null; }
qdepth() { squeue -u "$USER" -h 2>/dev/null | wc -l; }

depth=$(qdepth); free=$(( CAP - depth ))
echo "Queue depth=${depth}  free=${free} (need >=${HEADROOM_MIN})"
if (( free < HEADROOM_MIN )); then
  echo "Not enough headroom yet; try again later."
  exit 3
fi

REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
SUBMITTED=0

# ---- TF gaps (one job per missing arm) ----
TF_DIR="${SCRIPT_DIR}/TFActivity/Tsai"
TF_INPUT="${ACE_OUTPUT_ROOT}/DEG/Tsai/celltype_splits_${INTEGRATION}"
CELLTYPES_CSV="$(IFS=,; echo "${SIGNAL_CELLTYPES[*]}")"
for ARM in MaleBinaryAD MaleContAD MaleAceByAD; do
  if queued | grep -q "^tfact_${ARM}$"; then echo "TF ${ARM}: already queued"; continue; fi
  RESULTS_DIR="${ACE_OUTPUT_ROOT}/TFActivity/Tsai/results_${INTEGRATION}_${ARM}/${PHENOTYPE}"
  if [[ -s "${RESULTS_DIR}/tf_summary.csv" ]]; then echo "TF ${ARM}: output exists"; continue; fi
  mkdir -p "${RESULTS_DIR}"
  jid=$(sbatch --parsable ${PART} -n 8 --mem=96G -t 8:00:00 \
    -o "${ACE_OUTPUT_ROOT}/TFActivity/Tsai/logs/%j_tfact_${ARM}.out" \
    -e "${ACE_OUTPUT_ROOT}/TFActivity/Tsai/logs/%j_tfact_${ARM}.err" \
    --job-name="tfact_${ARM}" \
    --wrap="set -e; source ${REPO_ROOT}/config/paths.sh; export HDF5_USE_FILE_LOCKING=FALSE; \
      \"\${DECOUPLER_ENV}/bin/python\" ${TF_DIR}/tf_activity_analysis.py \
        --integration ${INTEGRATION} --phenotype ${PHENOTYPE} \
        --input-dir ${TF_INPUT} --pheno-csv \"\${ACE_SCORES_CSV}\" \
        --output-dir ${RESULTS_DIR} --run-progeny \
        --arm ${ARM} --sex-filter Male --celltypes ${CELLTYPES_CSV}" 2>/dev/null) \
    && { echo "TF ${ARM}: job ${jid}"; SUBMITTED=$((SUBMITTED+1)); } \
    || { echo "TF ${ARM}: submit failed (cap?) -- stop"; exit 4; }
done

# ---- SCENIC assoc gaps ----
SC_DIR="${SCRIPT_DIR}/SCENIC/Tsai"
BUILD_ROOT="${ACE_OUTPUT_ROOT}/SCENIC/Tsai/results_${INTEGRATION}/${PHENOTYPE}"
declare -A SC_GAPS=( [MaleContAD]="Oli OPC Inh" [MaleAceByAD]="In-PV_Basket Ast Mic Oli OPC Inh" )
for ARM in "${!SC_GAPS[@]}"; do
  for CT in ${SC_GAPS[$ARM]}; do
    if queued | grep -q "^scenic_assoc_${ARM}_${CT}$"; then echo "SCENIC ${ARM}/${CT}: queued"; continue; fi
    OUTDIR="${ACE_OUTPUT_ROOT}/SCENIC/Tsai/results_${INTEGRATION}_${ARM}/${PHENOTYPE}/${CT}"
    if [[ -s "${OUTDIR}/regression_results.csv" ]]; then echo "SCENIC ${ARM}/${CT}: output exists"; continue; fi
    H5CT="$(ct_to_h5ad "${CT}")"
    AUC="${BUILD_ROOT}/Male_${H5CT}/auc_matrix.csv"
    DEP=""
    # OPC AUCell may still be building (job from first launch); gate if missing
    if [[ ! -s "${AUC}" ]]; then
      bj=$(queued | grep -q "^scenic_build_OPC$" && squeue -u "$USER" -h -o "%i %j" 2>/dev/null | awk '$2=="scenic_build_OPC"{print $1}')
      [[ -n "${bj:-}" ]] && DEP="--dependency=afterok:${bj}"
    fi
    jid=$(sbatch --parsable ${PART} -n 2 --mem=16G -t 01:00:00 ${DEP} \
      -o "${ACE_OUTPUT_ROOT}/SCENIC/Tsai/logs/%j_scenic_assoc_${ARM}_${CT}.out" \
      -e "${ACE_OUTPUT_ROOT}/SCENIC/Tsai/logs/%j_scenic_assoc_${ARM}_${CT}.err" \
      --job-name="scenic_assoc_${ARM}_${CT}" \
      --wrap="set -e; source ${REPO_ROOT}/config/paths.sh; export HDF5_USE_FILE_LOCKING=FALSE; \
        if [[ ! -s '${AUC}' ]]; then echo 'ERROR: AUCell missing: ${AUC}'; exit 1; fi; \
        \"\${SCENIC_ANALYSIS_ENV}/bin/python\" ${SC_DIR}/scenic_associate.py \
          --auc-matrix '${AUC}' --pheno-csv \"\${ACE_SCORES_CSV}\" \
          --arm ${ARM} --phenotype ${PHENOTYPE} --cell-type ${CT} \
          --output-dir '${OUTDIR}'" 2>/dev/null) \
      && { echo "SCENIC ${ARM}/${CT}: job ${jid} ${DEP}"; SUBMITTED=$((SUBMITTED+1)); } \
      || { echo "SCENIC ${ARM}/${CT}: submit failed (cap?) -- stop"; exit 4; }
  done
done

echo ""
echo "Submitted ${SUBMITTED} gap job(s) this pass."
# exit 0 only when nothing remains to submit
remaining=0
for ARM in MaleBinaryAD MaleContAD MaleAceByAD; do queued | grep -q "^tfact_${ARM}$" || \
  [[ -s "${ACE_OUTPUT_ROOT}/TFActivity/Tsai/results_${INTEGRATION}_${ARM}/${PHENOTYPE}/tf_summary.csv" ]] || remaining=$((remaining+1)); done
for ARM in "${!SC_GAPS[@]}"; do for CT in ${SC_GAPS[$ARM]}; do
  queued | grep -q "^scenic_assoc_${ARM}_${CT}$" || \
  [[ -s "${ACE_OUTPUT_ROOT}/SCENIC/Tsai/results_${INTEGRATION}_${ARM}/${PHENOTYPE}/${CT}/regression_results.csv" ]] || remaining=$((remaining+1)); done; done
echo "Gaps still unsubmitted: ${remaining}"
[[ ${remaining} -eq 0 ]] && echo "ALL GAPS SUBMITTED" && exit 0
exit 5
