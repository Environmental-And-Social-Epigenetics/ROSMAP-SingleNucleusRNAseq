#!/bin/bash
#SBATCH -p pi_lhtsai,pi_manoli,ou_bcs_low,mit_normal
#SBATCH -n 8
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH -o functional_tot_adverse_male_%j.out
#SBATCH -e functional_tot_adverse_male_%j.err

set -euo pipefail

if [[ -n "${LAUNCHER_SCRIPT_DIR:-}" ]]; then
  SCRIPT_DIR="${LAUNCHER_SCRIPT_DIR}"
elif [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  SCRIPT_DIR="${SLURM_SUBMIT_DIR}"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

BASE_OVERLAP_ROOT="${ACE_MOUSE_OVERLAP_OUTPUT_ROOT:-${ACE_OUTPUT_ROOT}/MouseOverlap/Tsai}"
OUTPUT_ROOT="${ACE_MOUSE_FUNCTIONAL_MALE_OUTPUT_ROOT:-${BASE_OVERLAP_ROOT}/functional_tot_adverse_male}"
FIGURES_DIR="${OUTPUT_ROOT}/figures"

MOUSE_DEG="${ACE_MOUSE_FUNCTIONAL_MOUSE_DEG:-${BASE_OVERLAP_ROOT}/mouse_deg_normalized.csv}"
ORTHOLOG_TABLE="${ACE_MOUSE_FUNCTIONAL_ORTHOLOG_TABLE:-${BASE_OVERLAP_ROOT}/ortholog_mapping_used.csv}"
HUMAN_GSEA_DIR="${ACE_MOUSE_FUNCTIONAL_HUMAN_GSEA_DIR:-${ACE_OUTPUT_ROOT}/GSEA/Tsai/results_derived_batch/tot_adverse_exp/Male}"

if [[ ! -s "${MOUSE_DEG}" || ! -s "${ORTHOLOG_TABLE}" ]]; then
  echo "Base mouse overlap inputs are missing; running base overlap workflow first."
  bash "${SCRIPT_DIR}/run_mouse_overlap.sh"
fi

mkdir -p "${OUTPUT_ROOT}" "${FIGURES_DIR}"

set +u
activate_env "${GSEA_ANALYSIS_ENV}"
set -u

DEFAULT_THREADS="${SLURM_CPUS_PER_TASK:-${SLURM_NTASKS:-4}}"

CUSTOM_GMT_ARGS=()
if [[ -n "${ACE_MOUSE_FUNCTIONAL_CUSTOM_GMT:-}" ]]; then
  CUSTOM_GMT_ARGS=(--custom-gmt "${ACE_MOUSE_FUNCTIONAL_CUSTOM_GMT}")
fi

SKIP_ARGS=()
if [[ "${ACE_MOUSE_FUNCTIONAL_SKIP_EXISTING:-1}" == "1" ]]; then
  SKIP_ARGS=(--skip-existing)
fi

echo "=================================================="
echo "Mouse LNB Functional Concordance vs Male Human ACE"
echo "=================================================="
echo "  Mouse DEG:      ${MOUSE_DEG}"
echo "  Ortholog table: ${ORTHOLOG_TABLE}"
echo "  Human GSEA dir: ${HUMAN_GSEA_DIR}"
echo "  Output root:    ${OUTPUT_ROOT}"
echo "  Figures dir:    ${FIGURES_DIR}"
echo "  Conda env:      ${GSEA_ANALYSIS_ENV}"
echo "=================================================="

Rscript "${SCRIPT_DIR}/functional_tot_adverse_male.R" \
  --mouse-deg "${MOUSE_DEG}" \
  --ortholog-table "${ORTHOLOG_TABLE}" \
  --human-gsea-dir "${HUMAN_GSEA_DIR}" \
  --output-dir "${OUTPUT_ROOT}" \
  --figures-dir "${FIGURES_DIR}" \
  --permutations "${ACE_MOUSE_FUNCTIONAL_PERMUTATIONS:-1000}" \
  --threads "${ACE_MOUSE_FUNCTIONAL_THREADS:-${DEFAULT_THREADS}}" \
  --fdr-threshold "${ACE_MOUSE_FUNCTIONAL_FDR_THRESHOLD:-0.2}" \
  "${CUSTOM_GMT_ARGS[@]}" \
  "${SKIP_ARGS[@]}"

echo ""
echo "Done. Key outputs:"
echo "  ${OUTPUT_ROOT}/mouse_gsea_results.csv"
echo "  ${OUTPUT_ROOT}/human_gsea_reference.csv"
echo "  ${OUTPUT_ROOT}/human_ex_aggregate_reference.csv"
echo "  ${OUTPUT_ROOT}/functional_concordance_summary.csv"
echo "  ${OUTPUT_ROOT}/functional_concordance_pathways.csv"
echo "  ${FIGURES_DIR}"
