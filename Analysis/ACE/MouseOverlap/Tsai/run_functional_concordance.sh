#!/bin/bash
#SBATCH -p pi_lhtsai,pi_manoli,ou_bcs_low,mit_normal
#SBATCH -n 8
#SBATCH --mem=64G
#SBATCH -t 12:00:00
#SBATCH -o functional_concordance_%j.out
#SBATCH -e functional_concordance_%j.err

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

BASE_OUTPUT_ROOT="${ACE_MOUSE_OVERLAP_OUTPUT_ROOT:-${ACE_OUTPUT_ROOT}/MouseOverlap/Tsai}"
OUTPUT_ROOT="${ACE_MOUSE_FUNCTIONAL_OUTPUT_ROOT:-${BASE_OUTPUT_ROOT}/functional}"
FIGURES_DIR="${OUTPUT_ROOT}/figures"

MOUSE_DEG="${ACE_MOUSE_FUNCTIONAL_MOUSE_DEG:-${BASE_OUTPUT_ROOT}/mouse_deg_normalized.csv}"
HUMAN_DEG="${ACE_MOUSE_FUNCTIONAL_HUMAN_DEG:-${BASE_OUTPUT_ROOT}/human_deg_flat.csv}"
ORTHOLOG_TABLE="${ACE_MOUSE_FUNCTIONAL_ORTHOLOG_TABLE:-${BASE_OUTPUT_ROOT}/ortholog_mapping_used.csv}"

if [[ ! -s "${MOUSE_DEG}" || ! -s "${HUMAN_DEG}" || ! -s "${ORTHOLOG_TABLE}" ]]; then
  echo "Base mouse-human overlap inputs are missing; running base overlap workflow first."
  bash "${SCRIPT_DIR}/run_mouse_overlap.sh"
fi

mkdir -p "${OUTPUT_ROOT}" "${FIGURES_DIR}"
DEFAULT_FUNCTIONAL_THREADS="${SLURM_CPUS_PER_TASK:-${SLURM_NTASKS:-4}}"

set +u
activate_env "${GSEA_ANALYSIS_ENV}"
set -u

CUSTOM_GMT_ARGS=()
if [[ -n "${ACE_MOUSE_FUNCTIONAL_CUSTOM_GMT:-}" ]]; then
  CUSTOM_GMT_ARGS=(--custom-gmt "${ACE_MOUSE_FUNCTIONAL_CUSTOM_GMT}")
fi

SKIP_ARGS=()
if [[ "${ACE_MOUSE_FUNCTIONAL_SKIP_EXISTING:-0}" == "1" ]]; then
  SKIP_ARGS=(--skip-existing)
fi

echo "============================================="
echo "ACE Mouse LNB Functional Concordance - Tsai"
echo "============================================="
echo "  Mouse DEG:      ${MOUSE_DEG}"
echo "  Human DEG:      ${HUMAN_DEG}"
echo "  Ortholog table: ${ORTHOLOG_TABLE}"
echo "  Output root:    ${OUTPUT_ROOT}"
echo "  Figures dir:    ${FIGURES_DIR}"
echo "  Conda env:      ${GSEA_ANALYSIS_ENV}"
echo "============================================="

Rscript "${SCRIPT_DIR}/functional_concordance.R" \
  --mouse-deg "${MOUSE_DEG}" \
  --human-deg "${HUMAN_DEG}" \
  --ortholog-table "${ORTHOLOG_TABLE}" \
  --output-dir "${OUTPUT_ROOT}" \
  --figures-dir "${FIGURES_DIR}" \
  --databases "${ACE_MOUSE_FUNCTIONAL_DATABASES:-geneontology_Biological_Process_noRedundant,pathway_Reactome}" \
  --human-targets "${ACE_MOUSE_FUNCTIONAL_TARGETS:-Exc,Inh,In-PV_Basket,In-PV_Chandelier}" \
  --phenotypes "${ACE_MOUSE_FUNCTIONAL_PHENOTYPES:-tot_adverse_exp,ace_aggregate,early_hh_ses}" \
  --integrations "${ACE_MOUSE_FUNCTIONAL_INTEGRATIONS:-derived_batch,projid}" \
  --sexes "${ACE_MOUSE_FUNCTIONAL_SEXES:-Male,Female}" \
  --permutations "${ACE_MOUSE_FUNCTIONAL_PERMUTATIONS:-1000}" \
  --threads "${ACE_MOUSE_FUNCTIONAL_THREADS:-${DEFAULT_FUNCTIONAL_THREADS}}" \
  "${CUSTOM_GMT_ARGS[@]}" \
  "${SKIP_ARGS[@]}"

echo ""
echo "Done. Key outputs:"
echo "  ${OUTPUT_ROOT}/functional_gsea_results.csv"
echo "  ${OUTPUT_ROOT}/functional_concordance_summary.csv"
echo "  ${OUTPUT_ROOT}/functional_concordance_pathways.csv"
echo "  ${FIGURES_DIR}"
