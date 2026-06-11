#!/bin/bash
#SBATCH -p pi_lhtsai,pi_manoli,ou_bcs_normal,mit_normal
#SBATCH -n 4
#SBATCH --mem=32G
#SBATCH -t 4:00:00
#SBATCH -o %j.out
#SBATCH -e %j.err

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

set +u
activate_env "${NEBULA_ENV}"
set -u

export HDF5_USE_FILE_LOCKING=FALSE

OUTPUT_ROOT="${ACE_MOUSE_OVERLAP_OUTPUT_ROOT:-${ACE_OUTPUT_ROOT}/MouseOverlap/Tsai}"
MOUSE_DIR="${ACE_MOUSE_RESULTS_DIR:-${REPO_ROOT}/Analysis/ACE/Mouse_Results}"
HUMAN_DEG_ROOT="${ACE_HUMAN_DEG_ROOT:-${REPO_ROOT}/Analysis/ACE/DEG/Tsai}"
PROP_ROOT="${ACE_PROP_OUTPUT_ROOT:-${ACE_OUTPUT_ROOT}/CellTypeProportion/Tsai}"

if ! find "${HUMAN_DEG_ROOT}" -maxdepth 3 -name 'deseqAnalysisACE_*.rda' -print -quit 2>/dev/null | grep -q .; then
  FALLBACK_HUMAN_DEG_ROOT="${ACE_OUTPUT_ROOT}/DEG/Tsai"
  if find "${FALLBACK_HUMAN_DEG_ROOT}" -maxdepth 3 -name 'deseqAnalysisACE_*.rda' -print -quit 2>/dev/null | grep -q .; then
    HUMAN_DEG_ROOT="${FALLBACK_HUMAN_DEG_ROOT}"
  fi
fi

mkdir -p "${OUTPUT_ROOT}" "${OUTPUT_ROOT}/figures"

echo "============================================="
echo "ACE Mouse LNB Overlap - Tsai"
echo "============================================="
echo "  Repo root:       ${REPO_ROOT}"
echo "  Mouse DEG dir:   ${MOUSE_DIR}"
echo "  Human DEG root:  ${HUMAN_DEG_ROOT}"
echo "  Output root:     ${OUTPUT_ROOT}"
echo "  Ortholog source: ${ORTHOLOG_SOURCE:-MGI default}"
echo "  Conda env:       ${NEBULA_ENV}"
echo "============================================="

Rscript "${SCRIPT_DIR}/extract_human_deg.R" \
  --results-root "${HUMAN_DEG_ROOT}" \
  --output "${OUTPUT_ROOT}/human_deg_flat.csv"

"${NEBULA_ENV}/bin/python" "${SCRIPT_DIR}/normalize_mouse_deg.py" \
  --mouse-dir "${MOUSE_DIR}" \
  --output "${OUTPUT_ROOT}/mouse_deg_normalized.csv"

ORTHOLOG_RAW="${OUTPUT_ROOT}/HOM_MouseHumanSequence.rpt"
ORTHOLOG_USED="${OUTPUT_ROOT}/ortholog_mapping_used.csv"
ORTHOLOG_META="${OUTPUT_ROOT}/ortholog_mapping_metadata.json"

"${NEBULA_ENV}/bin/python" "${SCRIPT_DIR}/prepare_orthologs.py" \
  --source "${ORTHOLOG_SOURCE:-https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt}" \
  --raw-output "${ORTHOLOG_RAW}" \
  --output "${ORTHOLOG_USED}" \
  --metadata-output "${ORTHOLOG_META}"

"${NEBULA_ENV}/bin/python" "${SCRIPT_DIR}/overlap_analysis.py" \
  --mouse-deg "${OUTPUT_ROOT}/mouse_deg_normalized.csv" \
  --human-deg "${OUTPUT_ROOT}/human_deg_flat.csv" \
  --ortholog-table "${ORTHOLOG_USED}" \
  --output-dir "${OUTPUT_ROOT}" \
  --figures-dir "${OUTPUT_ROOT}/figures"

"${NEBULA_ENV}/bin/python" "${SCRIPT_DIR}/summarize_pathology_context.py" \
  --prop-root "${PROP_ROOT}" \
  --output-dir "${OUTPUT_ROOT}" || {
    echo "WARNING: pathology context summary failed; overlap outputs are still available." >&2
  }

echo ""
echo "Done. Key outputs:"
echo "  ${OUTPUT_ROOT}/overlap_summary.csv"
echo "  ${OUTPUT_ROOT}/overlap_genes.csv"
echo "  ${OUTPUT_ROOT}/figures"
