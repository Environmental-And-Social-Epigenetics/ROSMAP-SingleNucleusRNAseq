#!/bin/bash
# Merge per-cell-type GSEA summaries (gsea_summary_<CT>.csv) into a single
# gsea_summary.csv per arm. Run after the per-(arm x cell type) GSEA jobs finish.
#
# Usage: bash merge_gsea_summaries.sh   [ARM ...]   (default: all non-ANCOVA arms)

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"
source "${REPO_ROOT}/Analysis/ACE/_shared/arm_covariates.sh"

INTEGRATION="derived_batch"; PHENOTYPE="tot_adverse_exp"
OUTPUT_ROOT="${ACE_OUTPUT_ROOT}/GSEA/Tsai"

ARMS=("$@"); [[ ${#ARMS[@]} -eq 0 ]] && ARMS=("${NON_ANCOVA_ARMS[@]}")

for ARM in "${ARMS[@]}"; do
  d="${OUTPUT_ROOT}/results_${INTEGRATION}_${ARM}/${PHENOTYPE}"
  [[ -d "$d" ]] || { echo "  ${ARM}: no output dir, skip"; continue; }
  parts=( "$d"/gsea_summary_*.csv )
  if [[ ! -e "${parts[0]}" ]]; then echo "  ${ARM}: no per-CT summaries found"; continue; fi
  out="$d/gsea_summary.csv"
  # header from first part, then all data rows
  head -1 "${parts[0]}" > "$out"
  for p in "${parts[@]}"; do tail -n +2 "$p" >> "$out"; done
  n=$(( $(wc -l < "$out") - 1 ))
  cts=$(ls "$d"/gsea_summary_*.csv 2>/dev/null | sed -E 's|.*/gsea_summary_(.*)\.csv|\1|' | tr '\n' ' ')
  echo "  ${ARM}: merged ${n} rows from cell types: ${cts}-> $out"
done
