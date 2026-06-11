#!/bin/bash
# Top-level orchestrator: downstream analysis suite for the non-ANCOVA male ACE
# AD-model arms (MaleNoADadj, MaleNiaReagan, MaleBinaryAD, MaleContAD, MaleAceByAD)
# over the signal-bearing cell types (In-PV_Basket, Ast, Mic, Oli, OPC, Inh).
#
#   (A) DEG          -- already complete (input layer; not re-run here)
#   (B) GSEA         -- GSEA/Tsai/aceGseaT_male_arms.sh        (per-arm DEG -> WebGestaltR)
#   (C) Module enr.  -- hdWGCNA/Tsai/aceWgcnaT_male_arms.sh    (build-once / associate-per-arm)
#   (D) SCENIC       -- SCENIC/Tsai/aceScenicT_male_arms.sh    (build-once / associate-per-arm)
#   (E) TF activity  -- TFActivity/Tsai/aceTfActT_male_arms.sh (per-arm OLS)
#
# Each sub-launcher manages its own SLURM submission + build->associate
# dependencies. This wrapper just fires them (optionally a subset).
#
# Usage:
#   bash run_downstream_male_arms.sh                 # all of B,C,D,E
#   ANALYSES="B E" bash run_downstream_male_arms.sh  # only GSEA + TF
#   SMOKE_FLAG=--smoke bash run_downstream_male_arms.sh
#   SKIP_BUILD=1 bash run_downstream_male_arms.sh    # reuse existing SCENIC/WGCNA builds

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # .../Analysis/ACE
export SMOKE_FLAG="${SMOKE_FLAG:-}"
export SKIP_BUILD="${SKIP_BUILD:-0}"
ANALYSES="${ANALYSES:-B C D E}"

echo "=========================================================="
echo " ACE downstream suite -- non-ANCOVA male arms"
echo " Analyses: ${ANALYSES}"
echo " Smoke:    ${SMOKE_FLAG:-(none)}   Skip build: ${SKIP_BUILD}"
echo "=========================================================="

run_one() {
  local label="$1" script="$2"
  echo ""
  echo "----- ${label} : ${script} -----"
  bash "${script}"
}

for A in ${ANALYSES}; do
  case "${A}" in
    B) run_one "GSEA"          "${SCRIPT_DIR}/GSEA/Tsai/aceGseaT_male_arms.sh" ;;
    C) run_one "ModuleEnrich"  "${SCRIPT_DIR}/hdWGCNA/Tsai/aceWgcnaT_male_arms.sh" ;;
    D) run_one "SCENIC"        "${SCRIPT_DIR}/SCENIC/Tsai/aceScenicT_male_arms.sh" ;;
    E) run_one "TFActivity"    "${SCRIPT_DIR}/TFActivity/Tsai/aceTfActT_male_arms.sh" ;;
    A) echo ""; echo "----- DEG (A): using existing per-arm results; nothing to submit -----" ;;
    *) echo "WARN: unknown analysis code '${A}' (use A B C D E)" ;;
  esac
done

echo ""
echo "All requested sub-launchers fired. Monitor: squeue -u \$USER"
