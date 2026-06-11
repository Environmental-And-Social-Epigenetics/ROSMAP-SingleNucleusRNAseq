#!/bin/bash
# Cell-type confusion QC for ACE DEG hits.
#
# Lightweight (reads small CSVs only) -- runs fine on a login node / interactively.
# For an optional batch submission, uncomment the SBATCH block below and
# `sbatch celltype_confusion.sh`. NEVER use pi_tpoggio.
#
##SBATCH -p pi_lhtsai,pi_manoli
##SBATCH -n 2
##SBATCH --mem=8G
##SBATCH -t 00:30:00
##SBATCH -o %j_celltype_confusion.out
##SBATCH -e %j_celltype_confusion.err
#
# Usage:
#   bash celltype_confusion.sh
#   ARM=MaleNoADadj bash celltype_confusion.sh
#   PADJ=0.05 LFC=1.0 SOURCE=both bash celltype_confusion.sh

set -euo pipefail

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  SCRIPT_DIR="${SLURM_SUBMIT_DIR}"
else
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

set +u
activate_env "${NEBULA_ENV}"
set -u

export HDF5_USE_FILE_LOCKING=FALSE

ARM="${ARM:-MaleContAD}"
CONTRAST="${CONTRAST:-ACEmain}"
PHENOTYPE="${PHENOTYPE:-tot_adverse_exp}"
PADJ="${PADJ:-0.05}"
LFC="${LFC:-1.0}"
SPEC_MIN="${SPEC_MIN:-0.5}"
SOURCE="${SOURCE:-both}"
TS="$(date -u +%Y-%m-%dT%H:%M:%SZ)"

python "${SCRIPT_DIR}/celltype_confusion.py" \
  --arm "${ARM}" \
  --contrast "${CONTRAST}" \
  --phenotype "${PHENOTYPE}" \
  --padj "${PADJ}" \
  --lfc "${LFC}" \
  --spec-min "${SPEC_MIN}" \
  --source "${SOURCE}" \
  --timestamp "${TS}"
