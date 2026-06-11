#!/bin/bash
#SBATCH -J umap_v1_v2
#SBATCH -p pi_lhtsai,pi_manoli
#SBATCH -t 02:00:00
#SBATCH -n 4
#SBATCH --mem=80G
set -euo pipefail
REPO_ROOT=/orcd/data/lhtsai/001/mabdel03/ROSMAP_Code/Transcriptomics
source "${REPO_ROOT}/config/paths.sh"
set +u
source "${CONDA_INIT_SCRIPT}"
conda activate "${BATCHCORR_ENV}"
set -u
export HDF5_USE_FILE_LOCKING=FALSE
V1=/orcd/data/lhtsai/001/mabdel03/ROSMAP_Data/Single_Nucleus/Tsai/Processing_Outputs/03_Integrated/tsai_annotated.h5ad
V2=/orcd/data/lhtsai/001/mabdel03/ROSMAP_Data/Single_Nucleus/Tsai/Processing_Outputs/03_Integrated/tsai_annotated_v2.h5ad
OUT=/orcd/data/lhtsai/001/mabdel03/ROSMAP_Data/Single_Nucleus/Tsai/Processing_Outputs/03_Integrated/stage3b/umaps_v1_vs_v2
python "${REPO_ROOT}/Processing/Tsai/Pipeline/_umap_v1_v2.py" --v1 "${V1}" --v2 "${V2}" --out-dir "${OUT}"
