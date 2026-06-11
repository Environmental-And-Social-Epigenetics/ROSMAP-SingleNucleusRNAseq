#!/bin/bash
#SBATCH -p pi_lhtsai
#SBATCH -n 4
#SBATCH -o %j.out
#SBATCH -e %j.err

source /orcd/data/lhtsai/001/om2/mabdel03/miniforge3/etc/profile.d/conda.sh
conda activate /orcd/data/lhtsai/001/om2/mabdel03/conda_envs/nebulaAnalysis7

export HDF5_USE_FILE_LOCKING=FALSE
PYTHON=/orcd/data/lhtsai/001/om2/mabdel03/conda_envs/nebulaAnalysis7/bin/python3

cd /orcd/data/lhtsai/001/mabdel03/ROSMAP_Code/Transcriptomics/Analysis/ACE/DEG/Tsai

INTEGRATION="${1:?ERROR: integration argument required (batch or projid)}"
echo "Running pseudobulk for broad_Exc, integration: $INTEGRATION"
$PYTHON pseudobulk_broad_Exc.py --integration "$INTEGRATION"
