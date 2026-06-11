#!/bin/bash
#SBATCH -p pi_lhtsai,pi_manoli
#SBATCH -n 4
#SBATCH -o %j.out
#SBATCH -e %j.err

source /orcd/data/lhtsai/001/om2/mabdel03/miniforge3/etc/profile.d/conda.sh
conda activate /orcd/data/lhtsai/001/om2/mabdel03/conda_envs/nebulaAnalysis7

export HDF5_USE_FILE_LOCKING=FALSE

cd /orcd/data/lhtsai/001/mabdel03/ROSMAP_Code/Transcriptomics/Analysis/ACE/DEG/Tsai

INTEGRATION="${1:?ERROR: integration argument required (batch or projid)}"
echo "Running DESeq2 for broad_Exc, integration: $INTEGRATION"
Rscript run_deseq_broad_Exc.R --integration "$INTEGRATION"
