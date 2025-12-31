#!/bin/bash
#SBATCH -n 32                    # Number of cores requested
#SBATCH -t 47:00:00             # Runtime in minutes or
#SBATCH --output=/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/DeJager/Cellbender/Batch_Scripts_Outs/slurm-%j.out
#SBATCH --error=/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/DeJager/Cellbender/Batch_Scripts_Err/slurm-%j.err
#SBATCH --gres=gpu:a100:1
#SBATCH --mem=128G
#SBATCH --mail-user=mabdel03@mit.edu
#SBATCH --mail-type=FAIL

source /orcd/data/lhtsai/001/om2/mabdel03/miniforge3/etc/profile.d/conda.sh

conda activate /orcd/data/lhtsai/001/om2/mabdel03/conda_envs/Cellbender_env

cd /orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/DeJager/Preprocessed_Counts/190403-B4-A

cellbender remove-background --cuda --input /om/scratch/Sun/mabdel03/ROSMAP_SC/DeJager/Counts/190403-B4-A/outs/raw_feature_bc_matrix.h5 --fpr 0 --output /orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/DeJager/Preprocessed_Counts/190403-B4-A/processed_feature_bc_matrix.h5
