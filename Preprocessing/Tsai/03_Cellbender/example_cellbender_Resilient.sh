#!/bin/bash
#SBATCH -n 32                    # Number of cores requested
#SBATCH -t 47:00:00             # Runtime in minutes or
#SBATCH --output=/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/Tsai/Preprocessing/Cellbender/Resilient/Batch_Scripts_Outs/slurm-%j.out
#SBATCH --error=/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/Tsai/Preprocessing/Cellbender/Resilient/Batch_Scripts_Err/slurm-%j.err
#SBATCH --gres=gpu:a100:1
#SBATCH --mem=500G
#SBATCH --mail-user=mabdel03@mit.edu
#SBATCH --mail-type=FAIL

source /orcd/data/lhtsai/001/om2/mabdel03/miniforge3/etc/profile.d/conda.sh

conda activate /orcd/data/lhtsai/001/om2/mabdel03/conda_envs/Cellbender_env

cd /orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/Tsai/Preprocessing/Preprocessed_Counts/Resilient/10101589

cellbender remove-background --cuda --input /om/scratch/Mon/mabdel03/Tsai/Resilient/Counts/10101589/outs/raw_feature_bc_matrix.h5 --fpr 0 --output /orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/Tsai/Preprocessing/Preprocessed_Counts/Resilient/10101589/processed_feature_bc_matrix.h5
