#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=47:00:00
#SBATCH --output=/om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/03_Cellbender/Logs/Outs/slurm-%j.out
#SBATCH --error=/om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/03_Cellbender/Logs/Errs/slurm-%j.err
#SBATCH --gres=gpu:1
#SBATCH --mem=128G
#SBATCH --mail-user=mabdel03@mit.edu
#SBATCH --mail-type=FAIL

source /om2/user/mabdel03/anaconda/etc/profile.d/conda.sh

conda activate /om2/user/mabdel03/conda_envs/Cellbender_env

mkdir -p /om/scratch/Mon/mabdel03/Tsai_Data/Cellbender_Outputs/10260309
cd /om/scratch/Mon/mabdel03/Tsai_Data/Cellbender_Outputs/10260309

cellbender remove-background --cuda --input /om/scratch/Mon/mabdel03/Tsai_Data/Cellranger_Outputs/10260309/outs/raw_feature_bc_matrix.h5 --fpr 0 --output /om/scratch/Mon/mabdel03/Tsai_Data/Cellbender_Outputs/10260309/processed_feature_bc_matrix.h5

# Mark completion
echo "10260309" >> /om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/03_Cellbender/Tracking/cellbender_completed.txt
