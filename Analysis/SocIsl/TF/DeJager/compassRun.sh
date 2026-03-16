#!/bin/bash
#SBATCH -n 40                  	 # Number of cores requested
#SBATCH -t 24:00:00                # Runtime in hours
#SBATCH --mem=500G              # GB memory needed (memory PER CORE)
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-user=nkhera@college.harvard.edu
#SBATCH --mail-type=ALL

source /om2/user/mabdel03/anaconda/etc/profile.d/conda.sh

cd /om/scratch/Tue/mabdel03/SocialIsolation

conda activate /om2/user/mabdel03/conda_envs/nebulaAnalysis7


compass --data pseudobulk_female_Opc.tsv --num-processes 10 --species homo_sapiens
compass --data pseudobulk_female_Oli.tsv --num-processes 10 --species homo_sapiens
compass --data pseudobulk_female_Ast.tsv --num-processes 10 --species homo_sapiens
compass --data pseudobulk_female_Inh.tsv --num-processes 10 --species homo_sapiens
compass --data pseudobulk_female_Mic.tsv --num-processes 10 --species homo_sapiens
compass --data pseudobulk_female_Exc.tsv --num-processes 10 --species homo_sapiens
compass --data pseudobulk_male_Opc.tsv --num-processes 10 --species homo_sapiens
compass --data pseudobulk_male_Oli.tsv --num-processes 10 --species homo_sapiens
compass --data pseudobulk_male_Ast.tsv --num-processes 10 --species homo_sapiens
compass --data pseudobulk_male_Inh.tsv --num-processes 10 --species homo_sapiens
compass --data pseudobulk_male_Mic.tsv --num-processes 10 --species homo_sapiens
compass --data pseudobulk_male_Exc.tsv --num-processes 10 --species homo_sapiens
