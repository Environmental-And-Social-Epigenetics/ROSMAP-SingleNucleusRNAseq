#!/bin/bash
#SBATCH -n 1                  	 # Number of cores requested
#SBATCH -t 12:00:00           	  # Runtime in hours
#SBATCH -p short             	 # Partition (queue) to submit to
#SBATCH --mem=100G       		# GB memory needed (memory PER CORE)
#SBATCH --open-mode=append      # append adds to outfile, truncate deletes first
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file

# module load conda2/4.2.13
# conda activate /n/data2/bch/ped/raju/nkhera-renv
cd /net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/WGS/
source ~/.bashrc
conda activate ninaEnv
mamba activate ninaMambaEnv
# conda install bioconda::demuxlet
# conda install bioconda::freemuxlet
mkdir -p output

demuxlet --sam possorted_genome_bam.bam --vcf lifted_over.vcf.gz --out output/demux

# cd /n/data2/bch/ped/raju/nkhera-renv/
# singularity build Demuxafy.sif <spec>

# Rscript demuxlet.Rscript
deactivate

