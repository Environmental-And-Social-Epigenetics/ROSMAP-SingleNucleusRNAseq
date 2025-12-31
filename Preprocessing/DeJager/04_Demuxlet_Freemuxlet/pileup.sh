#!/bin/bash
#SBATCH -n 24                  	 # Number of cores requested
#SBATCH -t 24:00:00           	  # Runtime in hours
#SBATCH --mem=400G       		# GB memory needed (memory PER CORE)
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this filexs

cd /net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/WGS/

source ~/.bashrc

conda activate nina_env3

conda install -c bioconda samtools

samtools mpileup -f hg38.fa 190409-B5-B/outs/possorted_genome_bam.bam > sample.pileup

sbatch freemux.sh