#!/bin/bash
#SBATCH -n 24                  	 # Number of cores requested
#SBATCH -t 48:00:00           	  # Runtime in hours
#SBATCH --mem=400G       		# GB memory needed (memory PER CORE)
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this filexs

# module load openmind/singularity

screen

cd /net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/WGS/

# cd /net/vast-storage/scratch/vast/lhtsai/mabdel03/Docker_Images

module load openmind/singularity/3.10.4 
# singularity build demuxafy.simg docker-archive://drneavin-demuxafy.tar
# singularity shell vibsinglecellnf-popscle.simg

# singularity exec statgen-popscle.simg dsc-pileup --sam  /net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/WGS/190409-B5-B/outs/possorted_genome_bam.bam --out /net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/WGS/output/$pileup

unset SINGULARITY_VERIFY_CHECKS

#barcodes -> 190409-B5-B/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
# singularity exec Demuxafy.sif popscle_pileup.py --sam 190409-B5-B/outs/possorted_genome_bam.bam --group-list --out output/pileup
# singularity exec Demuxafy.sif popscle_pileup.py --sam 190409-B5-B/outs/possorted_genome_bam.bam --sm-list 190409-B5-B/outs/filtered_feature_bc_matrix/barcodes.tsv.gz --out output/pileup 
# --R hg38.fa
# samtools mpileup -f hg38.fa 190409-B5-B/outs/possorted_genome_bam.bam > sample.pileup

# ls -l /net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/WGS

srun -n 32 -t 48:00:00 --mem=500G --pty bash

singularity shell --bind /net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/WGS:/mnt Demuxafy.sif

# ls -l /mnt

cd /mnt

# popscle_pileup.py --sam /mnt/190409-B5-B/outs/possorted_genome_bam.bam --vcf /mnt/variants.vcf.gz --out /mnt/plp
#/mnt/190409-B5-B/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
#/mnt/barcodes/barcodes.tsv.gz
popscle freemuxlet --plp /mnt/plp/plp --group-list /mnt/190409-B5-B/outs/filtered_feature_bc_matrix/barcodes.tsv.gz --out /mnt/output/freemux --nsample 8

# popscle freemuxlet --plp /net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/WGS/plp --group-list /net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/WGS/barcodes --out output/freemux --nsample 8
# cd /n/data2/bch/ped/raju/nkhera-renv/
# singularity build Demuxafy.sif <spec>
# --tag-group "HY3M2CCXY" --tag-UMI "190409-B5-B"

# singularity exec Demuxafy.sif popscle freemuxlet --plp /mnt/plp --group-list /mnt/barcodes --out output/freemux --nsample 8


# Rscript demuxlet.Rscript
exit



# singularity exec souporcell_latest.sif souporcell_pipeline.py -i 190409-B5-B/outs/possorted_genome_bam.bam -b 190409-B5-B/outs/filtered_feature_bc_matrix/barcodes.tsv.gz -f hg38.fa -t 5 -o output/freemux -k 8
