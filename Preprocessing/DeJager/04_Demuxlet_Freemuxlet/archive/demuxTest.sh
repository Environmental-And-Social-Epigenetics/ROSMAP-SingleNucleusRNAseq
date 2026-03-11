#!/bin/bash
#SBATCH -n 80                  	 # Number of cores requested
#SBATCH -t 48:00:00           	  # Runtime in hours
#SBATCH --mem=400G       		# GB memory needed (memory PER CORE)
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-user=nkhera@college.harvard.edu
#SBATCH --mail-type=ALL

# module load openmind/singularity

# screen
source /etc/profile.d/modules.sh

cd /om/scratch/Sun/mabdel03/ROSMAP_SC/

mv filteredSNPFreq2.vcf.gz /net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/WGS/ROSMAP

mv filteredSNPFreq2.vcf.gz.tbi /net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/WGS/ROSMAP

cd /net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/

# module add openmind/singularity/3.10.4 
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
# srun -n 80 -t 48:00:00 --mem=900G --pty bash

# srun -n 50 -t 8:00:00 --mem=300G --pty bash
#
# srun -n 10 -t 0:20:00 --mem=100G --pty bash

# bcftools view -r chr1:100000-200000 /mnt/WGS/ROSMAP/concatenated_liftedROSMAP.vcf.gz -Oz -o /mnt/WGS/ROSMAP/small_test.vcf.gz

# singularity shell --bind /net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager:/mnt Demuxafy.sif
# # /,# ls -l /mnt
# /path/to/container.sif your_command

# singularity exec --pwd /net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager --bind /net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager:/mnt Demuxafy.sif popscle_pileup.py --sam /mnt/filteredMoreBAM.bam --vcf /mnt/WGS/ROSMAP/overlapping_SNPs.vcf.gz --group-list /mnt/processed_feature_bc_matrix_cell_barcodes_shuffled.csv --out /mnt/WGS/plp/plpDemux14 --sm-list /mnt/WGS/individualsTest200306-B16-A.txt 

singularity exec --pwd /net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager --bind /net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager:/mnt Demuxafy.sif popscle_pileup.py --sam /mnt/filteredMoreBAM2.bam --vcf /mnt/WGS/ROSMAP/filteredSNPFreq2.vcf.gz --group-list /mnt/processed_feature_bc_matrix_cell_barcodes.csv --out /mnt/WGS/plp/plpDemux16 --sm-list /mnt/WGS/individualsTest200306-B16-A.txt 

# singularity exec --pwd /net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager --bind /net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager:/mnt Demuxafy.sif popscle_pileup.py --sam /mnt/filteredMoreBAM2.bam --vcf /mnt/WGS/ROSMAP/filteredSNPFreq2.vcf.gz --group-list /mnt/processed_feature_bc_matrix_cell_barcodes.csv --out /mnt/WGS/plp/plpDemux15 --sm-list /mnt/WGS/individualsTest200306-B16-A.txt 

# singularity exec --pwd /net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager --bind /net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager:/mnt Demuxafy.sif popscle_pileup.py --sam /mnt/filteredBAM.bam --vcf /mnt/WGS/ROSMAP/filteredSNPFreq.vcf.gz --group-list /mnt/processed_feature_bc_matrix_cell_barcodes_shuffled2.csv --out /mnt/WGS/plp/plpDemux14 --sm-list /mnt/WGS/individualsTest200306-B16-A.txt 

# cd /mnt

# singularity exec --pwd /net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager --bind /net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager:/mnt Demuxafy.sif popscle_pileup.py --sam /mnt/filteredMoreBAM2.bam --vcf /mnt/WGS/ROSMAP/filteredSNPFreq2.vcf.gz --group-list /mnt/processed_feature_bc_matrix_cell_barcodes_shuffled.csv --out /mnt/WGS/plp/plpDemux15 --sm-list /mnt/WGS/individualsTest200306-B16-A.txt 


#combine vcf
# popscle_pileup.py --sam /mnt/filteredBAMCellSubset2.bam --vcf /mnt/WGS/ROSMAP/filteredFixedliftedROSMAPVCF200306-B16-A.vcf.gz --group-list /mnt/processed_feature_bc_matrix_cell_barcodes_shuffled2.csv --out /mnt/WGS/plp/plpDemux14 --sm-list /mnt/WGS/individualsTest200306-B16-A.txt 
# 2> error_log14.txt
# singularity exec --pwd /net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager --bind /net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager:/mnt Demuxafy.sif popscle demuxlet --plp /mnt/WGS/plp/plpDemux14 --vcf /mnt/WGS/ROSMAP/overlapping_SNPs.vcf.gz --field "PL" --group-list /mnt/processed_feature_bc_matrix_cell_barcodes_shuffled.csv --sm-list /mnt/WGS/individualsTest200306-B16-A.txt --out /mnt/WGS/output/demuxNewLibrary --alpha 0.05 --min-mac 1 --doublet-prior 0.1

singularity exec --pwd /net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager --bind /net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager:/mnt Demuxafy.sif popscle demuxlet --plp /mnt/WGS/plp/plpDemux16 --vcf /mnt/WGS/ROSMAP/filteredSNPFreq2.vcf.gz --field "PL" --group-list /mnt/processed_feature_bc_matrix_cell_barcodes.csv --sm-list /mnt/WGS/individualsTest200306-B16-A.txt --out /mnt/WGS/output/demuxNewLibraryLarger2 --alpha 0.05 --min-mac 1 --doublet-prior 0.1
# 2> error_logdemuxNewLibrary.txt 
# popscle_pileup.py --sam /mnt/191121-B6/outs/possorted_genome_bam.bam --vcf /mnt/WGS/ROSMAP/filteredFixedVCF.vcf.gz --group-list /mnt/TestLibrary_Rerun/Cellbender_Output/barcodesFiltered.csv --out /mnt/WGS/plp/plpDemux11 --sm-list /mnt/WGS/individualsTest191121-B6.txt 2> error_log11.txt


# # bcftools index yourfile.vcf.gz
# popscle_pileup.py --sam /mnt/191121-B6/outs/possorted_genome_bam.bam --vcf /mnt/WGS/variants191121-B6.vcf.gz --group-list /mnt/CellBender_Testing/Test_191121-B6/Output/processed_feature_bc_matrix_cell_barcodes.csv --out /mnt/WGS/plp/plpDemux7 --sm-list /mnt/WGS/individualsTest191121-B6.txt  2> error_log7.txt


# popscle_pileup.py --sam /mnt/191121-B6/outs/possorted_genome_bam.bam --vcf /mnt/WGS/ROSMAP/concatenated_liftedROSMAPTestLibrary191121-B6C.vcf.gz --group-list /mnt/CellBender_Testing/Test_191121-B6/Output/processed_feature_bc_matrix_cell_barcodes.csv --out /mnt/WGS/plp/plpDemux7 --sm-list /mnt/WGS/individualsTest191121-B6.txt 2> error_log7.txt
# #/mnt/190409-B5-B/outs/filtered_feature_bc_matrix/barcodes.tsv.gz
# #/mnt/barcodes/barcodes.tsv.gz
# # --tag-UMI "UB"
# #NoChr18

# #search for sample IDs in conversion chart

# popscle demuxlet --plp /mnt/WGS/plp/plpDemux10 --vcf /mnt/WGS/ROSMAP/filteredFixedVCF.vcf.gz  --field "PL" --group-list /mnt/CellBender_Testing/Test_191121-B6/Output/processed_feature_bc_matrix_cell_barcodes.csv --sm-list /mnt/WGS/individualsTest191121-B6.txt --out /mnt/WGS/output/demuxAfterChr2 --geno-error-offset 0.2 --doublet-prior 0.35  2> error_logDemuxChrAfter2.txt 
# #0.841
# popscle demuxlet --plp /mnt/WGS/plp/plpDemux10 --vcf /mnt/WGS/ROSMAP/filteredFixedVCF.vcf.gz  --field "PL" --group-list /mnt/CellBender_Testing/Test_191121-B6/Output/processed_feature_bc_matrix_cell_barcodes.csv --sm-list /mnt/WGS/individualsTest191121-B6.txt --out /mnt/WGS/output/demuxAfterChr3 --geno-error-offset 0.2 --doublet-prior 0.25  2> error_logDemuxChrAfter3.txt 
# #0.841
# popscle demuxlet --plp /mnt/WGS/plp/plpDemux10 --vcf /mnt/WGS/ROSMAP/filteredFixedVCF.vcf.gz  --field "PL" --group-list /mnt/CellBender_Testing/Test_191121-B6/Output/processed_feature_bc_matrix_cell_barcodes.csv --sm-list /mnt/WGS/individualsTest191121-B6.txt --out /mnt/WGS/output/demuxAfterChr4 --geno-error-offset 0.3 --doublet-prior 0.35  2> error_logDemuxChrAfter4.txt 
# #0.857
# popscle demuxlet --plp /mnt/WGS/plp/plpDemux10 --vcf /mnt/WGS/ROSMAP/filteredFixedVCF.vcf.gz  --field "PL" --group-list /mnt/CellBender_Testing/Test_191121-B6/Output/processed_feature_bc_matrix_cell_barcodes.csv --sm-list /mnt/WGS/individualsTest191121-B6.txt --out /mnt/WGS/output/demuxAfterChr5 --geno-error-offset 0.45 --doublet-prior 0.35  2> error_logDemuxChrAfter5.txt 
# #0.8629
# popscle demuxlet --plp /mnt/WGS/plp/plpDemux10 --vcf /mnt/WGS/ROSMAP/filteredFixedVCF.vcf.gz  --field "PL" --group-list /mnt/CellBender_Testing/Test_191121-B6/Output/processed_feature_bc_matrix_cell_barcodes.csv --sm-list /mnt/WGS/individualsTest191121-B6.txt --out /mnt/WGS/output/demuxAfterChr5 --geno-error-offset 0.4 --doublet-prior 0.35  2> error_logDemuxChrAfter5.txt 
# #0.8619
# popscle demuxlet --plp /mnt/WGS/plp/plpDemux10 --vcf /mnt/WGS/ROSMAP/filteredFixedVCF.vcf.gz  --field "PL" --group-list /mnt/CellBender_Testing/Test_191121-B6/Output/processed_feature_bc_matrix_cell_barcodes.csv --sm-list /mnt/WGS/individualsTest191121-B6.txt --out /mnt/WGS/output/demuxAfterChr6 --geno-error-offset 0.5 --doublet-prior 0.35 2> error_logDemuxChrAfter6.txt 
# #0.8627
# popscle demuxlet --plp /mnt/WGS/plp/plpDemux10 --vcf /mnt/WGS/ROSMAP/filteredFixedVCF.vcf.gz  --field "PL" --group-list /mnt/CellBender_Testing/Test_191121-B6/Output/processed_feature_bc_matrix_cell_barcodes.csv --sm-list /mnt/WGS/individualsTest191121-B6.txt --out /mnt/WGS/output/demuxAfterChr7 --geno-error-offset 0.45 --min-total 1000 --doublet-prior 0.35 2> error_logDemuxChrAfter7.txt 
# #0.8857206072230955
# popscle demuxlet --plp /mnt/WGS/plp/plpDemux10 --vcf /mnt/WGS/ROSMAP/filteredFixedVCF.vcf.gz  --field "PL" --group-list /mnt/CellBender_Testing/Test_191121-B6/Output/processed_feature_bc_matrix_cell_barcodes.csv --sm-list /mnt/WGS/individualsTest191121-B6.txt --out /mnt/WGS/output/demuxAfterChr7 --geno-error-offset 0.45 --min-total 1200 --doublet-prior 0.35 2> error_logDemuxChrAfter7.txt 
# #0.9158428910714379
# popscle demuxlet --plp /mnt/WGS/plp/plpDemux10 --vcf /mnt/WGS/ROSMAP/filteredFixedVCF.vcf.gz  --field "PL" --group-list /mnt/CellBender_Testing/Test_191121-B6/Output/processed_feature_bc_matrix_cell_barcodes.csv --sm-list /mnt/WGS/individualsTest191121-B6.txt --out /mnt/WGS/output/demuxAfterChr6 --geno-error-offset 0.5 --cap-BQ 15 --doublet-prior 0.35 2> error_logDemuxChrAfter6.txt 
# #0.8637
# popscle demuxlet --plp /mnt/WGS/plp/plpDemux10 --vcf /mnt/WGS/ROSMAP/filteredFixedVCF.vcf.gz  --field "PL" --group-list /mnt/CellBender_Testing/Test_191121-B6/Output/processed_feature_bc_matrix_cell_barcodes.csv --sm-list /mnt/WGS/individualsTest191121-B6.txt --out /mnt/WGS/output/demuxAfterChr8 --geno-error-offset 0.5 --cap-BQ 15 --min-total 1400 --doublet-prior 0.35 2> error_logDemuxChrAfter8.txt 
# #0.9453410519379127
# popscle demuxlet --plp /mnt/WGS/plp/plpDemux10 --vcf /mnt/WGS/ROSMAP/filteredFixedVCF.vcf.gz  --field "PL" --group-list /mnt/CellBender_Testing/Test_191121-B6/Output/processed_feature_bc_matrix_cell_barcodes.csv --sm-list /mnt/WGS/individualsTest191121-B6.txt --out /mnt/WGS/output/demuxAfterChr8 --geno-error-offset 0.5 --cap-BQ 15 --alpha 0.01 --doublet-prior 0.35 2> error_logDemuxChrAfter8.txt 
# #0.864365253971824
# popscle demuxlet --plp /mnt/WGS/plp/plpDemux10 --vcf /mnt/WGS/ROSMAP/filteredFixedVCF.vcf.gz  --field "PL" --group-list /mnt/CellBender_Testing/Test_191121-B6/Output/processed_feature_bc_matrix_cell_barcodes.csv --sm-list /mnt/WGS/individualsTest191121-B6.txt --out /mnt/WGS/output/demuxAfterChr8 --geno-error-offset 0.5 --cap-BQ 15 --alpha 0.8 --doublet-prior 0.35 2> error_logDemuxChrAfter8.txt 
# #0.822752672283886
# popscle demuxlet --plp /mnt/WGS/plp/plpDemux10 --vcf /mnt/WGS/ROSMAP/filteredFixedVCF.vcf.gz  --field "PL" --group-list /mnt/CellBender_Testing/Test_191121-B6/Output/processed_feature_bc_matrix_cell_barcodes.csv --sm-list /mnt/WGS/individualsTest191121-B6.txt --out /mnt/WGS/output/demuxAfterChr8 --geno-error-offset 0.5 --cap-BQ 15 --alpha 0.1 --doublet-prior 0.35 2> error_logDemuxChrAfter8.txt 
# #0.8644609881783641
# popscle demuxlet --plp /mnt/WGS/plp/plpDemux10 --vcf /mnt/WGS/ROSMAP/filteredFixedVCF.vcf.gz  --field "PL" --group-list /mnt/CellBender_Testing/Test_191121-B6/Output/processed_feature_bc_matrix_cell_barcodes.csv --sm-list /mnt/WGS/individualsTest191121-B6.txt --out /mnt/WGS/output/demuxAfterChr8 --geno-error-offset 0.5 --cap-BQ 15 --alpha 0.1 --min-mac 5 --doublet-prior 0.35 2> error_logDemuxChrAfter8.txt 
# #0.6096847126899851
# popscle demuxlet --plp /mnt/WGS/plp/plpDemux10 --vcf /mnt/WGS/ROSMAP/filteredFixedVCF.vcf.gz  --field "PL" --group-list /mnt/CellBender_Testing/Test_191121-B6/Output/processed_feature_bc_matrix_cell_barcodes.csv --sm-list /mnt/WGS/individualsTest191121-B6.txt --out /mnt/WGS/output/demuxAfterChr8 --geno-error-offset 0.5 --cap-BQ 15 --alpha 0.1 --min-callrate 0.4 --doublet-prior 0.35 2> error_logDemuxChrAfter8.txt 
# #0.8644609881783641
# popscle demuxlet --plp /mnt/WGS/plp/plpDemux10 --vcf /mnt/WGS/ROSMAP/filteredFixedVCF.vcf.gz  --field "PL" --group-list /mnt/CellBender_Testing/Test_191121-B6/Output/processed_feature_bc_matrix_cell_barcodes.csv --sm-list /mnt/WGS/individualsTest191121-B6.txt --out /mnt/WGS/output/demuxAfterChr8 --geno-error-offset 0.5 --cap-BQ 15 --alpha 0.1 --min-total 900 --doublet-prior 0.35 2> error_logDemuxChrAfter8.txt 
# #0.87439
# popscle demuxlet --plp /mnt/WGS/plp/plpDemux10 --vcf /mnt/WGS/ROSMAP/filteredFixedVCF.vcf.gz  --field "PL" --group-list /mnt/CellBender_Testing/Test_191121-B6/Output/processed_feature_bc_matrix_cell_barcodes.csv --sm-list /mnt/WGS/individualsTest191121-B6.txt --out /mnt/WGS/output/demuxAfterChr9 --geno-error-offset 0.25 --cap-BQ 15 --alpha 0.1 --min-mac 1 --min-total 900 --doublet-prior 0.35 2> error_logDemuxChrAfter9.txt 
# #0.874
# popscle demuxlet --plp /mnt/WGS/plp/plpDemux10 --vcf /mnt/WGS/ROSMAP/filteredFixedVCF.vcf.gz  --field "PL" --group-list /mnt/TestLibrary_Rerun/Cellbender_Output/barcodesFiltered.csv --sm-list /mnt/WGS/individualsTest191121-B6.txt --out /mnt/WGS/output/demuxAfterChr9 --alpha 0.05 --min-mac 1 --doublet-prior 0.1 2> error_logDemuxChrAfter9.txt 
# #0.8570418208900564
# popscle demuxlet --plp /mnt/WGS/plp/plpDemux10 --vcf /mnt/WGS/ROSMAP/filteredFixedVCF.vcf.gz  --field "PL" --group-list /mnt/TestLibrary_Rerun/Cellbender_Output/barcodesFiltered.csv --sm-list /mnt/WGS/individualsTest191121-B6.txt --out /mnt/WGS/output/demuxAfterChr9 --geno-error-offset 0.1 --alpha 0.05 --min-mac 1 --doublet-prior 0.15 2> error_logDemuxChrAfter9.txt 
# # --min-total 900
# # --geno-error-offset 0.25

# #R4078277,R4388056,R7090624,R1969233,R9101940,R6622577,R7641350,R1710143
# #77143621,50104134,60725338,43074402,36492755,72650337,50109927,84653463
# #SM-CJGNZ,MAP50104134,SM-CJK3D,SM-CTDTN,SM-CJGLZ,SM-CTDQY,SM-CJIY5,SM-CJJ1E
# popscle demuxlet --plp /mnt/WGS/plp/plpDemux --vcf /mnt/WGS/ROSMAP/concatenated_liftedROSMAPTestLibrary191121-B6.vcf.gz  --field "PL" --group-list /mnt/CellBender_Testing/Test_191121-B6/Output/processed_feature_bc_matrix_cell_barcodes.csv --sm-list /mnt/WGS/individualsTest191121-B6.txt --out /mnt/WGS/output/demuxNew  2> error_logDemuxNew.txt

# #new /mnt/TestLibrary_Rerun/Cellbender_Output/processed_feature_bc_matrix_cell_barcodes.csv

# popscle_pileup.py --sam /mnt/191121-B6/outs/possorted_genome_bam.bam --vcf /mnt/WGS/ROSMAP/small_test.vcf.gz --group-list /mnt/CellBender_Testing/Test_191121-B6/Output/processed_feature_bc_matrix_cell_barcodes.csv --out /mnt/WGS/plp/plpDemuxSmall --sm-list /mnt/WGS/individualsTest191121-B6.txt 
# popscle demuxlet --plp /mnt/WGS/plp/plpDemux4 --vcf /mnt/WGS/ROSMAP/subset_problematic_region.vcf.gz  --field "PL" --group-list /mnt/CellBender_Testing/Test_191121-B6/Output/processed_feature_bc_matrix_cell_barcodes.csv --sm-list /mnt/WGS/individualsTest191121-B6.txt --out /mnt/WGS/output/demuxRealTest
 # --nsample 8

# popscle freemuxlet --plp /net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/WGS/plp --group-list /net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/WGS/barcodes --out output/freemux --nsample 8
# cd /n/data2/bch/ped/raju/nkhera-renv/
# singularity build Demuxafy.sif <spec>
# --tag-group "HY3M2CCXY" --tag-UMI "190409-B5-B"

# singularity exec Demuxafy.sif popscle freemuxlet --plp /mnt/plp --group-list /mnt/barcodes --out output/freemux --nsample 8


# Rscript demuxlet.Rscript
exit



# singularity exec souporcell_latest.sif souporcell_pipeline.py -i 190409-B5-B/outs/possorted_genome_bam.bam -b 190409-B5-B/outs/filtered_feature_bc_matrix/barcodes.tsv.gz -f hg38.fa -t 5 -o output/freemux -k 8
