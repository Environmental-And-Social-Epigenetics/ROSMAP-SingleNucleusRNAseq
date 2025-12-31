"""
Script generate and submit batch scripts for counts for each library
This script will generate counts for all of the DeJager libraries
"""
import os
import pandas as pd
import subprocess
import time  # Import time module

df = pd.read_csv('/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/DeJager/FASTQs_Download/FASTQ_Download_CSVs/Synapse_FASTQ_IDs.csv')

#FASTQs location
in_root = '/om/scratch/Mon/mabdel03/FASTQs'

#Output directory to send counts to
out_root = '/om/scratch/Mon/mabdel03/Counts'

#Directory to save count batch scripts to 
parent_scripts = '/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/DeJager/Original_Counts/Counts_Batch_Scripts'

#Make folders for each library's counts outputs
for libID in set(df['LibraryID']):
	#Construct the full directory path
	dir_path = os.path.join(out_root, libID)

	# Check if the directory already exists
	if not os.path.exists(dir_path):
		os.mkdir(dir_path)
	else:
    		print(f"Directory {dir_path} already exists.")

script_counter = 0  # Initialize counter for bash scripts

#Iterate over each library
for libID in set(df['LibraryID']):
    fastqs = os.path.join(in_root, libID) #Directory w/ FASTQs
    files = os.listdir(fastqs) #list out FASTQs
    samples = [piece.split('_')[0]+'_'+piece.split('_')[1] for piece in files]
    unique_samples = list(set(samples)) #List of distinct samples

    if len(unique_samples) == 2: #if two samples, specify both
        sample = f'{unique_samples[0]},{unique_samples[1]}'
    else: #else just the one sample
        sample = unique_samples[0]

    
    filename = libID+'_count.sh' #name the file for the count batch script --> libID_count.sh
    f = open(os.path.join(parent_scripts, filename), 'x') #make the file
    out_dir = os.path.join(out_root, libID) 
    out = os.path.join(out_root, libID)


    #sbatch info + count command
    output = f"""#!/bin/bash
#SBATCH -t 47:00:00
#SBATCH -n 32
#SBATCH --mem=128G
#SBATCH --output=/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/DeJager/Original_Counts/Counts_outs/slurm-%j.out
#SBATCH --error=/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/DeJager/Original_Counts/Counts_errors/slurm-%j.err
#SBATCH --mail-user=mabdel03@mit.edu
#SBATCH --mail-type=FAIL
export PATH=/orcd/data/lhtsai/001/om2/mabdel03/apps/yard/cellranger-8.0.0:$PATH
cellranger count --create-bam true --include-introns true --nosecondary --r1-length 26 --id {libID} --transcriptome=/orcd/data/lhtsai/001/om2/mabdel03/yard/references/human/refdata-gex-GRCh38-2020-A --sample {sample} --fastqs {fastqs} --output-dir={out_dir}
"""
    f.write(output) #write to the file
    f.close()
    sbatch_command = f'sbatch {os.path.join(parent_scripts, filename)}' #submit sbatch command
    process = subprocess.Popen(sbatch_command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()

