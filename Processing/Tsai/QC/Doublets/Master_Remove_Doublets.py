import os
import subprocess


data_root = '/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/QC/Outliers'
out_root = '/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/QC/Doublets'
batch_scripts_dir = '/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Scripts/batch_scripts/QC/Doublets'
pipeline_scripts_dir = '/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Scripts/pipeline_scripts/QC/Doublets'

projids = os.listdir(data_root)

for projid in projids:

    final_dest = os.path.join(out_root, projid)

    os.mkdir(final_dest)

    batch_file = os.path.join(batch_scripts_dir, projid+'_doublet.sh')
    r_file = os.path.join(pipeline_scripts_dir, projid+'_doublet.r')


    input_data = os.path.join(data_root, projid+f'/{projid}_outlier_filtered.h5')
    output_data = os.path.join(final_dest, projid + '_doublet_filtered.h5')

    r_output = f"""library(reticulate)
# Load necessary libraries
library(zellkonverter)
library(scDblFinder)
library(SingleCellExperiment)

# Step 1: Read the AnnData object
adata <- readH5AD("{input_data}")

# Step 2: Ensure the 'counts' assay is present
if (!"counts" %in% assayNames(adata)) {{
  assay(adata, "counts") <- assay(adata, "X")
}}

# Step 3: Run scDblFinder on the SingleCellExperiment object
adata <- scDblFinder(adata)

# Step 4: Filter out doublets
# The 'scDblFinder.class' column is created by scDblFinder with 'doublet' and 'singlet' labels
# Keep only the singlets
adata_filtered <- adata[, colData(adata)$scDblFinder.class == "singlet"]

# Step 5: Save the filtered object back as an AnnData object
writeH5AD(adata_filtered, "{output_data}")

"""
    f = open(r_file, 'x') #make the file
    f.write(r_output) #write to the file
    f.close()

    batch_output = f"""#!/bin/bash
#SBATCH -t 47:00:00
#SBATCH -n 32
#SBATCH --mem=64G
#SBATCH --mail-user=mabdel03@mit.edu
#SBATCH --mail-type=FAIL

source /orcd/data/lhtsai/001/om2/mabdel03/miniforge3/etc/profile.d/conda.sh

conda activate /orcd/data/lhtsai/001/om2/mabdel03/conda_envs/single_cell_BP

Rscript {r_file}
"""
    f = open(batch_file, 'x') #make the file
    f.write(batch_output) #write to the file
    f.close()

    sbatch_command = f'sbatch {batch_file}' #submit sbatch command
    process = subprocess.Popen(sbatch_command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()



