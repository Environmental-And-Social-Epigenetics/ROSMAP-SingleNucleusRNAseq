library(reticulate)
# Load necessary libraries
library(zellkonverter)
library(scDblFinder)
library(SingleCellExperiment)

# Step 1: Read the AnnData object
adata <- readH5AD("/om2/user/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/QC/Outliers/14184286/14184286_outlier_filtered.h5")

# Step 2: Ensure the 'counts' assay is present
if (!"counts" %in% assayNames(adata)) {
  assay(adata, "counts") <- assay(adata, "X")
}

# Step 3: Run scDblFinder on the SingleCellExperiment object
adata <- scDblFinder(adata)

# Step 4: Filter out doublets
# The 'scDblFinder.class' column is created by scDblFinder with 'doublet' and 'singlet' labels
# Keep only the singlets
adata_filtered <- adata[, colData(adata)$scDblFinder.class == "singlet"]

# Step 5: Save the filtered object back as an AnnData object
writeH5AD(adata_filtered, "/om2/user/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/QC/Doublets/14184286/14184286_doublet_filtered.h5")

