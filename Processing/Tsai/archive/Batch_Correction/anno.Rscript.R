# Step 1: Install necessary packages
if (!requireNamespace("DecoupleR", quietly = TRUE)) {
  devtools::install_github("saezlab/decoupleR")
}
if (!requireNamespace("Seurat", quietly = TRUE)) {
  install.packages("Seurat")
}
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  install.packages("ComplexHeatmap")
}
if (!requireNamespace("hdf5r", quietly = TRUE)) {
  install.packages("hdf5r")
}

# Step 1: Install zellkonverter if necessary
if (!requireNamespace("zellkonverter", quietly = TRUE)) {
  install.packages("zellkonverter")
}

# Load required libraries
library(zellkonverter)
library(Seurat)

# Step 2: Load necessary libraries
library(DecoupleR)
library(Seurat)
library(ComplexHeatmap)
library(hdf5r)

# Assume 'seurat_obj' is your Seurat object equivalent to 'adata' in Python
# And 'markers_df' is a data frame of gene markers

library(zellkonverter)
library(Seurat)

# Step 2: Load the .h5ad file as a SingleCellExperiment object
adata <- zellkonverter::readH5AD('/om2/user/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/postPCA101524.h5ad')

# Step 3: Convert SingleCellExperiment to Seurat object
seurat_obj <- as.Seurat(adata, counts = "X")

# Step 3: Run ORA to get marker enrichment scores
ora_results <- run_ora(
  mat = seurat_obj@assays$RNA@data,
  net = markers_df,
  source = "cell_type",
  target = "gene",
  min_n = 3,
  verbose = TRUE
)

# Store ORA results in the Seurat object metadata
seurat_obj@meta.data$ora_estimate <- ora_results$estimate

# Step 4: Process ORA results, handling Inf values
max_val <- max(ora_results$estimate[is.finite(ora_results$estimate)], na.rm = TRUE)
ora_results$estimate[is.infinite(ora_results$estimate)] <- max_val

# Step 5: Rank sources for each cluster
# Identify clusters (assume 'leiden' clustering exists in the object)
ranked_clusters <- rank_sources_groups(
  ora_results,
  seurat_obj,
  groupby = "leiden",
  reference = "rest",
  method = "t_test_overestim_var"
)

# Extract top 3 cell types for each cluster
top_3_ctypes <- ranked_clusters %>%
  group_by(group) %>%
  top_n(n = 3, wt = estimate) %>%
  pull(cell_type)

# Create a dictionary for top cell type annotations per cluster
annotation_dict <- top_3_ctypes %>%
  group_by(group) %>%
  summarise(cell_type = first(cell_type)) %>%
  deframe()

# Step 6: Add cell type column to Seurat object
seurat_obj$cell_type <- sapply(seurat_obj$leiden, function(x) annotation_dict[[as.character(x)]])

# Step 7: Visualize UMAP plot with cell types
DimPlot(seurat_obj, reduction = "umap", group.by = "cell_type")

# Additional quality control visualizations
FeaturePlot(seurat_obj, features = c("nCount_RNA", "percent.mt", "scDblFinder_score", "scDblFinder_class"))

# Step 8: Matrix plot of enrichment scores
Heatmap(
  matrix = ora_results$estimate,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  name = "Z-scaled scores",
  col = colorRamp2(c(-3, 0, 3), c("blue", "white", "red")),
  show_row_names = TRUE,
  show_column_names = TRUE
)

# Assuming 'seurat_obj' is your Seurat object
# Step 2: Convert Seurat object to AnnData
adata2 <- as.SingleCellExperiment(seurat_obj)

# Step 3: Save as .h5ad
zellkonverter::writeH5AD(adata2, '/om2/user/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/postAnnotation.h5ad'))
