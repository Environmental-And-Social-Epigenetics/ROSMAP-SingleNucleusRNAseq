import anndata as ad
import scanpy as sc
import os
import numpy as np
import matplotlib.pyplot as plt

# file_path = '/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/postPCA101524.h5ad'
# # file_path = '/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/postFeatureSel3.h5ad'
# adata = ad.read_h5ad(file_path)

# Load source and target AnnData objects
file_path_source = '/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/batch_corrected.h5ad'
adata_source = ad.read_h5ad(file_path_source)

file_path_target = '/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/postClustering103124.h5ad'
adata_target = ad.read_h5ad(file_path_target)

# Make obs_names unique
adata_source.obs_names_make_unique()
adata_target.obs_names_make_unique()

# Extract and align batch information
sample_series = adata_source.obs['batch'][adata_source.obs_names.isin(adata_target.obs_names)]
adata_target.obs['batch'] = sample_series.reindex(adata_target.obs_names)

cluster_patient_counts = adata_target.obs.groupby(['leiden_res0_2', 'batch']).size().unstack(fill_value=0)


# Set the index name explicitly, if needed
cluster_patient_counts.index.name = 'Cluster'

# Display the DataFrame to check that clusters are the row names
print(cluster_patient_counts)
# Save the DataFrame to a CSV file
cluster_patient_counts.to_csv("cluster_patient_counts.csv")


# Display the result
print(cluster_patient_counts)

import pandas as pd

# Load the CSV file into a DataFrame
cluster_patient_counts = pd.read_csv("cluster_patient_counts.csv", index_col=0)

# Display the original DataFrame
print("Original DataFrame:")
print(cluster_patient_counts)

# Calculate the sum of each row (axis=1 means rows)
row_sums = cluster_patient_counts.sum(axis=1)

# Create a new DataFrame where each cell is divided by its respective row sum
normalized_df = cluster_patient_counts.div(row_sums, axis=0)

# Display the new normalized DataFrame
print("Normalized DataFrame:")
print(normalized_df)

# Optionally, save the new DataFrame to a CSV file with the clusters as row names
normalized_df.to_csv("normalized_cluster_patient_counts.csv", index=True)


import pandas as pd
import numpy as np
from itertools import combinations

import seaborn as sns
import matplotlib.pyplot as plt

# Load normalized data
normalized_df = pd.read_csv("normalized_cluster_patient_counts.csv", index_col=0)

# Plot a heatmap to visualize contributions of each patient across clusters
plt.figure(figsize=(20, 12))  # Increase figure size for readability
sns.heatmap(
    normalized_df, 
    cmap="YlGnBu", 
    annot=False,  # Remove annotations to reduce clutter
    fmt=".2f", 
    cbar_kws={'label': 'Patient Contribution'}
)
plt.xlabel("Patient")
plt.ylabel("Cluster")
plt.xticks(rotation=90)  # Rotate x-axis labels for better visibility
plt.title("Patient Contribution per Cluster")

# Save the figure as a PNG file with a higher DPI for quality
plt.savefig("patient_contribution_heatmap_larger.png", format="png", dpi=300)

plt.show()


# # # Iterate over all possible combinations of three columns
# # for combo in combinations(normalized_df.columns, 3):
# #     # Sum across each row for these three columns and check if the sum exceeds 0.95
# #     condition_met = (normalized_df[list(combo)].sum(axis=1) > 0.95)
# #     # Update the final mask where any combination meets the condition
# #     final_mask |= condition_met

# # # Filter the DataFrame to include only rows that meet the condition
# # rows_with_high_contribution = normalized_df[final_mask]

# # Display the resulting DataFrame

# import scanpy as sc
# import pandas as pd
# import numpy as np
# from itertools import combinations
# from sklearn.cluster import DBSCAN

# # Compute Leiden clusters
# sc.tl.leiden(adata, resolution=1, key_added="leiden_res0_2")

# # Define marker genes dictionary
# marker_genes_dict = {
#     "Excitatory": ["SNAP25", "SLC17A7", "CAMK2A"],
#     "Inhibitory": ["GAD1", "GAD2", "SLC6A1"],
#     "Oligo": ["MBP", "MOBP", "PLP1"],
#     "OPC": ["PDGFRA", "VCAN", "CSPG4"],
#     "Astro": ["AQP4", "GFAP", "ALDH1L1"],
#     "Micro": ["CSF1R", "C3", "CD74"]
# }

# # Calculate average expression of marker genes per cluster
# cluster_avg_exp = adata.to_df().groupby(adata.obs['leiden_res0_2']).mean()

# # Assign a cell type to each cluster based on marker gene expression
# cell_type_assignments = {}
# for cluster, genes in cluster_avg_exp.iterrows():
#     highest_score = 0
#     assigned_type = None
#     for cell_type, markers in marker_genes_dict.items():
#         score = genes[markers].mean()  # Average expression of marker genes
#         if score > highest_score:
#             highest_score = score
#             assigned_type = cell_type
#     cell_type_assignments[cluster] = assigned_type

# # Map cluster cell types back to adata.obs and normalized_df
# adata.obs['cell_type'] = adata.obs['leiden_res0_2'].map(cell_type_assignments)
# normalized_df['cell_type'] = adata.obs['cell_type'].values  # Assuming order matches

# # Add UMAP coordinates to the DataFrame for spatial analysis
# umap_coords = adata.obsm['X_umap']
# normalized_df['umap_x'] = umap_coords[:, 0]
# normalized_df['umap_y'] = umap_coords[:, 1]

# # Use DBSCAN to define spatial clusters in UMAP space
# db = DBSCAN(eps=0.5, min_samples=5).fit(umap_coords)  # Adjust `eps` for tighter/looser clusters
# normalized_df['umap_cluster'] = db.labels_

# # Identify rows with high patient contribution using the previous final_mask approach
# final_mask = np.zeros(normalized_df.shape[0], dtype=bool)
# for combo in combinations(normalized_df.columns[:-4], 3):  # Exclude cell_type, umap_x, umap_y, umap_cluster columns
#     condition_met = (normalized_df[list(combo)].sum(axis=1) > 0.95)
#     final_mask |= condition_met

# # Filter and merge rows by cell type and UMAP cluster for high-contribution rows
# high_contribution_df = normalized_df[final_mask]
# merged_df = high_contribution_df.groupby(['cell_type', 'umap_cluster']).sum()

# # Display the merged DataFrame
# print("Merged DataFrame by cell type and UMAP cluster with high patient contribution:")
# print(merged_df)

