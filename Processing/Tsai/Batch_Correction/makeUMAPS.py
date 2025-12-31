import anndata as ad
import scanpy as sc
import os
import numpy as np
import matplotlib.pyplot as plt

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

# Check if PCA has been computed, if not, compute PCA
if 'X_pca' not in adata_target.obsm:
    print("Computing PCA...")
    sc.tl.pca(adata_target)
else:
    print("PCA already computed.")

# Compute UMAP if not already computed
if 'X_umap' not in adata_target.obsm:
    print("Computing neighbors for UMAP...")
    sc.pp.neighbors(adata_target, n_pcs=40)  # Adjust `n_pcs` according to your data
    print("Computing UMAP...")
    sc.tl.umap(adata_target)
else:
    print("UMAP already computed.")

# adata_target.write("/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/updated_adata_target_with_umap.h5ad")

# Create UMAP plots for each patient
# Create UMAP plots for each patient, color-coded by leiden_res0_2 clusters
for patient_id in adata_target.obs['batch'].unique():
    # Subset the AnnData object for this specific patient
    patient_data = adata_target[adata_target.obs['batch'] == patient_id, :]
    
    # Plot UMAP, color by 'leiden_res0_2', and save as PNG
    filename = f"umap_patient_{patient_id}_leiden.png"
    print(f"Creating UMAP plot for Patient {patient_id}, color-coded by leiden_res0_2...")
    sc.pl.umap(patient_data, color='leiden_res0_2', title=f'UMAP for Patient {patient_id} (leiden clusters)', save=filename)
