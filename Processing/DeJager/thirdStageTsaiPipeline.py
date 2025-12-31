
import anndata as ad
import scanpy as sc
import os
import numpy as np

import os
import tempfile
import zarr

from functools import reduce
import scanpy as sc
import scanpy.external as sce

import scvi
import seaborn as sns
import torch
from rich import print
import anndata as ad
from pathlib import Path
import os
import scanpy as sc
import harmonypy as hm
import matplotlib.pyplot as plt

import numpy as np

scvi.settings.seed = 0
print("Last run with scvi-tools version:", scvi.__version__)

sc.set_figure_params(figsize=(6, 6), frameon=False)
sns.set_theme()
torch.set_float32_matmul_precision("high")
save_dir = tempfile.TemporaryDirectory()


import scanpy as sc
import anndata as ad
import harmonypy as hm
import numpy as np
import pandas as pd
import os

# Replace 'your_directory_path' with the path to the directory you want to scan
directory_path = '/net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/Preprocessed_Counts/'

# List all folders in the directory
libraries = [folder for folder in os.listdir(directory_path) if os.path.isdir(os.path.join(directory_path, folder))]

libraries.remove("qc_plots")
libraries.remove("figures")
libraries.remove("concatObjFinal.zarr")
libraries.remove("OutputP2")

adataList = []
for library in libraries:
	file_path = '/net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/Preprocessed_Counts/'+str(library)+'OutputP2.h5ad'
	if os.path.exists(file_path):
		adataL = ad.read_h5ad(file_path)
		if 'pct_counts_mt' not in adataL.obs:
			adataL.var['mt'] = adataL.var_names.str.startswith('MT-')
			adataL.var['ribo'] = adataL.var_names.str.startswith(('RPS', 'RPL'))
			adataL.var['hb'] = adataL.var_names.str.contains('HB')  # If hemoglobin genes are relevant
			sc.pp.calculate_qc_metrics(adataL, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)
		adataL.obs['batch'] = library
		adataList.append(adataL)


# Select a single patient by batch code
# selected_patient = adata.obs['batch'].unique()[0]  # Replace with your specific patient code
# adata = adata[adata.obs['batch'] == selected_patient].copy()

adata = ad.concat(adataList, merge="same", index_unique="-{object_id}")#join="outer",


print(adata.shape)

# Normalizing data
sc.pp.normalize_total(adata)  # per cell
sc.pp.log1p(adata)  # log1p

# Finding highly variable features
sc.pp.highly_variable_genes(adata, flavor='seurat')#,n_top_genes=2000
adata = adata[:, adata.var['highly_variable']]  # subsetting on this basis

# Scaling data
sc.pp.scale(adata)

# PCA
sc.tl.pca(adata, svd_solver='arpack')#n_comps=30

#batch correction
harmony_result = hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'batch')
adata.obsm['X_harmony'] = harmony_result.Z_corr.T 

# Neighbors, leiden clusters
sc.pp.neighbors(adata,use_rep='X_harmony')#n_pcs=30
sc.tl.leiden(adata, key_added="leiden_res0_2", resolution=0.2)
sc.tl.leiden(adata, key_added="leiden_res0_5", resolution=0.5)
sc.tl.leiden(adata, key_added="leiden_res1", resolution=1.0)

# UMAP
sc.tl.umap(adata)

# PCA Visualization
# sc.pl.pca(adata, color='batch', save="_preHarmony_PCA.pdf")

# Harmony Visualization
# sc.pl.embedding(adata, basis='X_harmony', color='batch', save="_postHarmony_PCA.pdf")

# UMAP Visualization
sc.pl.umap(adata,
           color=["pct_counts_mt", "log1p_total_counts", "log1p_n_genes_by_counts","pct_counts_in_top_20_genes"],
           # color=["leiden_res0_2", "leiden_res0_5", "leiden_res1"],
           legend_loc="on data",
           save='umapLeidenTotalPatient.png')
