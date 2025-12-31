
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

# Load the AnnData object
file_path = '/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/clearedObj2.h5ad'
adata = ad.read_h5ad(file_path)

# Select a single patient by batch code
selected_patient = adata.obs['batch'].unique()[0]  # Replace with your specific patient code
adata = adata[adata.obs['batch'] == selected_patient].copy()

# Normalizing data
sc.pp.normalize_total(adata)  # per cell
sc.pp.log1p(adata)  # log1p

# Finding highly variable features
sc.pp.highly_variable_genes(adata, flavor='seurat')
adata = adata[:, adata.var['highly_variable']]  # subsetting on this basis

# Scaling data
sc.pp.scale(adata)

# PCA
sc.tl.pca(adata, svd_solver='arpack')

# Batch correction (if applicable for a single patient)
# Harmony is used for batch correction across multiple batches,
# so it may not be meaningful to run Harmony on a single batch.
# However, if you're planning to explore embedding:
# harmony_result = hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'batch')
# adata.obsm['X_harmony'] = harmony_result.Z_corr.T 

# Neighbors, leiden clusters
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30) #, use_rep='X_harmony'
sc.tl.leiden(adata, key_added="leiden_res0_2", resolution=0.2)
sc.tl.leiden(adata, key_added="leiden_res0_5", resolution=0.5)
sc.tl.leiden(adata, key_added="leiden_res1", resolution=1.0)

# UMAP
sc.tl.umap(adata, min_dist=0.3)

# PCA Visualization
# sc.pl.pca(adata, color='batch', save="_preHarmony_PCA.pdf")

# Harmony Visualization
# sc.pl.embedding(adata, basis='X_harmony', color='batch', save="_postHarmony_PCA.pdf")

# UMAP Visualization
sc.pl.umap(adata,
           color=["leiden_res0_2", "leiden_res0_5", "leiden_res1"],
           legend_loc="on data",
           save='umapLeidenOnePatient.png')
