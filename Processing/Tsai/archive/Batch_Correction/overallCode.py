
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


# Load the AnnData object
file_path = '/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/clearedObj2.h5ad'
adata = ad.read_h5ad(file_path)

import numpy as np
import pandas as pd

# counts=adata.X

# total_elements = counts.size
# total_zeros = (counts == 0).sum().sum()  # Sum of zeros
# sparsity = total_zeros / total_elements

# depth_per_sample = counts.sum(axis=0)  # Sum across rows (axis=0)

# averageDepth = np.mean(depth_per_sample)

# print("Sparsity:", sparsity)
# print("Avg Depth per sample:", averageDepth)

# normalizing data
sc.pp.normalize_total(adata)  #per cell
sc.pp.log1p(adata) #log1p

# finding highly variable features
sc.pp.highly_variable_genes(adata, flavor='seurat')
adata = adata[:, adata.var['highly_variable']]  # subsetting on this basis

# scaling data
sc.pp.scale(adata)

# sc.pp.filter_genes(adata, min_cells=5) #PLAY AROUND WITH THIS!!

# PCA
sc.tl.pca(adata, svd_solver='arpack')

#batch correction with harmony
harmony_result = hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'batch')
adata.obsm['X_harmony'] = harmony_result.Z_corr.T 

# neighbors, leiden clusters
sc.pp.neighbors(adata, use_rep='X_harmony', n_neighbors=10, n_pcs=30)
sc.tl.leiden(adata, key_added="leiden_res0_2", resolution=0.2)
sc.tl.leiden(adata, key_added="leiden_res0_5", resolution=0.5)
sc.tl.leiden(adata, key_added="leiden_res1", resolution=1.0)

# UMAP
sc.tl.umap(adata, min_dist=0.3)

# PCA
sc.pl.pca(adata, color='batch', save="_preHarmony_PCA.pdf") 

# visualize harmony!
sc.pl.embedding(adata, basis='X_harmony', color='batch', save="_postHarmony_PCA.pdf")

# visualize UMAP with leiden clusters


umap = sc.pl.umap(adata,
  color=["leiden_res0_2", "leiden_res0_5", "leiden_res1"],
  legend_loc="on data",
  save='umapLeiden2.png')
#this is my shown output

#saving object
adata.write('/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/postFeatureSel3New.h5ad')

#NEXT PART
import pandas as pd
import anndata as ad
from scipy import sparse
from sklearn.preprocessing import normalize
import decoupler as dc
import numpy as np
import pyreadr

import scanpy as sc
from rpy2.robjects import r
from rpy2.robjects import pandas2ri


#creating and formatting markers DF

pandas2ri.activate() # Load the RDS file
rds_data = r['readRDS']('Brain_Human_PFC_Markers_Mohammadi2020.rds')
markers_df = pandas2ri.rpy2py(rds_data)
cell_type_names = ["Ast", "Endo", "Ex-L5", "Ex-L5/6","Ex-L5/6-CC","Ex-NRGN","Ex-L4","Ex-L2/3"]

cell_types = list(markers_df)

data = []

for cell_type_name, genes in zip(cell_type_names, cell_types):
	for gene in genes:
		data.append({'source': cell_type_name, 'target': gene, 'weight': 1.0})

markers_df = pd.DataFrame(data)

# load in h5ad
adata = ad.read_h5ad('/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/postFeatureSel3New.h5ad')

# adata.raw = adata

if 'ora_estimate' in adata.obsm:
    del adata.obsm['ora_estimate']

 #Running ORA/whole process there

dc.run_ora(adata,markers_df,source = "source",target = "target", use_raw=False)

acts = dc.get_acts(adata, obsm_key='ora_estimate')

# We need to remove inf and set them to the maximum value observed for pvals=0
acts_v = acts.X.ravel()
max_e = np.nanmax(acts_v[np.isfinite(acts_v)])
acts.X[~np.isfinite(acts.X)] = max_e

sc.pl.umap(acts, color=['Ast', 'leiden_res0_2'], cmap='RdBu_r')
sc.pl.violin(acts, keys=['Ast'], groupby='leiden_res0_2')

df = dc.rank_sources_groups(acts, groupby='leiden_res0_2', reference='rest', method='t-test_overestim_var')

n_ctypes = 3
ctypes_dict = df.groupby('group').head(n_ctypes).groupby('group')['names'].apply(lambda x: list(x)).to_dict()

sc.pl.matrixplot(acts, ctypes_dict, 'leiden_res0_2', dendrogram=True, standard_scale='var',
                 colorbar_title='Z-scaled scores', cmap='RdBu_r')


annotation_dict = df.groupby('group').head(1).set_index('group')['names'].to_dict()

#Cell type insertion

adata.obs['cell_type'] = [annotation_dict[clust] for clust in adata.obs['leiden_res0_2']]

# visualizing
sc.pl.umap(adata, color='cell_type')

#formatting columns

adata.obs['cell_type'] = adata.obs['cell_type'].str.replace(r'[-/]', '_', regex=True)

if 'ora_pvals' in adata.obsm:
    if isinstance(adata.obsm['ora_pvals'], pd.DataFrame):
        adata.obsm['ora_pvals'].columns = [
            col.replace('/', '_') for col in adata.obsm['ora_pvals'].columns
        ]
    elif isinstance(adata.obsm['ora_pvals'], dict):
        adata.obsm['ora_pvals'] = {
            key.replace('/', '_'): value for key, value in adata.obsm['ora_pvals'].items()
        }

adata.obsm['ora_estimate'] = adata.obsm['ora_estimate'].applymap(
    lambda x: x.replace('/', '_') if isinstance(x, str) else x
)
adata.obsm['ora_estimate'].columns = [
    col.replace('/', '_') for col in adata.obsm['ora_estimate'].columns
]

adata.write_h5ad('/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/postAnnotation111624.h5ad')

import anndata as ad
import scanpy as sc
import os
import numpy as np
import matplotlib.pyplot as plt

# source + target adata - for batch info transfer
file_path_source = '/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/batch_corrected.h5ad'
adata_source = ad.read_h5ad(file_path_source)

file_path_target = '/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/postClustering103124.h5ad'
adata_target = ad.read_h5ad(file_path_target)


adata_source.obs_names_make_unique()
adata_target.obs_names_make_unique()

# get and align batch info!!!
sample_series = adata_source.obs['batch'][adata_source.obs_names.isin(adata_target.obs_names)]
adata_target.obs['batch'] = sample_series.reindex(adata_target.obs_names)

# PCA!
if 'X_pca' not in adata_target.obsm:
    print("Computing PCA...")
    sc.tl.pca(adata_target)
else:
    print("PCA already computed.")

# UMAP!
if 'X_umap' not in adata_target.obsm:
    print("Computing neighbors for UMAP...")
    sc.pp.neighbors(adata_target, n_pcs=40)  # Adjust `n_pcs` according to your data
    print("Computing UMAP...")
    sc.tl.umap(adata_target)
else:
    print("UMAP already computed.")


# UMAP plots for each patient, color-coded by leiden_res0_2 clusters
for patient_id in adata_target.obs['batch'].unique():
    patient_data = adata_target[adata_target.obs['batch'] == patient_id, :]
    
    filename = f"umap_patient_{patient_id}_leiden.png"
    print(f"Creating UMAP plot for Patient {patient_id}, color-coded by leiden_res0_2...")
    sc.pl.umap(patient_data, color='leiden_res0_2', title=f'UMAP for Patient {patient_id} (leiden clusters)', save=filename)



