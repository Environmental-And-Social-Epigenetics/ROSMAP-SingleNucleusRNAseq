
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
import SingleCellExperiment as sce

import numpy as np

scvi.settings.seed = 0
print("Last run with scvi-tools version:", scvi.__version__)

sc.set_figure_params(figsize=(6, 6), frameon=False)
sns.set_theme()
torch.set_float32_matmul_precision("high")
save_dir = tempfile.TemporaryDirectory()

# last_store = '/net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/Preprocessed_Counts/OutputP2'
last_storeZarr = '/net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/Preprocessed_Counts/concatObjFinal.zarr'

# library_ids=[]
# obj_list=[]
# for file in os.listdir(last_store):
#     if os.path.isfile(os.path.join(last_store, file)):
#         library_id = file.split('O')[0]  # Modify this line according to your file naming scheme
#         library_ids.append(library_id)

#load adata in

adata = ad.read_zarr(last_storeZarr)

# defining the list of patients to exclude, filtering anndata
exclude_patients = [
    "2518573", "82317494", "83034844", "94430339",
    "95919181", "52311825", "50101523", "65499271"
]

adata = adata[~adata.obs['batch'].isin(exclude_patients)].copy()

# Check to confirm removal
print(f"Remaining patients: {adata.obs['batch'].unique()}")


#
sce.pp.harmony_integrate(adata, 'batch')

# scvi.model.SCVI.setup_anndata(
#     adata,
#     batch_key="batch"
# )


# autoencoder = scvi.model.SCVI(adata)
# autoencoder.train()

# adata.obsm["X_scVI"] = autoencoder.get_latent_representation()  
# adata.X = autoencoder.get_normalized_expression() #store normalized matrix (this is the decoded matrix, should in theory be the same as the one that was input


# Clear all data in obs, var, obsm, and varm for storage
# adata.obs = adata.obs[[]]     # Clears obs DataFrame without removing the structure
# adata.var = adata.var[[]]     # Clears var DataFrame without removing the structure
# adata.obsm = {}               # Clears all embeddings in obsm
# adata.varm = {}               # Clears all embeddings in varm
# del adata.raw
from scipy.sparse import csr_matrix

adata.X = adata.X.astype('float32')  # Float32 precision
adata.X = csr_matrix(adata.X)  # Sparse matrix

for key in adata.layers.keys():
    adata.layers[key] = csr_matrix(adata.layers[key])

# Save the modified AnnData object
adata.write_h5ad('/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/clearedObj2.h5ad', compression="gzip", compression_opts=4)


#load in batch correction obj where i cleared out extra layers and used normalized scvi matrix as X
# file_path = '/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/clearedObj2.h5ad' 
# file_path = '/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/postFeatureSel3.h5ad'
# adata = ad.read_h5ad(file_path)


#feature selection
sc.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=5000)

#saving object
adata.write('/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/postFeatureSel3111624.h5ad')

adata = adata[:, adata.var.highly_variable]

#scale, pca, neighbors, leiden
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=100, n_pcs=20)
sc.tl.leiden(adata,resolution=0.03,n_iterations=20)
adata.write('/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/postPCA111624.h5ad')


import anndata as ad
import scanpy as sc
import os
import numpy as np

file_path = '/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/postPCA111624.h5ad'
# file_path = '/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/postFeatureSel3.h5ad'
# adata = ad.read_h5ad(file_path)
#read in adata

#umap and leiden

sc.pp.neighbors(adata, n_pcs=30)
sc.tl.umap(adata)

sc.tl.leiden(adata)

leiden0 = sc.tl.leiden(adata, key_added="leiden_res0_2", resolution=0.2)
leiden1 = sc.tl.leiden(adata, key_added="leiden_res0_5", resolution=0.5)
leiden2 = sc.tl.leiden(adata, key_added="leiden_res1", resolution=1.0)

umap = sc.pl.umap(adata,
	color=["leiden_res0_2", "leiden_res0_5", "leiden_res1"],
	legend_loc="on data",
	save='umapLeiden2.png')

#leiden 02 picked!

sc.tl.dendrogram(adata, groupby='leiden_res0_2')

marker_genes_dict = {
    "Excitatory":["SNAP25","SLC17A7","CAMK2A"],"Inhibitory":["GAD1","GAD2","SLC6A1"],
    "Oligo":["MBP","MOBP","PLP1"],"OPC": ["PDGFRA","VCAN","CSPG4"],"Astro":["AQP4","GFAP","ALDH1L1"],
    "Micro":["CSF1R","C3","CD74"]
}
sc.pl.dotplot(adata, marker_genes_dict, "leiden_res0_2", dendrogram=True, save='dotplotMarkerGenes02.png')


# marker_genes_dict = {
#     "Excitatory":["SNAP25","SLC17A7","CAMK2A"],"Inhibitory":["GAD1","GAD2","SLC6A1"],
#     "Oligo":["MBP","MOBP","PLP1"],"OPC": ["PDGFRA","VCAN","CSPG4"],"Astro":["AQP4","GFAP","ALDH1L1"],
#     "Micro":["CSF1R","C3","CD74"]
# }
# sc.pl.dotplot(adata, marker_genes_dict, "leiden_res0_5", dendrogram=True, save='dotplotMarkerGenes05.png')


# sc.tl.dendrogram(adata, groupby='leiden_res0_2')

# marker_genes_dict = {
#     "Excitatory":["SNAP25","SLC17A7","CAMK2A"],"Inhibitory":["GAD1","GAD2","SLC6A1"],
#     "Oligo":["MBP","MOBP","PLP1"],"OPC": ["PDGFRA","VCAN","CSPG4"],"Astro":["AQP4","GFAP","ALDH1L1"],
#     "Micro":["CSF1R","C3","CD74"]
# }
# sc.pl.dotplot(adata, marker_genes_dict, "leiden_res1", dendrogram=True, save='dotplotMarkerGenes1.png')


#dotplots and dendrograms

adata.write('/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/postClustering111624.h5ad')


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
# adata = ad.read_h5ad('/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/postClustering111624.h5ad')

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

# Visualize
sc.pl.umap(adata, color='cell_type')

#Debugging for invalid structures and formatting issues

adata.obs['cell-type'] = adata.obs['cell-type'].str.replace(r'[-/]', '_', regex=True)

adata.obsm['ora_estimate'] = adata.obsm['ora_estimate'].str.replace(r'[-/]', '_', regex=True)

print(type(adata.obs))
print(adata.obs.head())


if isinstance(adata.obsm, dict):
    raise ValueError("Unexpected dictionary structure in adata.obsm.")



for key in adata.obsm.keys():
    if isinstance(adata.obsm[key], dict):
        adata.obsm[key] = pd.DataFrame.from_dict(adata.obsm[key])
    elif not hasattr(adata.obsm[key], 'shape'):
        adata.obsm[key] = np.array(adata.obsm[key])


if 'ora_estimate' in adata.obsm:
    print("Type of ora_estimate:", type(adata.obsm['ora_estimate']))
    if isinstance(adata.obsm['ora_estimate'], dict):
        print("Error: ora_estimate is a dictionary, not compatible with AnnData.")
    elif not hasattr(adata.obsm['ora_estimate'], 'shape'):
        print("Error: ora_estimate has no shape attribute.")
    else:
        print("ora_estimate shape:", adata.obsm['ora_estimate'].shape)


if 'cell_type' in adata.obs:
    print("Type of cell_type:", type(adata.obs['cell_type']))
    if isinstance(adata.obs['cell_type'], pd.Series):
        print("cell_type is correctly stored as a Series.")
    else:
        print("Warning: cell_type may not be in the correct format:", type(adata.obs['cell_type']))

import pandas as pd

if 'ora_estimate' in adata.obsm and isinstance(adata.obsm['ora_estimate'], dict):
    print("Converting ora_estimate from dict to DataFrame.")
    adata.obsm['ora_estimate'] = pd.DataFrame.from_dict(adata.obsm['ora_estimate'])

if 'ora_pvals' in adata.obsm and isinstance(adata.obsm['ora_pvals'], dict):
    print("Converting ora_pvals from dict to DataFrame.")
    adata.obsm['ora_pvals'] = pd.DataFrame.from_dict(adata.obsm['ora_pvals'])


if 'cell_type' in adata.obs and not isinstance(adata.obs['cell_type'], pd.Series):
    print("Converting cell_type to Series.")
    adata.obs['cell_type'] = pd.Series(adata.obs['cell_type'], index=adata.obs.index)


print(adata.obsm)

adata.obsm['X_pca'] = adata.obsm['X_pca']
adata.obsm['X_umap'] = adata.obsm['X_umap']
adata.obsm['ora_estimate'] = adata.obsm['ora_estimate']
adata.obsm['ora_pvals'] = adata.obsm['ora_pvals']

for key, value in adata.obs.items():
    print(f"obs Key: {key}, Type: {type(value)}")

for key, value in adata.var.items():
    print(f"var Key: {key}, Type: {type(value)}")

for key, value in adata.uns.items():
    print(f"uns Key: {key}, Type: {type(value)}")

# keys_to_remove = [
#     "_scvi", "dendrogram_leiden_res0_2", "dendrogram_leiden_res0_5", "dendrogram_leiden_res1",
#     "hvg", "leiden", "leiden_res0_2", "leiden_res0_5", "leiden_res1", "neighbors", "pca", "umap"
# ]

# # Loop through the keys and delete them if they exist
# for key in keys_to_remove: 
#     if key in adata.uns:
#         del adata.uns[key]

# print("Specified keys have been removed from adata.uns.")


adata.write_h5ad('/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/postAnnotation111624.h5ad')