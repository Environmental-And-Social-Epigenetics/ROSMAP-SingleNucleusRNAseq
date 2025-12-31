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
# Activate automatic conversion between R and pandas
pandas2ri.activate() # Load the RDS file
rds_data = r['readRDS']('Brain_Human_PFC_Markers_Mohammadi2020.rds')
# Convert to a pandas DataFrame
markers_df = pandas2ri.rpy2py(rds_data)

cell_type_names = ["Ast", "Endo", "Ex-L5", "Ex-L5/6","Ex-L5/6-CC","Ex-NRGN","Ex-L4","Ex-L2/3"]

cell_types = list(markers_df)

data = []

for cell_type_name, genes in zip(cell_type_names, cell_types):
	for gene in genes:
		data.append({'source': cell_type_name, 'target': gene, 'weight': 1.0})

markers_df = pd.DataFrame(data)

# Step 2: Load .h5ad file
adata = ad.read_h5ad('/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/postClustering103124.h5ad')

# adata.raw = adata

if 'ora_estimate' in adata.obsm:
    del adata.obsm['ora_estimate']

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

adata.obs['cell_type'] = [annotation_dict[clust] for clust in adata.obs['leiden_res0_2']]

# Visualize
sc.pl.umap(adata, color='cell_type')
# Replace slashes with underscores in the values of the DataFrame
# adata.obsm['ora_estimate'] = adata.obsm['ora_estimate'].applymap(
#     lambda x: x.replace('/', '_') if isinstance(x, str) else x
# )
# adata.obsm['ora_estimate'].columns = [
#     col.replace('/', '_') for col in adata.obsm['ora_estimate'].columns
# ]
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


# Check the type and structure of ora_estimate
if 'ora_estimate' in adata.obsm:
    print("Type of ora_estimate:", type(adata.obsm['ora_estimate']))
    if isinstance(adata.obsm['ora_estimate'], dict):
        print("Error: ora_estimate is a dictionary, not compatible with AnnData.")
    elif not hasattr(adata.obsm['ora_estimate'], 'shape'):
        print("Error: ora_estimate has no shape attribute.")
    else:
        print("ora_estimate shape:", adata.obsm['ora_estimate'].shape)


# Check if cell_type exists in obs and validate it
if 'cell_type' in adata.obs:
    print("Type of cell_type:", type(adata.obs['cell_type']))
    if isinstance(adata.obs['cell_type'], pd.Series):
        print("cell_type is correctly stored as a Series.")
    else:
        print("Warning: cell_type may not be in the correct format:", type(adata.obs['cell_type']))

import pandas as pd

# Convert ora_estimate if it's a dictionary
if 'ora_estimate' in adata.obsm and isinstance(adata.obsm['ora_estimate'], dict):
    print("Converting ora_estimate from dict to DataFrame.")
    # Assuming each key in the dictionary is a column and each value is a list or array of values
    adata.obsm['ora_estimate'] = pd.DataFrame.from_dict(adata.obsm['ora_estimate'])

if 'ora_pvals' in adata.obsm and isinstance(adata.obsm['ora_pvals'], dict):
    print("Converting ora_pvals from dict to DataFrame.")
    # Assuming each key in the dictionary is a column and each value is a list or array of values
    adata.obsm['ora_pvals'] = pd.DataFrame.from_dict(adata.obsm['ora_pvals'])


# Convert cell_type to a pandas Series if it's not already
if 'cell_type' in adata.obs and not isinstance(adata.obs['cell_type'], pd.Series):
    print("Converting cell_type to Series.")
    adata.obs['cell_type'] = pd.Series(adata.obs['cell_type'], index=adata.obs.index)


print(adata.obsm)

# Ensure no nesting or duplicate entries
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


adata.write_h5ad('/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/postAnnotation110224.h5ad')