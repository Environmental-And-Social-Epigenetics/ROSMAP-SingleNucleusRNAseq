import os
import tempfile
import zarr

from functools import reduce
import scanpy as sc
import scvi
import seaborn as sns
import torch
from rich import print
import anndata as ad
from pathlib import Path

import numpy as np

scvi.settings.seed = 0
print("Last run with scvi-tools version:", scvi.__version__)

sc.set_figure_params(figsize=(6, 6), frameon=False)
sns.set_theme()
torch.set_float32_matmul_precision("high")
save_dir = tempfile.TemporaryDirectory()

last_store = '/net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/Preprocessed_Counts/OutputP2'
last_storeZarr = '/net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/Preprocessed_Counts/concatObj.zarr'

obj_list = []
library_ids=[]
for file in os.listdir(last_store):
    if os.path.isfile(os.path.join(last_store, file)):
        adata = sc.read(os.path.join(last_store, file))
        adata_filtered = adata[adata.obs['scDblFinder_class'] != 'doublet']
        # Extract library ID from filename or metadata
        library_id = file.split('O')[0]  # Modify this line according to your file naming scheme
        adata_filtered.obs['batch'] = library_id
        adata_filtered.var_names_make_unique()
        obj_list.append(adata_filtered)
        library_ids.append(library_id)
        if os.path.exists(last_storeZarr):
            adata_existing = ad.read_zarr(last_storeZarr)
            adata_combined = ad.concat([adata_existing, adata_filtered])
            del adata_existing
        else:
            adata_combined = adata_filtered
        # Write the combined data back to Zarr
        adata_combined.write_zarr(last_storeZarr)
        # adata_filtered.write_zarr(last_storeZarr, append=True)
        del adata
        del adata_filtered
        del adata_combined
        torch.cuda.empty_cache()

# adata = ad.concat(obj_list)


# intermediate_path='/net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/Preprocessed_Counts/IntermediateAdataObj.h5ad'

# for i in range(1, len(obj_list)):
#     # Load the next file in chunks
#     next_adata = obj_list[i]

#     # Concatenate with the current AnnData object
#     adata = adata.concatenate(next_adata, batch_key="library_id")

#     # Save intermediate result to disk
#     adata.write_h5ad(intermediate_path)
    
#     # Clear variables to free memory
#     del next_adata

# adata.write_h5ad("/net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/Preprocessed_Counts/finalLargeAdataObj.h5ad")

# adata_concat = ad.concat(obj_list, axis=0)

# adata = obj_list[0].concatenate(*obj_list[1:], batch_key="library_id")

# batch_key = 'library_id'

# label_key = 'cell_type'



#REENTER

# # adata = sc.read_zarr(zarr_store)


# scvi.model.SCVI.setup_anndata(
#     adata,
#     layer="counts",
#     batch_key="batch"
# )
# # Normalization and log transformation

# autoencoder = scvi.model.SCVI(adata)
# autoencoder.train()

# adata.obsm["X_scVI"] = autoencoder.get_latent_representation()  
# adata.obsm["X_normalized_scVI"] = autoencoder.get_normalized_expression()
# # Save the integrated data
# adata.write_h5ad('/net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/Preprocessed_Counts/Integrated_Data.h5ad')
