"""
Script to aggregate all Tsai ACE QC outputs into one object
Then to batch correct this aggregate object with scvi
"""

import os
import scanpy as sc
import anndata as ad 
import scvi
import seaborn as sns
import torch



in_root = '/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/QC/Doublets'

out_root = '/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction'

projids = os.listdir('/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/QC/Doublets')

obj_list = []
projid_list = []

#For all patients w/ an scDblFinder output, load output in
#load objects into a list of objects, also track projids
for projid in projids:
    try:
        file_path = os.path.join(in_root, projid, projid + '_doublet_filtered.h5')
        obj_list.append(sc.read(file_path))
        projid_list.append(projid)
    except FileNotFoundError:
        #if Doublet removal subdirectory is empty, print the projid
        print(f'{projid} does not have an scDblFinder output available!')
        
for obj, projid in zip(obj_list, projid_list):
    #batch is the projid; add batch tags for all objects
	obj.obs['batch'] = [projid for i in range(len(obj.obs))] 

#concatenate objects into one big anndata object
adata = ad.concat(obj_list)

del obj_list

#setup scvi model with layer as raw counts and batch as "batch" obs slot (containing projids)
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    batch_key="batch"
)

#train vae scvi model
vae = scvi.model.SCVI(adata)
vae.train()

#save model outputs in .obsm slot
adata.obsm["X_scVI"] = vae.get_latent_representation() #latent representation, e.g the batch corrected output
adata.X = vae.get_normalized_expression() #store normalized matrix (this is the decoded matrix, should in theory be the same as the one that was input

adata.write_h5ad('/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/batch_corrected.h5ad')





