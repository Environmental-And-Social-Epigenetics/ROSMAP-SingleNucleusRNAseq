
import anndata as ad
import scanpy as sc
import os
import anndata2ri
import numpy as np

file_path = '/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/postPCA101524.h5ad'
# file_path = '/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/postFeatureSel3.h5ad'
adata = ad.read_h5ad(file_path)

import scanpy as sc

#NEXT STEPS

# umap = sc.pl.umap(adata, color="total_counts", save='umapAdata.png')

import scanpy as sc
import os

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

sc.tl.dendrogram(adata, groupby='leiden_res0_2')

marker_genes_dict = {
    "Excitatory":["SNAP25","SLC17A7","CAMK2A"],"Inhibitory":["GAD1","GAD2","SLC6A1"],
    "Oligo":["MBP","MOBP","PLP1"],"OPC": ["PDGFRA","VCAN","CSPG4"],"Astro":["AQP4","GFAP","ALDH1L1"],
    "Micro":["CSF1R","C3","CD74"]
}
sc.pl.dotplot(adata, marker_genes_dict, "leiden_res0_2", dendrogram=True, save='dotplotMarkerGenes02.png')

sc.tl.dendrogram(adata, groupby='leiden_res0_2')

marker_genes_dict = {
    "Excitatory":["SNAP25","SLC17A7","CAMK2A"],"Inhibitory":["GAD1","GAD2","SLC6A1"],
    "Oligo":["MBP","MOBP","PLP1"],"OPC": ["PDGFRA","VCAN","CSPG4"],"Astro":["AQP4","GFAP","ALDH1L1"],
    "Micro":["CSF1R","C3","CD74"]
}
sc.pl.dotplot(adata, marker_genes_dict, "leiden_res0_5", dendrogram=True, save='dotplotMarkerGenes05.png')


sc.tl.dendrogram(adata, groupby='leiden_res0_2')

marker_genes_dict = {
    "Excitatory":["SNAP25","SLC17A7","CAMK2A"],"Inhibitory":["GAD1","GAD2","SLC6A1"],
    "Oligo":["MBP","MOBP","PLP1"],"OPC": ["PDGFRA","VCAN","CSPG4"],"Astro":["AQP4","GFAP","ALDH1L1"],
    "Micro":["CSF1R","C3","CD74"]
}
sc.pl.dotplot(adata, marker_genes_dict, "leiden_res1", dendrogram=True, save='dotplotMarkerGenes1.png')


#next

adata.write('/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/postClustering103124.h5ad')

# import scanpy as sc
# import pandas as pd
# import decoupler as dc
# import numpy as np

# import os
# import numpy as np
# import scanpy as sc
# import h5py

# root = '/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Sample_Five/Annotated_Filtered'

# # Run ORA to get marker enrichment scores
# dc.run_ora(
#     mat=adata,
#     net=markers_df,
#     source='cell_type',
#     target='gene',
#     min_n=3,
#     verbose=True,
#     use_raw=False
# )

# acts = dc.get_acts(adata, obsm_key='ora_estimate') #new anndata object with activity scores in its counts matrix (adata.X)

# # We need to remove inf and set them to the maximum value observed for pvals=0
# acts_v = acts.X.ravel()
# max_e = np.nanmax(acts_v[np.isfinite(acts_v)])
# acts.X[~np.isfinite(acts.X)] = max_e

# df = dc.rank_sources_groups(acts, groupby='leiden', reference='rest', method='t-test_overestim_var') #determines top predicted celltype per cluster

# #extract the top 3 predicted celltypes per cluster
# n_ctypes = 3
# ctypes_dict = df.groupby('group').head(n_ctypes).groupby('group')['names'].apply(lambda x: list(x)).to_dict() #Cluster number to celltypes dictionary

# sc.pl.matrixplot(acts, ctypes_dict, 'leiden', dendrogram=True, standard_scale='var',
# 	colorbar_title='Z-scaled scores', cmap='RdBu_r')

# annotation_dict = df.groupby('group').head(1).set_index('group')['names'].to_dict() #create dictionary of top 1 predicted cell type per cluster

# # Add cell type column based on annotation
# adata.obs['cell_type'] = [annotation_dict[clust] for clust in adata.obs['leiden']]

# # Visualize
# sc.pl.umap(adata, color='cell_type')
# annotatedUmap = sc.pl.umap(
#     adata,
#     color=["total_counts", "pct_counts_mt", "scDblFinder_score", "scDblFinder_class"],
#     save="annotatedUmapPostAnnotation.png"  # Use the 'save' parameter to save the figure
# )

# # Write the adata to the file
# adata.write_h5ad('/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/postAnnotation101524.h5ad')