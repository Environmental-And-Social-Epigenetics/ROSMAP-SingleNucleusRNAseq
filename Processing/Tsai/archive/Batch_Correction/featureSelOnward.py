
import anndata as ad
import scanpy as sc
import os
import anndata2ri
import numpy as np

#load in batch correction obj where i cleared out extra layers and used normalized scvi matrix as X
file_path = '/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/clearedObj.h5ad' 
# file_path = '/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/postFeatureSel3.h5ad'
adata = ad.read_h5ad(file_path)


#feature selection
sc.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes=5000)

#saving object
adata.write('/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/postFeatureSel3101524.h5ad')

adata = adata[:, adata.var.highly_variable]

#scale, pca, neighbors, leiden
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=100, n_pcs=20)
sc.tl.leiden(adata,resolution=0.03,n_iterations=20)
adata.write('/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/postPCA101524.h5ad')



















#IGNORE
# sc.pl.umap(adata, color="total_counts")

# adata.write('/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/postUMAP.h5ad')



# import rpy2.robjects as ro
# from rpy2.robjects import numpy2ri, pandas2ri
# from rpy2.robjects.packages import importr

# # Activate numpy and pandas conversion for rpy2
# numpy2ri.activate()
# pandas2ri.activate()

# # Load the AnnData object
# file_path = '/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/optimized_data092324.h5ad'
# adata = ad.read_h5ad(file_path)

# ro.pandas2ri.activate()
# anndata2ri.activate()

# deviance_pkg = importr('scry')

# # Ensure the AnnData object uses compatible data types
# obj_name = f'adata'
# ro.globalenv[obj_name] = adata

# # Use R magic command to run devianceFeatureSelection


# ro.r('''
# library(SingleCellExperiment)
# sce <- devianceFeatureSelection(get(obj_name), assay="X")
# binomial_deviance <- rowData(sce)$binomial_deviance
# ''')

# # Convert binomial_deviance to numpy array
# binomial_deviance = np.array(ro.r['binomial_deviance'])

# idx = binomial_deviance.argsort()[-4000:] 
# mask = np.zeros(adata.var_names.shape, dtype=bool)
# mask[idx] = True

# adata.var["highly_deviant"] = mask
# adata.var["binomial_deviance"] = binomial_deviance

# sc.pp.highly_variable_genes(adata, layer="X")

# # Ensure observation names are unique
# # obj.obs_names_make_unique()

# # # Remove other matrices and keep only the counts matrix
# # obj.layers.clear()
# # obj.raw = None
# # obj.obsm.clear()
# # obj.varm.clear()
# # obj.obsp.clear()

# # obj.varp.clear()





# # obj.X = obj.X.astype(np.float32)

# # Convert to sparse if it's not already sparse
# # if not sp.issparse(obj.X):

# # obj.write_h5ad('optimized_data092324.h5ad')




# # Extract the expression matrix as a numpy array


# # Perform the standard scanpy processing

# # Save the object after feature selection
# root = '/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/'
# name = 'DejagerLibrarySubjectPostFeatureSel.h5ad'
# adata.write_h5ad(os.path.join(root, name))


# # Continue with PCA, UMAP, clustering, and saving the updated object
# sc.settings.verbosity = 0
# sc.settings.set_figure_params(dpi=80, facecolor="white", frameon=False)

# # Set highly variable genes
# adata.var["highly_variable"] = adata.var["highly_deviant"]

# # Perform PCA
# sc.pp.pca(adata, svd_solver="arpack", use_highly_variable=True)
# sc.pl.pca_scatter(adata, color="total_counts")

# # Perform neighbors and UMAP
# sc.pp.neighbors(adata)
# sc.tl.umap(adata)
# sc.pl.umap(adata, color="total_counts")

# # Save the pre-dimension reduction filtered object
# name = 'DejagerLibrarySubjectOGPreDimRedFilter.h5ad'
# adata.write_h5ad(os.path.join(root, name))

# # Perform clustering
# sc.pp.neighbors(adata, n_pcs=30)
# sc.tl.umap(adata)

# sc.tl.leiden(adata, key_added="leiden_res0_25", resolution=0.25)
# sc.tl.leiden(adata, key_added="leiden_res0_5", resolution=0.5)
# sc.tl.leiden(adata, key_added="leiden_res1", resolution=1.0)

# # Plot UMAP with Leiden clusters
# sc.pl.umap(
#     adata, 
#     color=["leiden_res0_25", "leiden_res0_5", "leiden_res1"], 
#     legend_loc="on data"
# )

# # Save the object after clustering
# name = 'DejagerLibrarySubjectPostClustering.h5ad'
# adata.write_h5ad(os.path.join(root, name))

# Optionally save new_obj if necessary
# root = '/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/'
# name = 'DejagerLibrarySubjectPostDimRedFiltering.h5ad'
# new_obj.write_h5ad(os.path.join(root, name))
