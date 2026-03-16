
# import anndata as ad
# import scanpy as sc
# import os
# os.environ["HDF5_USE_FILE_LOCKING"] = "TRUE"
# import numpy as np

# import tempfile
# import zarr

# from functools import reduce
# import scanpy.external as sce

# import seaborn as sns
# import torch
# from rich import print
# from pathlib import Path
# import harmonypy as hm
# import matplotlib.pyplot as plt

# import pandas as pd

# # # #list all folders wanted

# libraries = ['200313-B22-B', '200225-B10-B', '200316-B24-B', '200305-B15-B', '190409-B5-B', '200317-B26-B', '200306-B16-B', '191219-B9-B', '200713-B33-B', '200312-B20-B', '201007-B57-B', '200810-B47-B', '200317-B27-B', '190403-B4-B', '200707-B30-B', '200226-B11-B', '200730-B41-B', '201007-B58-B', '200303-B14-B', '200701-B28-B', '200313-B23-B', '200715-B35-B', '200804-B42-B', '200708-B31-B', '200310-B18-B', '200810-B46-B', '200930-B55-B', '201022-B61-B', '201024-B59-B', '200728-B39-B', '201002-B56-B', '200702-B29-B', '200806-B44-B', '200311-B19-B', '191122-B6-R7090624-alone', '190403-B4-A', '191121-B6', '191213-B7-A', '191213-B7-B', '191217-B8-A', '191217-B8-B', '200309-B17-B', '200312-B21-B', '200316-B24-A', '200316-B25-A', '200317-B26-A', '200317-B27-A', '200707-B30-A', '200708-B31-A', '200714-B34-B', '200720-B36-A', '200720-B36-B', '200721-B37-A', '200722-B38-A', '200722-B38-B', '200729-B40-A', '200729-B40-B', '200804-B42-A', '200805-B43-A', '200805-B43-B', '200806-B44-A', '200810-B45-A', '200810-B46-A', '200810-B47-A', '200825-B48-A', '200825-B48-B', '200826-B49-A', '200908-B50-A', '200908-B50-B', '200909-B51-B', '200910-B52-B', '200915-B53-A', '200915-B53-B', '200916-B54-A', '201002-B56-A', '201021-B60-A', '201028-B62-A', '201207-B63-A', '201207-B63-B']

# # #Loading in files + cell assignment
# # mapping_csv = pd.read_csv("/om/scratch/Mon/mabdel03/SocialIsolation/newCellAnno.csv")
# mapping_csv = pd.read_csv("/net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/Preprocessed_Counts/cellAssignSocIsl.csv")
# mapping_csv.rename(columns={"Cell Barcode": "barcode"}, inplace=True)
# mapping_csv.rename(columns={"Assigned Patient": "patient_id"}, inplace=True)
# mapping_csv['barcode'] = mapping_csv['barcode'].astype(str)
# mapping_csv['Library'] = mapping_csv['Library'].astype(str)
# mapping_csv['patient_id'] = mapping_csv['patient_id'].astype(str)

# adataList = []
# patient_list=[]
# for library in libraries:
#     file_path = f'/om/scratch/Mon/mabdel03/SocialIsolation/{library}OutputP2.h5ad'
#     if os.path.exists(file_path):
#         adataL = ad.read_h5ad(file_path)
#         print(adataL)
#         # QC metrics
#         if 'pct_counts_mt' not in adataL.obs:
#             adataL.var['mt'] = adataL.var_names.str.startswith('MT-')
#             adataL.var['ribo'] = adataL.var_names.str.startswith(('RPS', 'RPL'))
#             adataL.var['hb'] = adataL.var_names.str.contains('HB')
#             sc.pp.calculate_qc_metrics(adataL, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)
#         adataL.obs['batch'] = library
#         if library == "60725338" or library == "MAP60725338":
#             # Special case — assign all cells to patient 2
#             adataL.obs['patient_id'] = '2'
#         else:
#             lib_map = mapping_csv[mapping_csv['Library'] == library].copy()
#             lib_map = lib_map.set_index('barcode')
#             adataL.obs.index = adataL.obs.index.astype(str)
#             adataL.obs['patient_id'] = adataL.obs.index.map(lib_map['patient_id'])
#         adataL = adataL[~adataL.obs['patient_id'].isna()].copy()
#         print("new")
#         print(adataL)
#         library_barcodes = set(mapping_csv[mapping_csv['Library'] == library]['barcode'])
#         adata_barcodes = set(adataL.obs_names)
#         overlap = adata_barcodes.intersection(library_barcodes)
#         for i in adataL.obs['patient_id'].unique():
#             if i not in patient_list:
#                 patient_list.append(i)
#             else:
#                 adataL = adataL[adataL.obs["patient_id"] != i].copy()
#         adataList.append(adataL)

# CSVPatient = pd.read_csv('dataset_652_basic_03-23-2022.csv')

# #subset for patients with this data
# adata = ad.concat(adataList, merge="same", index_unique="-{object_id}")#join="outer",

# qc_match = pd.read_csv("WGS_sample_QC_info.csv")
# unique_ones = adata.obs['patient_id'].unique()
# for idx, patient_id in adata.obs['patient_id'].items():
# 	if patient_id.startswith("SM"):
# 		ind = qc_match.index[qc_match.iloc[:, 1] == patient_id].tolist()[0]
# 		adata.obs.at[idx,'patient_id']=str(qc_match.iloc[ind, 0])
# 		print(adata.obs.at[idx,'patient_id'])

# ##Debug check!
# # unique_ones = adata.obs['patient_id'].unique()
# # for i in unique_ones:
# # 	print(i)
# # print(len(unique_ones))

# for i in adata.obs['patient_id'].unique():
# 	print(CSVPatient.loc[CSVPatient['projid'].astype(str) == str(i),'social_isolation_avg'])
# 	# print(CSVPatient.loc[CSVPatient['projid'].astype(str) == str(i), 'social_isolation_avg'])
# 	value = CSVPatient.loc[CSVPatient['projid'].astype(str) == str(i), 'social_isolation_avg']
# 	if not value.empty and pd.notna(float(value.values[0])):
# 	    print("YAY: " + str(i))
# 	else:
# 	    adata=adata[adata.obs['patient_id']!=i]

# adata.raw = ad.AnnData(X=adata.X.copy(), var=adata.var.copy())

# sc.pp.normalize_total(adata)  # per cell
# sc.pp.log1p(adata)  # log1p

# # finding highly variable features, and also making sure m20 genes are included
# sc.pp.highly_variable_genes(adata, flavor='seurat')#,n_top_genes=2000
# # gene_list = [
# #     "SCD", "LDLR", "STARD4", "ACAT2", "FDFT1", "LSS", "MMAB", "MVK",
# #     "DHCR24", "SREBF2", "INSIG1", "FDPS", "HMGCR", "MSMO1", "SQLE",
# #     "MVD", "IDI1", "HSD17B7", "HMGCS1", "ZNF730", "SC5D"
# # ]

# # hvg_mask = adata.var['highly_variable']
# # gene_mask = adata.var_names.isin(gene_list)
# # combined_mask = hvg_mask | gene_mask
# # if adata.raw is not None:
# #     adata.raw = adata.raw.to_adata()[:, combined_mask].copy()
# # adata = adata[:, combined_mask]

# sc.pp.scale(adata)

# # PCA
# sc.tl.pca(adata, svd_solver='arpack')#n_comps=30

# harmony_result = hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'batch')
# adata.obsm['X_harmony'] = harmony_result.Z_corr.T 

# # Neighbors, leiden clusters
# sc.pp.neighbors(adata,use_rep='X_harmony')#n_pcs=30
# sc.tl.leiden(adata, key_added="leiden_res0_2", resolution=0.2)
# sc.tl.leiden(adata, key_added="leiden_res0_5", resolution=0.5)
# sc.tl.leiden(adata, key_added="leiden_res1", resolution=1.0)

# # UMAP
# sc.tl.umap(adata)

# # UMAP Visualization
# sc.pl.umap(adata,
#            # color=["pct_counts_mt", "log1p_total_counts", "log1p_n_genes_by_counts","pct_counts_in_top_20_genes"],
#            color=["leiden_res0_2", "leiden_res0_5", "leiden_res1"],
#            legend_loc="on data",
#            save='umapLeidenTotalPatientFinal.png')


# # UMAP visualization
# sc.pl.umap(adata,
#            color=["pct_counts_mt", "log1p_total_counts", "log1p_n_genes_by_counts","pct_counts_in_top_20_genes"],
#            # color=["leiden_res0_2", "leiden_res0_5", "leiden_res1"],
#            legend_loc="on data",
#            save='umapLeidenTotalPatientQCFinal.png')

# # saving object
# adata.write('/om/scratch/Mon/mabdel03/SocialIsolation/totalAdata2.h5ad')

# # cell type assignment
# import pandas as pd
# import anndata as ad
# from scipy import sparse
# from sklearn.preprocessing import normalize
# import decoupler as dc
# import numpy as np
# import pyreadr

# import scanpy as sc
# from rpy2.robjects import r
# from rpy2.robjects import pandas2ri

# #creating and formatting markers DF

# pandas2ri.activate() # Load the RDS file
# rds_data = r['readRDS']('Brain_Human_PFC_Markers_Mohammadi2020.rds')
# markers_df = pandas2ri.rpy2py(rds_data)
# cell_type_names = list(markers_df.names)

# cell_types = list(markers_df)
# cell_type_names = [ct for ct in cell_type_names if ct != 'Ex_NRGN']

# data = []

# for cell_type_name, genes in zip(cell_type_names, cell_types):
# 	for gene in genes:
# 		data.append({'source': cell_type_name, 'target': gene, 'weight': 1.0})

# markers_df = pd.DataFrame(data)

# # load in h5ad
# adata = ad.read_h5ad('/om/scratch/Mon/mabdel03/SocialIsolation/totalAdata2.h5ad')

# adata = adata[adata.obs["leiden_res0_5"] != 0].copy()

# # adata.raw = adata
# if 'ora_estimate' in adata.obsm:
#     del adata.obsm['ora_estimate']

#  #Running ORA/whole process there

# dc.run_ora(adata,markers_df,source = "source",target = "target", use_raw=False)

# acts = dc.get_acts(adata, obsm_key='ora_estimate')

# # We need to remove inf and set them to the maximum value observed for pvals=0
# acts_v = acts.X.ravel()
# max_e = np.nanmax(acts_v[np.isfinite(acts_v)])
# acts.X[~np.isfinite(acts.X)] = max_e

# # acts_df = pd.DataFrame(acts.X, index=acts.obs_names, columns=acts.var_names)

# # acts_df.to_csv("enrichment_scoresDejager.csv")

# sc.pl.umap(acts, color=['Ast', 'leiden_res0_5'], cmap='RdBu_r')
# sc.pl.violin(acts, keys=['Ast'], groupby='leiden_res0_5')

# df = dc.rank_sources_groups(acts, groupby='leiden_res0_5', reference='rest', method='t-test_overestim_var')

# n_ctypes = 3
# ctypes_dict = df.groupby('group').head(n_ctypes).groupby('group')['names'].apply(lambda x: list(x)).to_dict()
# # with open("ctypes_dictDejager.txt", "w") as f:
# #     for key, value in ctypes_dict.items():
# #         f.write(f"{key}: {value}\n")

# sc.pl.matrixplot(acts, ctypes_dict, 'leiden_res0_5', dendrogram=True, standard_scale='var',
#                  colorbar_title='Z-scaled scores', cmap='RdBu_r')


# annotation_dict = df.groupby('group').head(1).set_index('group')['names'].to_dict()

# #Cell type insertion

# adata.obs['cell_type'] = [annotation_dict[clust] for clust in adata.obs['leiden_res0_5']]

# # visualizing
# sc.pl.umap(adata, color='cell_type')

# #formatting columns

# adata.obs['cell_type'] = adata.obs['cell_type'].str.replace(r'[-/]', '_', regex=True)

# if 'ora_pvals' in adata.obsm:
#     if isinstance(adata.obsm['ora_pvals'], pd.DataFrame):
#         adata.obsm['ora_pvals'].columns = [
#             col.replace('/', '_') for col in adata.obsm['ora_pvals'].columns
#         ]
#     elif isinstance(adata.obsm['ora_pvals'], dict):
#         adata.obsm['ora_pvals'] = {
#             key.replace('/', '_'): value for key, value in adata.obsm['ora_pvals'].items()
#         }

# adata.obsm['ora_estimate'] = adata.obsm['ora_estimate'].applymap(
#     lambda x: x.replace('/', '_') if isinstance(x, str) else x
# )
# adata.obsm['ora_estimate'].columns = [
#     col.replace('/', '_') for col in adata.obsm['ora_estimate'].columns
# ]



# adata.write('/om/scratch/Mon/mabdel03/SocialIsolation/totalAdataAnno2.h5ad')

#cell type graph
import anndata as ad
import scanpy as sc
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
# sc.pl.umap(adata,
#            # color=["pct_counts_mt", "log1p_total_counts", "log1p_n_genes_by_counts","pct_counts_in_top_20_genes"],
#            color=["cell_type"],
#            legend_loc="on data",
#            save='umapLeidenTotalPatientCellType.png')

adata_target = ad.read_h5ad('/om/scratch/Mon/mabdel03/SocialIsolation/totalAdataAnno2.h5ad')
patients = ["32383679","66406040","51400993"]

adata_target = adata_target[~adata_target.obs['patient_id'].isin(patients)]

for i in adata_target.obs['patient_id'].unique():
	print(adata_target[adata_target.obs['patient_id']==i])

adata_target = adata_target[~((adata_target.obs['patient_id'] == "50107583") & (adata_target.obs['cell_type'] == 'OPC'))].copy()

adata=adata_target.copy()
gene_list = [
    "SCD", "LDLR", "STARD4", "ACAT2", "FDFT1", "LSS", "MMAB", "MVK",
    "DHCR24", "SREBF2", "INSIG1", "FDPS", "HMGCR", "MSMO1", "SQLE",
    "MVD", "IDI1", "HSD17B7", "HMGCS1", "ZNF730", "SC5D"
]

hvg_mask = adata.var['highly_variable']
gene_mask = adata.var_names.isin(gene_list)
combined_mask = hvg_mask | gene_mask
if adata.raw is not None:
    adata.raw = adata.raw.to_adata()[:, combined_mask].copy()
adata = adata[:, combined_mask]

adata.write('/om/scratch/Mon/mabdel03/SocialIsolation/totalAdataAnno3HVG.h5ad')
adata_target.write('/om/scratch/Mon/mabdel03/SocialIsolation/totalAdataAnno3.h5ad')
# # PCA!
# if 'X_pca' not in adata_target.obsm:
#     print("Computing PCA...")
#     sc.tl.pca(adata_target)
# else:
#     print("PCA already computed.")

# # UMAP!
# if 'X_umap' not in adata_target.obsm:
#     print("Computing neighbors for UMAP...")
#     sc.pp.neighbors(adata_target, n_pcs=40)  # Adjust `n_pcs` according to your data
#     print("Computing UMAP...")
#     sc.tl.umap(adata_target)
# else:
#     print("UMAP already computed.")

# import math
# #UMAPs for each patient
# ncols = 6
# nrows = 6
# plots_per_image = ncols * nrows


# CSVPatient = pd.read_csv('dataset_652_basic_03-23-2022.csv')


# for i in adata_target.obs['patient_id'].unique():
# 	print(CSVPatient.loc[CSVPatient['projid'].astype(str) == str(i),'social_isolation_avg'])
# 	print(i)


# #batch plot UMAPs
# for i in range(0, len(adata_target.obs['patient_id'].unique()), plots_per_image):
#     fig, axes = plt.subplots(nrows, ncols, figsize=(20, 20))
#     axes = axes.flatten()
#     patient_chunk = adata_target.obs['patient_id'].unique()[i:i + plots_per_image]
#     for j, patient_id in enumerate(patient_chunk):
#         ax = axes[j]
#         patient_data = adata_target[adata_target.obs['patient_id'] == patient_id, :]
#         legend_option = 'on data' if j == 0 else False
#         sc.pl.umap(
#             patient_data,
#             color='leiden_res0_2',
#             title=str(patient_id),
#             ax=ax,
#             show=False,
#             legend_loc=legend_option,
#         )
#     for k in range(len(patient_chunk), len(axes)):
#         fig.delaxes(axes[k])
    
#     plt.tight_layout()
#     plt.savefig(f"umap_tile_patients_{i // plots_per_image + 1}.png", dpi=300)
#     plt.close()