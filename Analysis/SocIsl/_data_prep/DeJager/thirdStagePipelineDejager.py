
import anndata as ad
import scanpy as sc
import os
os.environ["HDF5_USE_FILE_LOCKING"] = "TRUE"
import numpy as np

import tempfile
import zarr

from functools import reduce
import scanpy.external as sce

import seaborn as sns
import torch
from rich import print
from pathlib import Path
import harmonypy as hm
import matplotlib.pyplot as plt

import pandas as pd

# # # Replace 'your_directory_path' with the path to the directory you want to scan
# directory_path = '/om/scratch/Fri/mabdel03/Subfolder/'

# # List all folders in the directory
# # libraries = [folder for folder in os.listdir("/net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/Preprocessed_Counts/") if os.path.isdir(os.path.join("/net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/Preprocessed_Counts/", folder))]
# # remove_libraries = {
# #     '200930-B55-B', '200311-B19-B', '200701-B28-B', '200810-B46-B',
# #     '200708-B31-B', '201022-B61-B', '200806-B44-B', '200715-B35-B',
# #     '200804-B42-B', '201007-B58-B', '200313-B23-B', '200702-B29-B',
# #     '200303-B14-B', '201024-B59-B', '200730-B41-B', '200310-B18-B',
# #     '200728-B39-B'
# # }

libraries = [
    '200916-B54-A', '200313-B22-B', '200225-B10-B', '200316-B24-B', '200810-B45-A',
    '200305-B15-B', '190409-B5-B', '200317-B26-B', '200306-B16-B', '191219-B9-B',
    '200713-B33-B', '200312-B20-B', '200309-B17-B', '201007-B57-B', '200810-B47-B',
    '200317-B27-B', '190403-B4-B', '200707-B30-B', '200226-B11-B', '200730-B41-B',
    '201007-B58-B', '200303-B14-B', '200701-B28-B', '200313-B23-B', '200715-B35-B',
    '200804-B42-B', '200708-B31-B', '200310-B18-B', '200810-B46-B', '200930-B55-B',
    '201022-B61-B', '201024-B59-B', '200728-B39-B', '201002-B56-B', '200702-B29-B',
    '200806-B44-B', '200311-B19-B', '191122-B6-R7090624-alone'
]

# # # # libraries = [lib for lib in libraries if lib not in remove_libraries]

# # # # libraries.remove("qc_plots")
# # # # libraries.remove("figures")
# # # # libraries.remove("concatObjFinal.zarr")
# # # # libraries.remove("OutputP2")
# # # # count=0
# # # mapping_csv = pd.read_csv("cell_to_patient_assignments5.csv")
# # # mapping_csv.rename(columns={"Cell Barcode": "barcode"}, inplace=True)
# # # mapping_csv.rename(columns={"Assigned Patient": "patient_id"}, inplace=True)
# # # mapping_csv['barcode'] = mapping_csv['barcode'].astype(str)
# # # mapping_csv['Library'] = mapping_csv['Library'].astype(str)
# # # mapping_csv['patient_id'] = mapping_csv['patient_id'].astype(str)

mapping_csv = pd.read_csv("/net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/Preprocessed_Counts/cell_to_patient_assignments5.csv")
header_row = pd.DataFrame([mapping_csv.columns.tolist()], columns=mapping_csv.columns)
mapping_csv = pd.concat([header_row, mapping_csv], ignore_index=True)
# mapping_csv.rename(columns={"Cell Barcode": "barcode"}, inplace=True)
# mapping_csv.rename(columns={"Assigned Patient": "patient_id"}, inplace=True)
mapping_csv.columns = ["barcode","patient_id","Library"]
mapping_csv['barcode'] = mapping_csv['barcode'].astype(str)
mapping_csv['Library'] = mapping_csv['Library'].astype(str)
mapping_csv['patient_id'] = mapping_csv['patient_id'].astype(str)


mapping_csv2 = pd.read_csv("/net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/Preprocessed_Counts/cell_to_patient_assignments7.csv")
header_row2 = pd.DataFrame([mapping_csv2.columns.tolist()], columns=mapping_csv2.columns)
mapping_csv2 = pd.concat([header_row2, mapping_csv2], ignore_index=True)
# mapping_csv.rename(columns={"Cell Barcode": "barcode"}, inplace=True)
# mapping_csv.rename(columns={"Assigned Patient": "patient_id"}, inplace=True)
mapping_csv2.columns = ["barcode","patient_id","Library"]
mapping_csv2['barcode'] = mapping_csv2['barcode'].astype(str)
mapping_csv2['Library'] = mapping_csv2['Library'].astype(str)
mapping_csv2['patient_id'] = mapping_csv2['patient_id'].astype(str)
existing_libs = set(mapping_csv["Library"])
mapping_csv2 = mapping_csv2[~mapping_csv2["Library"].isin(existing_libs)]


mapping_csv3 = pd.read_csv("/net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/Preprocessed_Counts/cell_to_patient_assignments12.csv")
header_row3 = pd.DataFrame([mapping_csv3.columns.tolist()], columns=mapping_csv3.columns)
mapping_csv3 = pd.concat([header_row3, mapping_csv3], ignore_index=True)
# mapping_csv.rename(columns={"Cell Barcode": "barcode"}, inplace=True)
# mapping_csv.rename(columns={"Assigned Patient": "patient_id"}, inplace=True)
mapping_csv3.columns = ["barcode","patient_id","Library"]
mapping_csv3['barcode'] = mapping_csv3['barcode'].astype(str)
mapping_csv3['Library'] = mapping_csv3['Library'].astype(str)
mapping_csv3['patient_id'] = mapping_csv3['patient_id'].astype(str)
existing_libs = set(mapping_csv2["Library"])
mapping_csv3 = mapping_csv3[~mapping_csv3["Library"].isin(existing_libs)]
existing_libs = set(mapping_csv["Library"])
mapping_csv3 = mapping_csv3[~mapping_csv3["Library"].isin(existing_libs)]


mapping_csv = pd.concat([mapping_csv, mapping_csv2, mapping_csv3], ignore_index=True)

adataList = []
patient_list=[]


for library in libraries:
    print(library)
    # root = "/net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/Preprocessed_Counts/"+library
    # libraryPath = os.path.join(root, 'processed_feature_bc_matrix_filtered.h5')
    libraryPath="/om/scratch/Fri/mabdel03/Subfolder/"+library+"OutputP2.h5ad"
    if os.path.exists(libraryPath):
       library_adata =  ad.read_h5ad(libraryPath)
        # library_adata = sc.read_10x_h5(libraryPath)
        original_adata = library_adata
        if len(library_adata.obs_names) !=0:
            library_adata.var_names_make_unique()
            if library == "191122-B6-R7090624-alone":
                # Special case — assign all cells to patient 2
                library_adata.obs['patient_id'] = '2'
            else:
                lib_map = mapping_csv[mapping_csv['Library'] == library].copy()
                lib_map = lib_map.set_index('barcode')
                library_adata.obs.index = library_adata.obs.index.astype(str)
                library_adata.obs['patient_id'] = library_adata.obs.index.map(lib_map['patient_id'])
            library_adata = library_adata[~library_adata.obs['patient_id'].isna()].copy()
            # mitochondrial genes
            # library_adata.var["mt"] = library_adata.var_names.str.startswith("MT-")
            # library_adata.var["ribo"] = library_adata.var_names.str.startswith(("RPS", "RPL"))
            # library_adata.var["hb"] = library_adata.var_names.str.contains(("^HB[^(P)]"))
            # # QC Metrics
            # sc.pp.calculate_qc_metrics(library_adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)
            # metrics = ["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt","pct_counts_in_top_20_genes"]
            for i in library_adata.obs['patient_id'].unique():
                if i not in patient_list:
                    patient_list.append(i)
                else:
                    library_adata = library_adata[library_adata.obs["patient_id"] != i].copy()
            library_adata.obs['batch'] =library_adata.obs['patient_id']
            print(library_adata)
            adataList.append(library_adata)


# # # Select a single patient by batch code
# # # selected_patient = adata.obs['batch'].unique()[0]  # Replace with your specific patient code
# # # adata = adata[adata.obs['batch'] == selected_patient].copy()
for adataL in adataList:
    adataL.obs['patient_id'] = adataL.obs['patient_id'].astype(str)

adata = ad.concat(adataList, merge="same")


qc_match = pd.read_csv("WGS_sample_QC_info.csv")
unique_ones = adata.obs['patient_id'].unique()
for idx, patient_id in adata.obs['patient_id'].items():
    if patient_id.startswith("SM"):
        ind = qc_match.index[qc_match.iloc[:, 1] == patient_id].tolist()[0]
        adata.obs.at[idx,'patient_id']=str(qc_match.iloc[ind, 0])

adata.obs_names_make_unique()

desiredPatients =[10205244, 15218541, 15176592, 11165535, 10246987, 11353057, 10277308,10394182, 50405330, 10202345, 11200645, 10490993, 50101659, 70883578, 50405042, 10218339, 10101327, 10101589, 10248033, 10684501, 11157783, 11326252, 11327005,11606935, 18414513, 21362537, 22789958, 43485807, 44019405, 66754397,90267190, 50402693,21000180, 20252720, 20933324, 21176459, 20173942, 21293107, 98953007, 20358955, 20500815, 21000630, 81852640, 20242958, 20152393, 54122640, 14641578, 81874628, 20898476, 89614402, 60725338, 32383679, 20983400, 21157370, 20153984, 20236838,21411459, 63874408, 50106578, 50108886, 20221273, 20254588, 21140119, 61827429,75675336, 69432088,2]
    string_list = [str(item) for item in desiredPatients]

adata=adata[adata.obs['patient_id'].isin(string_list)]
print("AD!!!!")
print(adata)
# adata =  ad.read_h5ad('/om/scratch/Fri/mabdel03/Subfolder/totalAdata040825.h5ad')
# adata.X = adata.raw.X.copy()
# adata2 = ad.read_h5ad('/om/scratch/Fri/mabdel03/Subfolder/totalAdataAnno040825.h5ad')
# adata = adata[:, adata.var['highly_variable']]
# adata2 = adata2[:, adata2.var['highly_variable']]
# adata2 = adata2[adata2.obs['cell_type'] != 'Ex_NRGN']
# adata = adata[adata.obs_names.isin(adata2.obs_names)].copy()
# adata.raw = ad.AnnData(X=adata.X.copy(), var=adata.var.copy())

# adata2.write_h5ad('/om/scratch/Fri/mabdel03/Subfolder/totalAdataAnno040825.h5ad')
# adata.write_h5ad('/om/scratch/Fri/mabdel03/Subfolder/totalAdataAnno040825.h5ad')
# print(adata.obs['patient_id'].value_counts())
# print("Shape2")
# print(adata.shape)
# # Load the CSV file (adjust delimiter if needed)

# mapping_csv = pd.read_csv("/net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/Preprocessed_Counts/cellAssignSocIsl.csv")

# mapping_csv.rename(columns={"Cell Barcode": "barcode"}, inplace=True)

# cells = mapping_csv.iloc[:, 0].values
# cells = cells.astype(str)
# patients = mapping_csv.iloc[:, 1].values
# libraries =  mapping_csv.iloc[:, 2].values
# adata.obs_names = adata.obs_names.str.replace(r'-1-\{object_id\}\d+$', '-1', regex=True)


# adata.obs['patient_id'] = np.full(len(adata.obs_names), "", dtype=str)

# for i in range(0,len(adata.obs_names)):
#   assignedPatient=""
#   matches = np.char.find(cells, str(adata.obs_names[i])) >= 0
#   possibleIndices = np.where(matches)[0]
#   # print(possibleIndices)
#   for k in possibleIndices:
#     if str(libraries[k]) == str(adata.obs.loc[adata.obs_names[i], "batch"]):
#      assignedPatient=str(patients[k])
#      break
#   adata.obs.loc[adata.obs_names[i], 'patient_id'] = assignedPatient

# adata.obs_names_make_unique()
# adata = adata[adata.obs['patient_id'] != '', :]


# print("patients")
# numbers = [
#     482428, 2525608, 3430444, 7265221, 7253015, 8109170, 10100574, 10100862, 10222853, 
#     10221262, 10222853, 10277308, 10516762, 11259716, 11392518, 11631558, 14452889, 
#     15115927, 15844425, 15712086, 17219510, 20105242, 20124321, 20139850, 20156469, 
#     20237131, 20225925, 20240514, 20267709, 20283691, 20348895, 20380417, 20594407, 
#     20646778, 20865035, 20892128, 20906493, 20911508, 21145740, 21145876, 21183160, 
#     21403733, 23004922, 22868024, 24747976, 35072859, 35286551, 35286551, 36492755, 
#     39484737, 43074402, 50104134, 50108048, 50108750, 50109927, 50109477, 50301675, 
#     51400993, 50405330, 51668135, 51520126, 60725338, 77143621, 76867532, 72650337, 
#     78353027, 83984043, 84653463, 84732827, 82624422, 87410220, 86903794
# ]
# numbers = [
#     33411712,
#     50305165,
#     20225925,
#     20798913,
#     50104846,
#     20453625,
#     20500815,
#     10478041,
#     10101291,
#     91018909,
#     6073025,
#     10490993,
#     50302978
# ]

# adata = adata[~adata.obs['patient_id'].isin(numbers), :]

# # #filter by patient id
adata.raw = ad.AnnData(X=adata.X.copy(), var=adata.var.copy())


# # Add gene names as a column if needed
# # adataK.var['gene_name'] = adataK.var.index
df = pd.DataFrame(adata.var_names, columns=["gene"])

# # # # Save as CSV
df.to_csv("var_names.csv", index=False)
# # adata.raw.var = adata.var.copy()  # Copy gene metadata to raw
# # adata.raw.obs = adata.obs.copy()  # Copy obs metadata to raw


# # adataListA = []
# # for library in libraries:
# # 	file_path = '/net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/Preprocessed_Counts/'+str(library)+'OutputP2A.h5ad'
# # 	if os.path.exists(file_path):
# # 		adataL = ad.read_h5ad(file_path)
# # 		if 'pct_counts_mt' not in adataL.obs:
# # 			adataL.var['mt'] = adataL.var_names.str.startswith('MT-')
# # 			adataL.var['ribo'] = adataL.var_names.str.startswith(('RPS', 'RPL'))
# # 			adataL.var['hb'] = adataL.var_names.str.contains('HB')  # If hemoglobin genes are relevant
# # 			sc.pp.calculate_qc_metrics(adataL, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)
# # 		adataL.obs['batch'] = library
# # 		adataListA.append(adataL)


# # Select a single patient by batch code
# # selected_patient = adata.obs['batch'].unique()[0]  # Replace with your specific patient code
# # adata = adata[adata.obs['batch'] == selected_patient].copy()

# adataA = ad.concat(adataListA, merge="same", index_unique="-{object_id}")#join="outer",


# adataListB = []
# for library in libraries:
# 	file_path = '/net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/Preprocessed_Counts/'+str(library)+'OutputP2B.h5ad'
# 	if os.path.exists(file_path):
# 		adataL = ad.read_h5ad(file_path)
# 		if 'pct_counts_mt' not in adataL.obs:
# 			adataL.var['mt'] = adataL.var_names.str.startswith('MT-')
# 			adataL.var['ribo'] = adataL.var_names.str.startswith(('RPS', 'RPL'))
# 			adataL.var['hb'] = adataL.var_names.str.contains('HB')  # If hemoglobin genes are relevant
# 			sc.pp.calculate_qc_metrics(adataL, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)
# 		adataL.obs['batch'] = library
# 		adataListB.append(adataL)


# # Select a single patient by batch code
# # selected_patient = adata.obs['batch'].unique()[0]  # Replace with your specific patient code
# # adata = adata[adata.obs['batch'] == selected_patient].copy()

# adataB = ad.concat(adataListB, merge="same", index_unique="-{object_id}")#join="outer",

# adataListC = []
# for library in libraries:
# 	file_path = '/net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/Preprocessed_Counts/'+str(library)+'OutputP2C.h5ad'
# 	if os.path.exists(file_path):
# 		adataL = ad.read_h5ad(file_path)
# 		if 'pct_counts_mt' not in adataL.obs:
# 			adataL.var['mt'] = adataL.var_names.str.startswith('MT-')
# 			adataL.var['ribo'] = adataL.var_names.str.startswith(('RPS', 'RPL'))
# 			adataL.var['hb'] = adataL.var_names.str.contains('HB')  # If hemoglobin genes are relevant
# 			sc.pp.calculate_qc_metrics(adataL, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)
# 		adataL.obs['batch'] = library
# 		adataListC.append(adataL)


# # Select a single patient by batch code
# # selected_patient = adata.obs['batch'].unique()[0]  # Replace with your specific patient code
# # adata = adata[adata.obs['batch'] == selected_patient].copy()

# adataC = ad.concat(adataListC, merge="same", index_unique="-{object_id}")#join="outer",

# Assuming adata.X contains raw counts initially


# Normalizing data
# sc.pp.normalize_total(adataA)  # per cell
# sc.pp.log1p(adataA)  # log1p

# sc.pp.normalize_total(adataB)  # per cell
# sc.pp.log1p(adataB)  # log1p

# sc.pp.normalize_total(adataC)  # per cell
# sc.pp.log1p(adataC)  # log1p
# adata = adata[:, (adata.X != 0).sum(axis=0) > 0]

sc.pp.normalize_total(adata)  # per cell
sc.pp.log1p(adata)  # log1p
print("norm")

# # Finding highly variable features
# # sc.pp.highly_variable_genes(adataA, flavor='seurat')#,n_top_genes=2000
# # adataA = adataA[:, adataA.var['highly_variable']]  # subsetting on this basis
# # sc.pp.highly_variable_genes(adataB, flavor='seurat')#,n_top_genes=2000
# # adataB = adataB[:, adataB.var['highly_variable']]  # subsetting on this basis
# # sc.pp.highly_variable_genes(adataC, flavor='seurat')#,n_top_genes=2000
# # adataC = adataC[:, adataC.var['highly_variable']]  # subsetting on this basis
sc.pp.highly_variable_genes(adata, flavor='seurat')#,n_top_genes=2000
print("highlyvar")

# # adataK = adataK[:, adata.var['highly_variable']]  # subsetting on this basis
# # adata.raw=adataK

# # Scaling data
# # sc.pp.scale(adataA)
# # sc.pp.scale(adataB)
# # sc.pp.scale(adataC)
sc.pp.scale(adata)

# # PCA

sc.tl.pca(adata, svd_solver='arpack')#n_comps=30
print("pca")

# # print(adata.obs['patient_id'].value_counts())

# # print(adata.obs['patient_id'].value_counts())
# # print(adata.obs['patient_id'].value_counts())
# # print(adata.obs['patient_id'].value_counts())
# # #batch correction

# # print(adata.shape)
# # print(adata.obs['patient_id'].value_counts())
# # print(adata.obs['patient_id'].isna().sum())
# # print(np.isnan(adata.obsm['X_pca']).sum())
# # print(np.all(adata.obsm['X_pca'] == 0))



harmony_result = hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'patient_id')
adata.obsm['X_harmony'] = harmony_result.Z_corr.T 

# Neighbors, leiden clusters
sc.pp.neighbors(adata, use_rep='X_harmony')
sc.tl.leiden(adata, key_added="leiden_res0_2", resolution=0.2)
sc.tl.leiden(adata, key_added="leiden_res0_5", resolution=0.5)
sc.tl.leiden(adata, key_added="leiden_res1", resolution=1.0)

print("neigh")
# # UMAP

sc.tl.umap(adata)
print("Shape3")
print(adata.shape)


# UMAP Visualization
sc.pl.umap(adata,
           # color=["pct_counts_mt", "log1p_total_counts", "log1p_n_genes_by_counts","pct_counts_in_top_20_genes"],
           color=["leiden_res0_2", "leiden_res0_5", "leiden_res1"],
           legend_loc="on data",
           save='umapLeidenTotalPatientFinal.png')


# UMAP Visualization
sc.pl.umap(adata,
           color=["pct_counts_mt", "log1p_total_counts", "log1p_n_genes_by_counts","pct_counts_in_top_20_genes"],
           # color=["leiden_res0_2", "leiden_res0_5", "leiden_res1"],
           legend_loc="on data",
           save='umapLeidenTotalPatientQCFinal.png')
print("UNIIIWIWIE")
print(adata.obs['leiden_res0_5'].unique())


# #saving object
adata.write('/om/scratch/Fri/mabdel03/Subfolder/totalAdata040825.h5ad')

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
cell_type_names = list(markers_df.names)

cell_types = list(markers_df)
cell_type_names = [ct for ct in cell_type_names if ct != 'Ex_NRGN']

data = []

for cell_type_name, genes in zip(cell_type_names, cell_types):
	for gene in genes:
		data.append({'source': cell_type_name, 'target': gene, 'weight': 1.0})

markers_df = pd.DataFrame(data)

# load in h5ad
if 'ora_estimate' in adata.obsm:
    del adata.obsm['ora_estimate']

 #Running ORA/whole process there

dc.run_ora(adata,markers_df,source = "source",target = "target", use_raw=False)

acts = dc.get_acts(adata, obsm_key='ora_estimate')

sc.pl.umap(acts, color=['Ast', 'leiden_res0_5'], cmap='RdBu_r')
sc.pl.violin(acts, keys=['Ast'], groupby='leiden_res0_5')

df = dc.rank_sources_groups(acts, groupby='leiden_res0_5', reference='rest', method='t-test_overestim_var')

n_ctypes = 3
ctypes_dict = df.groupby('group').head(n_ctypes).groupby('group')['names'].apply(lambda x: list(x)).to_dict()
# with open("ctypes_dictDejager.txt", "w") as f:
#     for key, value in ctypes_dict.items():
#         f.write(f"{key}: {value}\n")
print("ANUEUR")
acts_v = acts.X.ravel()
max_e = np.nanmax(acts_v[np.isfinite(acts_v)])
acts.X[~np.isfinite(acts.X)] = max_e

print(np.any(~np.isfinite(acts.X)))

group_means = acts.to_df().groupby(acts.obs['leiden_res0_5']).mean()
print(group_means.isna().sum().sum())

corr = np.corrcoef(group_means)
print(np.allclose(corr, corr.T, equal_nan=True))

# sc.pl.matrixplot(acts, ctypes_dict, 'leiden_res0_5', dendrogram=True, standard_scale='var',
#                  colorbar_title='Z-scaled scores', cmap='RdBu_r')


annotation_dict = df.groupby('group').head(1).set_index('group')['names'].to_dict()

#Cell type insertion

adata.obs['cell_type'] = [annotation_dict[clust] for clust in adata.obs['leiden_res0_5']]

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



adata.write('/om/scratch/Fri/mabdel03/Subfolder/totalAdataAnno040825.h5ad')

import anndata as ad
import scanpy as sc
import os
import numpy as np
import matplotlib.pyplot as plt

# # source + target adata - for batch info transfer

# import anndata as ad
# import scanpy as sc
# import os
# import numpy as np
# import matplotlib.pyplot as plt

#
adata=ad.read_h5ad('/om/scratch/Fri/mabdel03/Subfolder/totalAdataAnno040825.h5ad')
# file_path_target = '/om/scratch/Fri/mabdel03/Subfolder/totalAdataAnno040825.h5ad'


# # adata = ad.read_h5ad(file_path_target)

# # adata = adata[:, adata.var['highly_variable']]

# adata.write('/om/scratch/Fri/mabdel03/Subfolder/totalAdataAnno040825.h5ad')


# print(adata)

# patient_ids = [21000180, 20252720, 20933324, 21293107, 98953007, 20358955, 21000630, 81852640,20152393, 50104846, 20254902, 54122640, 20195344, 81874628, 20706215, 20120255,20232474, 89614402, 60725338, 32383679, 20983400, 21157370, 20153984, 20236838,21411459, 63874408, 50106578, 50108886, 20221273, 20254588, 21140119, 61827429,75675336, 69432088]
# patient_ids=[str(i) for i in patient_ids]
# print(adata.obs['patient_id'].unique())
# adata=adata[adata.obs['patient_id'].isin(patient_ids)]

# sc.pl.umap(adata,
#            # color=["pct_counts_mt", "log1p_total_counts", "log1p_n_genes_by_counts","pct_counts_in_top_20_genes"],
#            color=["cell_type"],
#            legend_loc="on data",
#            save='umapLeidenTotalPatientCellTypeDejagerFemale.png')


# # # UMAP Visualization
# sc.pl.umap(adata,
#            color=["pct_counts_mt", "log1p_total_counts", "log1p_n_genes_by_counts","pct_counts_in_top_20_genes"],
#            # color=["leiden_res0_2", "leiden_res0_5", "leiden_res1"],
#            legend_loc="on data",
#            save='umapLeidenTotalPatientQCFinalFemale.png')

# adata = ad.read_h5ad(file_path_target)
# print(adata)

# patient_ids = [10205244, 15218541, 15176592, 11165535, 10246987, 11353057, 10277308, 10394182, 50405330, 10202345, 11200645, 10490993, 50101659, 70883578, 50405042, 10218339,46547648, 43485807, 44019405, 50402693, 50108912, 82317494, 11157783, 10248033,11327005, 46000440, 50106730, 21362537, 90267190, 10684501, 66754397, 10101327]
# patient_ids=[str(i) for i in patient_ids]
# print(adata.obs['patient_id'].unique())
# adata=adata[adata.obs['patient_id'].isin(patient_ids)]

# sc.pl.umap(adata,
#            # color=["pct_counts_mt", "log1p_total_counts", "log1p_n_genes_by_counts","pct_counts_in_top_20_genes"],
#            color=["cell_type"],
#            legend_loc="on data",
#            save='umapLeidenTotalPatientCellTypeDejagerMale.png')


# # # UMAP Visualization
# sc.pl.umap(adata,
#            color=["pct_counts_mt", "log1p_total_counts", "log1p_n_genes_by_counts","pct_counts_in_top_20_genes"],
#            # color=["leiden_res0_2", "leiden_res0_5", "leiden_res1"],
#            legend_loc="on data",
#            save='umapLeidenTotalPatientQCFinalMale.png')


# # sc.pl.umap(adata_target,color=["cell_type"],legend_loc="on data",save='umapumapLeidenTotalPatientCellTypeFinal.png')
# # sc.pl.umap(adata_target,color=["patient_id"],legend_loc="on data",save='umapumapLeidenTotalPatientByPatient.png')

# # get and align batch info!!!
# # sample_series = adata_source.obs['batch'][adata_source.obs_names.isin(adata_target.obs_names)]
# # adata_target.obs['batch'] = sample_series.reindex(adata_target.obs_names)

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