
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


libraries = ['482428', '2518573', '2899847', '3283241', '3713990', '3889845', '6107196', '6804844', '7253015', '7265221', '7311370', '8109170', '8132197', '9650662', '9841821', '14184286', '14498577', '15218541', '16068769', '17929065', '18414513', '18659212', '18920002', '20907534', '21362537', '22396591', '22776575', '22789958', '22868024', '23004922', '23690880', '24039289', '24680888', '24747976', '25300551', '26569730', '26631069', '27586957', '29933130', '31509843', '31726180', '32383679', '32705437', '33411712', '33501827', '34542628', '34726040', '34962204', '35941263', '36492755', '36830117', '37030589', '37065652', '37178462', '37436329', '37527863', '38264019', '39484737', '39989287', '41285665', '41773404', '42988567', '43074402', '43485807', '44019405', '44299049', '44671043', '45115248', '45566083', '46000440', '46246604', '46251007', '48023497', '49676048', '50100806', '50100932', '50101523', '50102114', '50102376', '50103129', '50105301', '50105699', '50106280', '50106442', '50106730', '50107583', '50107619', '50107745', '50107907', '50108886', '50108912', '50109477', '50109639', '50109927', '50111971', '50197261', '50300084', '50300408', '50300822', '50301099', '50301125', '50301387', '50301675', '50302266', '50302554', '50400259', '50400385', '50400835', '50401002', '50402567', '50402693', '50402729', '50403446', '50404299', '50405330', '50406057', '50409406', '50410319', '50500136', '50500550', '51520126', '51624179', '51815338', '51826377', '52311825', '52940446', '53772202', '54324848', '57978756', '60725338', '60747316', '60848460', '61029627', '61142759', '61344957', '61827429', '62404688', '62574447', '63188799', '64291829', '64336939', '64493027', '65001607', '65002324', '65292866', '65499271', '65652206', '65736039', '66406040', '66754397', '66924745', '67185070', '67429065', '68015667', '68525196', '68539908', '68745332', '69982533', '70153803', '70625336', '70816595', '70883578', '70917649', '71626373', '71648351', '72650337', '72777797', '74284255', '74753465', '75169847', '75990666', '76733461', '76867532', '77143621', '77239958', '77891596', '78353027', '78452313', '80790863', '81086436', '82317494', '82684406', '83034844', '83173570', '83984043', '84417209', '84642424', '84653463', '85065193', '85171938', '86934089', '87410220', '87645445', '87779516', '89001223', '89164957', '89546375', '89903942', '90267190', '90780976', '90942860', '91018909', '91444029', '91707643', '92393245', '92629514', '93815598', '94092977', '94144536', '94430339', '94974890', '95491648', '96095092', '98953007', '86408244', '50102088']
# batch_deidentified = pd.read_csv("batch_mapping_deidentified.tsv", sep = '\t')
# batchiden = pd.read_csv('full_batch_mapping.tsv', sep='\t')

# opc = ad.read_h5ad("opc.h5ad")
# subjects = opc.obs['subject'].unique()

# subjectMatches = {}

# for i in subjects:
# 	result = batch_deidentified[batch_deidentified["subject"] == i]
# 	subjectMatches[i]=[]
# 	print(type(result))
# 	for idx, row in result.iterrows():
# 		subjectMatches[i].append(row["dataset"])

# subjectMatchesSecond = {}


# meta = pd.read_csv('scPFC_432_withACEandSIvariables.csv')
# for k in subjectMatches:
# 	subjectMatchesSecond[k]=[]
# 	for l in subjectMatches[k]:
# 		result = batchiden[batchiden["dataset"] == l]
# 		for idx, row in result.iterrows():
# 			subjectMatchesSecond[k].append(row["projid"])

# subjectMatchesFinal={}
# for k in subjectMatchesSecond:
# 	subjectMatchesFinal[k]=subjectMatchesSecond[k][0]

# opc.obs["patient_id"] = opc.obs["subject"].map(subjectMatchesFinal)

# op2 = opc.raw.to_adata()
# op2.var_names = op2.var['_index']
# opc.var_names= op2.var_names


# op2.var = op2.var.drop(columns="_index")

# opc.raw=op2

# opc.write("opc2.h5ad")

# adataList = []
# import scipy.sparse as sp

# for library in libraries:
#     file_path = f'/om/scratch/Mon/mabdel03/SocialIsolation/Tsai/{library}OutputP2.h5ad'
#     if os.path.exists(file_path):
#         adataL = ad.read_h5ad(file_path)
#         # QC metrics
#         if 'pct_counts_mt' not in adataL.obs:
#             adataL.var['mt'] = adataL.var_names.str.startswith('MT-')
#             adataL.var['ribo'] = adataL.var_names.str.startswith(('RPS', 'RPL'))
#             adataL.var['hb'] = adataL.var_names.str.contains('HB')
#             sc.pp.calculate_qc_metrics(adataL, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)
#         adataL.obs['batch'] = library
#         adataL.obs['patient_id'] = library
#         adataList.append(adataL)
#     # count+=1
#     # if count>50:
#     #     break

# CSVPatient = pd.read_csv('scPFC_432_withACEandSIvariables.csv')

# maleAdataList = []
# femaleAdataList = []


# adata = ad.concat(adataList, merge="same")
# adata.obs_names_make_unique()

# for i in adataList:
#     del i

# del adataList

# for i in adata.obs['patient_id'].unique():
#     k = i
#     print(i)
#     row = CSVPatient[CSVPatient["projid"] == int(k)]
#     if not row.empty:
#         msex_value = row['msex'].values[0]
#         if msex_value==1:
#             maleAdataList.append(i)
#             print("ya")
#         elif msex_value==0:
#             femaleAdataList.append(i)

# # # adataMale = adata[adata.obs['patient_id'].isin(maleAdataList)].copy()
# # # adataFemale = adata[adata.obs['patient_id'].isin(femaleAdataList)].copy()

# # # adataMale.write("maleAdata.h5ad")
# # # adataFemale.write("femaleAdata.h5ad")

# adata = adata[adata.obs['patient_id'].isin(maleAdataList + femaleAdataList)].copy()

# adata.raw = ad.AnnData(X=adata.X.copy(), var=adata.var.copy())


# sc.pp.normalize_total(adata)  # per cell
# sc.pp.log1p(adata)  # log1p

# print("norm!")

# # finding highly variable features, and also making sure m20 genes are included
# sc.pp.highly_variable_genes(adata, flavor='seurat')#,n_top_genes=2000

# print("hvg!")
# adata = adata[:, adata.var['highly_variable']]
# sc.pp.scale(adata)


# # PCA
# sc.tl.pca(adata, svd_solver='arpack')#n_comps=30
# print("pca!")

# harmony_result = hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'batch')
# adata.obsm['X_harmony'] = harmony_result.Z_corr.T 

# # Neighbors, leiden clusters
# sc.pp.neighbors(adata,use_rep='X_harmony')#n_pcs=30
# sc.tl.leiden(adata, key_added="leiden_res0_2", resolution=0.2)
# sc.tl.leiden(adata, key_added="leiden_res0_5", resolution=0.5)
# sc.tl.leiden(adata, key_added="leiden_res1", resolution=1.0)
# print("leiden!")

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


# print("umap!")

# # # adata =  ad.read_h5ad('/om/scratch/Wed/mabdel03/Subfolder/totalAdata040825.h5ad')
# # # adata.X = adata.raw.X.copy()
# # # adata2 = ad.read_h5ad('/om/scratch/Wed/mabdel03/Subfolder/totalAdataAnno040825.h5ad')
# # # adata = adata[:, adata.var['highly_variable']]
# # # adata2 = adata2[:, adata2.var['highly_variable']]
# # # adata2 = adata2[adata2.obs['cell_type'] != 'Ex_NRGN']
# # # adata = adata[adata.obs_names.isin(adata2.obs_names)].copy()

# # # adata2.write_h5ad('/om/scratch/Wed/mabdel03/Subfolder/totalAdataAnno040825.h5ad')
# # # adata.write_h5ad('/om/scratch/Wed/mabdel03/Subfolder/totalAdataAnno040825.h5ad')
# # # print(adata.obs['patient_id'].value_counts())
# # # print("Shape2")
# # # print(adata.shape)
# # # # Load the CSV file (adjust delimiter if needed)

# # # mapping_csv = pd.read_csv("cell_to_patient_assignments1.csv")

# # # mapping_csv.rename(columns={"Cell Barcode": "barcode"}, inplace=True)

# # # cells = mapping_csv.iloc[:, 0].values
# # # cells = cells.astype(str)
# # # patients = mapping_csv.iloc[:, 1].values
# # # libraries =  mapping_csv.iloc[:, 2].values
# # # adata.obs_names = adata.obs_names.str.replace(r'-1-\{object_id\}\d+$', '-1', regex=True)



# # # adata.obs['patient_id'] = np.full(len(adata.obs_names), "", dtype=str)

# # # for i in range(0,len(adata.obs_names)):
# # #   assignedPatient=""
# # #   matches = np.char.find(cells, str(adata.obs_names[i])) >= 0
# # #   possibleIndices = np.where(matches)[0]
# # #   # print(possibleIndices)
# # #   for k in possibleIndices:
# # #     if str(libraries[k]) == str(adata.obs.loc[adata.obs_names[i], "batch"]):
# # #      assignedPatient=str(patients[k])
# # #      break
# # #   adata.obs.loc[adata.obs_names[i], 'patient_id'] = assignedPatient

# # # adata.obs_names_make_unique()
# # # adata = adata[adata.obs['patient_id'] != '', :]

# # # # numbers = [
# # # #     482428, 2525608, 3430444, 7265221, 7253015, 8109170, 10100574, 10100862, 10222853, 
# # # #     10221262, 10222853, 10277308, 10516762, 11259716, 11392518, 11631558, 14452889, 
# # # #     15115927, 15844425, 15712086, 17219510, 20105242, 20124321, 20139850, 20156469, 
# # # #     20237131, 20225925, 20240514, 20267709, 20283691, 20348895, 20380417, 20594407, 
# # # #     20646778, 20865035, 20892128, 20906493, 20911508, 21145740, 21145876, 21183160, 
# # # #     21403733, 23004922, 22868024, 24747976, 35072859, 35286551, 35286551, 36492755, 
# # # #     39484737, 43074402, 50104134, 50108048, 50108750, 50109927, 50109477, 50301675, 
# # # #     51400993, 50405330, 51668135, 51520126, 60725338, 77143621, 76867532, 72650337, 
# # # #     78353027, 83984043, 84653463, 84732827, 82624422, 87410220, 86903794
# # # # ]
# # # # numbers = [
# # # #     33411712,
# # # #     50305165,
# # # #     20225925,
# # # #     20798913,
# # # #     50104846,
# # # #     20453625,
# # # #     20500815,
# # # #     10478041,
# # # #     10101291,
# # # #     91018909,
# # # #     6073025,
# # # #     10490993,
# # # #     50302978
# # # # ]

# # # # adata = adata[~adata.obs['patient_id'].isin(numbers), :]

# # # # # #filter by patient id
# # # # adata.raw = ad.AnnData(X=adata.X.copy(), var=adata.var.copy())


# # # # # Add gene names as a column if needed
# # # # # adataK.var['gene_name'] = adataK.var.index
# # df = pd.DataFrame(adata.var_names, columns=["gene"])

# # # # # # # Save as CSV
# # df.to_csv("var_names.csv", index=False)
# # # # # adata.raw.var = adata.var.copy()  # Copy gene metadata to raw
# # # # # adata.raw.obs = adata.obs.copy()  # Copy obs metadata to raw


# # # # # adataListA = []
# # # # # for library in libraries:
# # # # #     file_path = '/net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/Preprocessed_Counts/'+str(library)+'OutputP2A.h5ad'
# # # # #     if os.path.exists(file_path):
# # # # #         adataL = ad.read_h5ad(file_path)
# # # # #         if 'pct_counts_mt' not in adataL.obs:
# # # # #             adataL.var['mt'] = adataL.var_names.str.startswith('MT-')
# # # # #             adataL.var['ribo'] = adataL.var_names.str.startswith(('RPS', 'RPL'))
# # # # #             adataL.var['hb'] = adataL.var_names.str.contains('HB')  # If hemoglobin genes are relevant
# # # # #             sc.pp.calculate_qc_metrics(adataL, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)
# # # # #         adataL.obs['batch'] = library
# # # # #         adataListA.append(adataL)


# # # # # Select a single patient by batch code
# # # # # selected_patient = adata.obs['batch'].unique()[0]  # Replace with your specific patient code
# # # # # adata = adata[adata.obs['batch'] == selected_patient].copy()

# # # # adataA = ad.concat(adataListA, merge="same", index_unique="-{object_id}")#join="outer",


# # # # adataListB = []
# # # # for library in libraries:
# # # #     file_path = '/net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/Preprocessed_Counts/'+str(library)+'OutputP2B.h5ad'
# # # #     if os.path.exists(file_path):
# # # #         adataL = ad.read_h5ad(file_path)
# # # #         if 'pct_counts_mt' not in adataL.obs:
# # # #             adataL.var['mt'] = adataL.var_names.str.startswith('MT-')
# # # #             adataL.var['ribo'] = adataL.var_names.str.startswith(('RPS', 'RPL'))
# # # #             adataL.var['hb'] = adataL.var_names.str.contains('HB')  # If hemoglobin genes are relevant
# # # #             sc.pp.calculate_qc_metrics(adataL, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)
# # # #         adataL.obs['batch'] = library
# # # #         adataListB.append(adataL)


# # # # # Select a single patient by batch code
# # # # # selected_patient = adata.obs['batch'].unique()[0]  # Replace with your specific patient code
# # # # # adata = adata[adata.obs['batch'] == selected_patient].copy()

# # # # adataB = ad.concat(adataListB, merge="same", index_unique="-{object_id}")#join="outer",

# # # # adataListC = []
# # # # for library in libraries:
# # # #     file_path = '/net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/Preprocessed_Counts/'+str(library)+'OutputP2C.h5ad'
# # # #     if os.path.exists(file_path):
# # # #         adataL = ad.read_h5ad(file_path)
# # # #         if 'pct_counts_mt' not in adataL.obs:
# # # #             adataL.var['mt'] = adataL.var_names.str.startswith('MT-')
# # # #             adataL.var['ribo'] = adataL.var_names.str.startswith(('RPS', 'RPL'))
# # # #             adataL.var['hb'] = adataL.var_names.str.contains('HB')  # If hemoglobin genes are relevant
# # # #             sc.pp.calculate_qc_metrics(adataL, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)
# # # #         adataL.obs['batch'] = library
# # # #         adataListC.append(adataL)


# # # # # Select a single patient by batch code
# # # # # selected_patient = adata.obs['batch'].unique()[0]  # Replace with your specific patient code
# # # # # adata = adata[adata.obs['batch'] == selected_patient].copy()

# # # # adataC = ad.concat(adataListC, merge="same", index_unique="-{object_id}")#join="outer",

# # # # Assuming adata.X contains raw counts initially


# # # # Normalizing data
# # # # sc.pp.normalize_total(adataA)  # per cell
# # # # sc.pp.log1p(adataA)  # log1p

# # # # sc.pp.normalize_total(adataB)  # per cell
# # # # sc.pp.log1p(adataB)  # log1p

# # # # sc.pp.normalize_total(adataC)  # per cell
# # # # sc.pp.log1p(adataC)  # log1p

# # adata.raw = ad.AnnData(X=adata.X.copy(), var=adata.var.copy())

# sc.pp.normalize_total(adata)  # per cell
# sc.pp.log1p(adata)  # log1p

# # # # # Finding highly variable features
# # # # # sc.pp.highly_variable_genes(adataA, flavor='seurat')#,n_top_genes=2000
# # # # # adataA = adataA[:, adataA.var['highly_variable']]  # subsetting on this basis
# # # # # sc.pp.highly_variable_genes(adataB, flavor='seurat')#,n_top_genes=2000
# # # # # adataB = adataB[:, adataB.var['highly_variable']]  # subsetting on this basis
# # # # # sc.pp.highly_variable_genes(adataC, flavor='seurat')#,n_top_genes=2000
# # # # # adataC = adataC[:, adataC.var['highly_variable']]  # subsetting on this basis
# sc.pp.highly_variable_genes(adata, flavor='seurat')#,n_top_genes=2000

# adata = adata[:, adata.var['highly_variable']]
# # # print(adata.shape)
# # # print(adata.raw.shape)
# # # print(len(adata.var_names))

# # # # # adataK = adataK[:, adata.var['highly_variable']]  # subsetting on this basis
# # # # # adata.raw=adataK

# # # # # Scaling data
# # # # sc.pp.scale(adataA)
# # # # sc.pp.scale(adataB)
# # # # sc.pp.scale(adataC)
# sc.pp.scale(adata)

# # # # PCA
# sc.tl.pca(adata, svd_solver='arpack')#n_comps=30

# print("pcaaaa!")
# # # # print(adata.obs['patient_id'].value_counts())

# # # # print(adata.obs['patient_id'].value_counts())
# # # # print(adata.obs['patient_id'].value_counts())
# # # # print(adata.obs['patient_id'].value_counts())
# # # # #batch correction

# # # # print(adata.shape)
# # # # print(adata.obs['patient_id'].value_counts())
# # # # print(adata.obs['patient_id'].isna().sum())
# # # # print(np.isnan(adata.obsm['X_pca']).sum())
# # # # print(np.all(adata.obsm['X_pca'] == 0))



# harmony_result = hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'batch')
# adata.obsm['X_harmony'] = harmony_result.Z_corr.T 

# print("harmony!")
# # # Neighbors, leiden clusters
# sc.pp.neighbors(adata,use_rep='X_harmony')#n_pcs=30
# sc.tl.leiden(adata, key_added="leiden_res0_2", resolution=0.2)
# sc.tl.leiden(adata, key_added="leiden_res0_5", resolution=0.5)
# sc.tl.leiden(adata, key_added="leiden_res1", resolution=1.0)

# print("neighbors!")
# # # # UMAP

# sc.tl.umap(adata)
# # print("Shape3")
# # print(adata.shape)


# # # UMAP Visualization
# sc.pl.umap(adata,
#            # color=["pct_counts_mt", "log1p_total_counts", "log1p_n_genes_by_counts","pct_counts_in_top_20_genes"],
#            color=["leiden_res0_2", "leiden_res0_5", "leiden_res1"],
#            legend_loc="on data",
#            save='umapLeidenTotalPatientFinal.png')


# # UMAP Visualization
# sc.pl.umap(adata,
#            color=["pct_counts_mt", "log1p_total_counts", "log1p_n_genes_by_counts","pct_counts_in_top_20_genes"],
#            # color=["leiden_res0_2", "leiden_res0_5", "leiden_res1"],
#            legend_loc="on data",
#            save='umapLeidenTotalPatientQCFinal.png')



# # # #saving object
# adata.write('/om/scratch/Mon/mabdel03/SocialIsolation/Tsai/tsaiTotalAdata040825.h5ad')

# #NEXT PART
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




print("cell anno next!")
#creating and formatting markers DF

#remove Ex_NRGN and also non overlapping patients

# subset = ['83984043', '77143621', '72777797', '76867532', '51624179', '50100806', '50197261', '50402729', '72650337', '50302554', '53772202', '6107196', '89546375', '61344957', '50302266', '61142759', '24680888', '50107619', '94430339', '68015667', '50301099', '75990666', '50500136', '96095092', '41285665', '78353027', '67429065', '89903942', '64291829', '32383679', '14498577', '36830117', '50301125', '50106442', '50301387', '62404688', '67185070', '66754397', '85065193', '2899847', '43485807', '83173570', '74753465', '7265221', '36492755', '22789958', '71626373', '50301675', '90942860', '31726180', '23004922', '29933130', '23690880', '17929065', '50109639', '35941263', '65652206', '65001607', '66406040', '52311825', '82317494', '44671043', '75169847', '82684406', '85171938', '92629514', '20907534', '50500550', '54324848', '14184286', '50109477', '50400835', '77891596', '87645445', '70883578', '32705437', '68539908', '63188799', '15218541', '68525196', '50105699', '94144536', '50108912', '50404299', '60848460', '24747976', '50400259', '89001223', '70625336', '39989287', '80790863', '37436329', '22868024', '64493027', '18659212', '61029627', '50107745', '44299049', '50106730', '50102114', '92393245', '50402693', '69982533', '50101523', '37178462', '9841821', '8132197', '50107907', '45115248', '9650662', '18920002', '46251007', '50300084', '50300408', '50103129', '50406057', '50405330', '8109170', '65736039', '91444029', '31509843', '83034844', '66924745', '46000440', '3713990', '50100932', '37527863', '77239958', '48023497', '50106280', '94974890', '50401002', '60725338', '22776575', '43074402', '46246604', '7311370', '39484737', '44019405', '7253015', '98953007', '93815598', '70816595', '37065652', '34962204', '50410319', '51826377', '74284255', '87410220', '61827429', '3889845', '89164957', '24039289', '91018909', '50102376', '94092977', '50402567', '52940446', '84417209', '76733461', '16068769', '33411712', '50300822', '90267190', '65292866', '70917649', '91707643', '21362537', '51520126', '37030589', '50107583', '50400385', '60747316', '70153803', '6804844', '51815338', '41773404', '33501827', '81086436', '18414513', '50109927', '45566083', '34726040', '87779516', '57978756', '90780976', '26631069', '50105301', '84642424', '86934089', '78452313', '95491648', '68745332', '50403446']

pandas2ri.activate() # Load the RDS file
rds_data = r['readRDS']('Brain_Human_PFC_Markers_Mohammadi2020.rds')
markers_df = pandas2ri.rpy2py(rds_data)
cell_type_names = list(markers_df.names)

cell_types = list(markers_df)
# load in h5ad
adata = ad.read_h5ad('/om/scratch/Mon/mabdel03/SocialIsolation/Tsai/tsaiTotalAdata040825.h5ad')
adata2 = ad.read_h5ad('/om/scratch/Mon/mabdel03/SocialIsolation/Tsai/tsaiTotalAdataAnno040825.h5ad')
adata = adata[adata.obs_names.isin(adata2.obs_names)]
data = []

for cell_type_name, genes in zip(cell_type_names, cell_types):
    for gene in genes:
        data.append({'source': cell_type_name, 'target': gene, 'weight': 1.0})

markers_df = pd.DataFrame(data)


# adata.raw = adata
if 'ora_estimate' in adata.obsm:
    del adata.obsm['ora_estimate']

 #Running ORA/whole process there

dc.run_ora(adata,markers_df,source = "source",target = "target", use_raw=False)

acts = dc.get_acts(adata, obsm_key='ora_estimate')

# all_marker_genes = pd.unique(markers_df['target'])
# X = adata.raw.X if adata.raw is not None else adata.X

# marker_gene_expr = pd.DataFrame(
#     X[:, adata.var_names.isin(all_marker_genes)].toarray(),  # convert sparse to dense if needed
#     index=adata.obs_names,
#     columns=adata.var_names[adata.var_names.isin(all_marker_genes)]
# )


# marker_gene_expr.to_csv("marker_gene_expressionDejager.csv")

# We need to remove inf and set them to the maximum value observed for pvals=0
acts_v = acts.X.ravel()
max_e = np.nanmax(acts_v[np.isfinite(acts_v)])
acts.X[~np.isfinite(acts.X)] = max_e

# acts_df = pd.DataFrame(acts.X, index=acts.obs_names, columns=acts.var_names)

# acts_df.to_csv("enrichment_scoresDejager.csv")

sc.pl.umap(acts, color=['Ast', 'leiden_res0_5'], cmap='RdBu_r')
sc.pl.violin(acts, keys=['Ast'], groupby='leiden_res0_5')

df = dc.rank_sources_groups(acts, groupby='leiden_res0_5', reference='rest', method='t-test_overestim_var')

n_ctypes = 3
ctypes_dict = df.groupby('group').head(n_ctypes).groupby('group')['names'].apply(lambda x: list(x)).to_dict()
# with open("ctypes_dictDejager.txt", "w") as f:
#     for key, value in ctypes_dict.items():
#         f.write(f"{key}: {value}\n")

sc.pl.matrixplot(acts, ctypes_dict, 'leiden_res0_5', dendrogram=True, standard_scale='var',
                 colorbar_title='Z-scaled scores', cmap='RdBu_r')


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

adata.write('/om/scratch/Mon/mabdel03/SocialIsolation/Tsai/tsaiTotalAdataAnno040825.h5ad')


adK = adata.raw.to_adata()
adK = adK[adK.obs["cell_type"].notna()].copy()
adK = adK[~adK.obs["cell_type"].isin(["Ex_NRGN"])].copy()
adK.obs["patient_id"] = adK.obs["patient_id"].astype(str)


adata = adata[adata.obs["cell_type"].notna()].copy()
adata = adata[~adata.obs["cell_type"].isin(["Ex_NRGN"])].copy()
adata.obs["patient_id"] = adata.obs["patient_id"].astype(str)
import numpy as np
import pandas as pd
from scipy import sparse

from collections import defaultdict

adata.raw= adK
adata.write("tsaiTotalAdataAnno040825.h5ad")

adNewExc = adata[adata.obs["cell_type"].isin(["Ex_L4_5", "Ex_L5_6", "Ex_L5_6_CC", "Ex_L2_3", "Ex_L4","Ex_L5"])].copy()
adNewInh = adata[adata.obs["cell_type"].isin(["In_VIP", "In_PV (Basket)", "In_PV (Chandelier)", "In_SST", "In_Rosehip"])].copy()
adNewOli = adata[adata.obs["cell_type"].isin(["Oli"])].copy()
adNewOpc = adata[adata.obs["cell_type"].isin(["OPC"])].copy()
adNewAst = adata[adata.obs["cell_type"].isin(["Ast"])].copy()
adNewMic = adata[adata.obs["cell_type"].isin(["Mic"])].copy()

adNewExc.write("/om/scratch/Mon/mabdel03/SocialIsolation/Tsai/excAnno.h5ad")
adNewInh.write("/om/scratch/Mon/mabdel03/SocialIsolation/Tsai/inhAnno.h5ad")
adNewOli.write("/om/scratch/Mon/mabdel03/SocialIsolation/Tsai/oliAnno.h5ad")
adNewOpc.write("/om/scratch/Mon/mabdel03/SocialIsolation/Tsai/opcAnno.h5ad")
adNewAst.write("/om/scratch/Mon/mabdel03/SocialIsolation/Tsai/astAnno.h5ad")
adNewMic.write("/om/scratch/Mon/mabdel03/SocialIsolation/Tsai/micAnno.h5ad")
# adata1 = adata[adata.obs['leiden_res0_5']!=12].copy()
# adata1= adata1[adata1.obs['leiden_res0_5']!="12"].copy()

# adata = ad.read_h5ad('/om/scratch/Mon/mabdel03/SocialIsolation/Tsai/tsaiTotalAdata040825.h5ad')

# adataadata.obs[]

# data = []

# cell_type_names = [name for name in cell_type_names if name.lower() != remove_name.lower()]
# cell_types = [
#     ct for name, ct in zip(cell_type_names, cell_types)
#     if name.lower() != remove_name.lower()
# ]

# for cell_type_name, genes in zip(cell_type_names, cell_types):
#     for gene in genes:
#         data.append({'source': cell_type_name, 'target': gene, 'weight': 1.0})

# markers_df = pd.DataFrame(data)


# # adata.raw = adata
# if 'ora_estimate' in adata.obsm:
#     del adata.obsm['ora_estimate']

#  #Running ORA/whole process there

# dc.run_ora(adata,markers_df,source = "source",target = "target", use_raw=False)

# acts = dc.get_acts(adata, obsm_key='ora_estimate')

# # all_marker_genes = pd.unique(markers_df['target'])
# # X = adata.raw.X if adata.raw is not None else adata.X

# # marker_gene_expr = pd.DataFrame(
# #     X[:, adata.var_names.isin(all_marker_genes)].toarray(),  # convert sparse to dense if needed
# #     index=adata.obs_names,
# #     columns=adata.var_names[adata.var_names.isin(all_marker_genes)]
# # )


# # marker_gene_expr.to_csv("marker_gene_expressionDejager.csv")

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


# adata.write('/om/scratch/Mon/mabdel03/SocialIsolation/Tsai/tsaiTotalAdataAnno040825.h5ad')

import anndata as ad
import scanpy as sc
import os
import numpy as np
import matplotlib.pyplot as plt
sc.pl.umap(adata,
           # color=["pct_counts_mt", "log1p_total_counts", "log1p_n_genes_by_counts","pct_counts_in_top_20_genes"],
           color=["cell_type"],
           legend_loc="on data",
           save='umapLeidenTotalPatientCellTypeTsai.png')


# # # source + target adata - for batch info transfer

# import anndata as ad
# import scanpy as sc
# import os
# import numpy as np
# import matplotlib.pyplot as plt

# file_path_target = '/om/scratch/Wed/mabdel03/Subfolder/totalAdataAnno040825.h5ad'
# adata = ad.read_h5ad(file_path_target)
# print(adata)

# patient_ids = [21000180, 20252720, 20933324, 21293107, 98953007, 20358955, 21000630, 81852640,20152393, 50104846, 20254902, 54122640, 20195344, 81874628, 20706215, 20120255,20232474, 89614402, 60725338, 32383679, 20983400, 21157370, 20153984, 20236838,21411459, 63874408, 50106578, 50108886, 20221273, 20254588, 21140119, 61827429,75675336, 69432088]
# patient_ids=[str(i) for i in patient_ids]
# print(adata.obs['patient_id'].unique())
# adata=adata[adata.obs['patient_id'].isin(patient_ids)]



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