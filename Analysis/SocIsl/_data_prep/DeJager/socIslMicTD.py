# DEPRECATED: Legacy data preparation script from Openmind era.
# The current Processing pipeline (Processing/*/Pipeline/) supersedes this.
# Paths may still reference decommissioned Openmind directories.
#
import anndata as ad
import pandas as pd
# variables = ['In_PV (Basket)', 'In_PV (Chandelier)']
#, 'Ex_L4_5''Ex_L5/6-CC', 'Ex_NRGN','Ast', 'Endo', 'Ex_L2_3', 'Ex_L4', 'Ex_L5', 'Ex_L5_6', 'In_Rosehip', 'In_SST', 'In_VIP', 'Mic', 'Oli', 'OPC'
from scipy.stats import ranksums
from scipy import sparse
import rdata
import scanpy as sc
import numpy as np
import os
import glob
import pickle
import numpy as np
import loompy as lp
import re
import asyncio
from dask.diagnostics import ProgressBar
from dask.distributed import LocalCluster, Client
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2
from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_rss
import seaborn as sns
import matplotlib.pyplot as plt
# print(adataOG.shape)
from statsmodels.stats.multitest import multipletests

import statsmodels.formula.api as smf
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

file_path_target = 'totalAdataAnno3D.h5ad'

adataOG = ad.read_h5ad(file_path_target)
variables = ['Ast']
def fixed_micropool(adata, patient_col="patient_id", pool_size=50):
	pooled_X = []
	pooled_obs = []
	obs_cols = [c for c in adata.obs.columns if c != patient_col]
	adata = ad.AnnData(X = adata.raw.to_adata().X,obs = adata.obs.copy(),var = adata.raw.to_adata().var.copy())
	for patient in adata.obs[patient_col].unique():
		idx = np.where(adata.obs[patient_col] == patient)[0]
		n_cells = len(idx)
		n_pools = n_cells // pool_size
		for i in range(n_pools):
			pool_idx = idx[i*pool_size : (i+1)*pool_size]
			pooled_counts = adata[pool_idx].X.sum(axis=0)
			if hasattr(pooled_counts, "toarray"):  # if sparse
				pooled_counts = pooled_counts.toarray().ravel()
			else:
				pooled_counts = np.array(pooled_counts).ravel()
			pooled_X.append(pooled_counts)
			pooled_info = {
				"patient_id": patient,
				"pool_id": f"{patient}_pool{i}"
			}
			obs_subset = adata.obs.iloc[pool_idx]
			for col in obs_cols:
				if pd.api.types.is_numeric_dtype(obs_subset[col]):
					pooled_info[col] = obs_subset[col].mean()
				else:
					pooled_info[col] = obs_subset[col].mode().iloc[0]
			pooled_obs.append(pooled_info)
	pooled_X = np.vstack(pooled_X)
	pooled_obs = pd.DataFrame(pooled_obs)
	adata_pooled = sc.AnnData(X=pooled_X)
	adata_pooled.var_names = adata.var_names
	adata_pooled.obs = pooled_obs.set_index("pool_id")
	return adata_pooled

for x in variables:
	if __name__ == '__main__':
		#todo: put files into folder
		cellTypeToSubset =x
		#defining paths to files
		rootDirectory = "/om/scratch/Mon/mabdel03/SocialIsolation/" 
		os.chdir(rootDirectory)
		folderPath="/om/scratch/Mon/mabdel03/SocialIsolation/" #change path to folder
		db_paths = [folderPath+"hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
		folderPath+"hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"]
		dbs = [RankingDatabase(fname=p, name=os.path.basename(p)) for p in db_paths]
		f_tfs = folderPath+"hg.txt"
		tf_names = load_tf_names(f_tfs)
		CSVPatient = pd.read_csv('dataset_652_basic_03-23-2022.csv')
		adata_filtered=adataOG
		if cellTypeToSubset=="Exc":
			adata_filtered = adataOG[adataOG.obs['cell_type'].isin(["Ex_L4_5", "Ex_L5_6", "Ex_L5_6_CC", "Ex_L2_3", "Ex_L4","Ex_L5"])]
		elif cellTypeToSubset=="Inh":
			adata_filtered = adataOG[adataOG.obs['cell_type'].isin(["In_VIP", "In_PV (Basket)", "In_PV (Chandelier)", "In_SST", "In_Rosehip"])]
		else:
			adata_filtered = adataOG[adataOG.obs['cell_type']==cellTypeToSubset]
		adata = adata_filtered.to_memory()
		metaDF=adata.obs
		metaDF['patient_id'] = metaDF['patient_id'].astype(str)
		CSVPatient['projid'] = CSVPatient['projid'].astype(str)
		metaDF = metaDF.merge(
		CSVPatient[['projid', 'social_isolation_avg','msex']],how='left',left_on='patient_id',right_on='projid')
		metaDF.drop(columns=['projid'], inplace=True)
		adata = adata[metaDF["msex"] == 1] #subset for sex
		print(adata)
		adata = fixed_micropool(adata, patient_col="patient_id", pool_size=50) #micropooling
		metaDF=adata.obs #regenerate metadata
		metaDF['patient_id'] = metaDF['patient_id'].astype(str)
		CSVPatient['projid'] = CSVPatient['projid'].astype(str)
		metaDF = metaDF.reset_index().merge(
		CSVPatient[['projid', 'social_isolation_avg','msex']],how='left',left_on='patient_id',right_on='projid')
		metaDF.drop(columns=['projid'], inplace=True)
		metaDF = metaDF.set_index('pool_id')
		metaDF.to_csv("male_metaDF"+x+".csv")
		if adata.shape[0] > 0:
			X_sparse = adata.X if sparse.issparse(adata.X) else sparse.csr_matrix(adata.X)
			print(adata)
			X_sparse=X_sparse.T
			col_sums = X_sparse.sum(axis=0).A1 if hasattr(X_sparse, "A1") else np.array(X_sparse.sum(axis=0)).ravel()
			col_sums[col_sums == 0] = 1
			scaling_factors = 1e6 / col_sums
			D = sparse.diags(scaling_factors)
			X_sparseNew = X_sparse.dot(D)
			df1 = pd.DataFrame.sparse.from_spmatrix(X_sparseNew,index=adata.var_names,columns=adata.obs_names)
			df1.to_csv("dejagermale_matrix"+x+".tsv", sep="\t")


for x in variables:
	if __name__ == '__main__':
		#todo: put files into folder
		cellTypeToSubset =x
		#defining paths to files
		rootDirectory = "/om/scratch/Mon/mabdel03/SocialIsolation/" 
		os.chdir(rootDirectory)
		folderPath="/om/scratch/Mon/mabdel03/SocialIsolation/" #change path to folder
		db_paths = [folderPath+"hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
		folderPath+"hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"]
		dbs = [RankingDatabase(fname=p, name=os.path.basename(p)) for p in db_paths]
		f_tfs = folderPath+"hg.txt"
		tf_names = load_tf_names(f_tfs)
		CSVPatient = pd.read_csv('dataset_652_basic_03-23-2022.csv')
		adata_filtered=adataOG
		if cellTypeToSubset=="Exc":
			adata_filtered = adataOG[adataOG.obs['cell_type'].isin(["Ex_L4_5", "Ex_L5_6", "Ex_L5_6_CC", "Ex_L2_3", "Ex_L4","Ex_L5"])]
		elif cellTypeToSubset=="Inh":
			adata_filtered = adataOG[adataOG.obs['cell_type'].isin(["In_VIP", "In_PV (Basket)", "In_PV (Chandelier)", "In_SST", "In_Rosehip"])]
		else:
			adata_filtered = adataOG[adataOG.obs['cell_type']==cellTypeToSubset]
		adata = adata_filtered.to_memory()
		metaDF=adata.obs
		metaDF['patient_id'] = metaDF['patient_id'].astype(str)
		CSVPatient['projid'] = CSVPatient['projid'].astype(str)
		metaDF = metaDF.merge(
		CSVPatient[['projid', 'social_isolation_avg','msex']],how='left',left_on='patient_id',right_on='projid')
		metaDF.drop(columns=['projid'], inplace=True)
		adata = adata[metaDF["msex"] == 0] #subset for sex
		print(adata)
		adata = fixed_micropool(adata, patient_col="patient_id", pool_size=50) #micropooling
		metaDF=adata.obs #regenerate metadata
		metaDF['patient_id'] = metaDF['patient_id'].astype(str)
		CSVPatient['projid'] = CSVPatient['projid'].astype(str)
		metaDF = metaDF.reset_index().merge(
		CSVPatient[['projid', 'social_isolation_avg','msex']],how='left',left_on='patient_id',right_on='projid')
		metaDF.drop(columns=['projid'], inplace=True)
		metaDF = metaDF.set_index('pool_id')
		metaDF.to_csv("female_metaDF"+x+".csv")
		if adata.shape[0] > 0:
			X_sparse = adata.X if sparse.issparse(adata.X) else sparse.csr_matrix(adata.X)
			print(adata)
			X_sparse=X_sparse.T
			col_sums = X_sparse.sum(axis=0).A1 if hasattr(X_sparse, "A1") else np.array(X_sparse.sum(axis=0)).ravel()
			col_sums[col_sums == 0] = 1
			scaling_factors = 1e6 / col_sums
			D = sparse.diags(scaling_factors)
			X_sparseNew = X_sparse.dot(D)
			df1 = pd.DataFrame.sparse.from_spmatrix(X_sparseNew,index=adata.var_names,columns=adata.obs_names)
			df1.to_csv("dejagerfemale_matrix"+x+".tsv", sep="\t")