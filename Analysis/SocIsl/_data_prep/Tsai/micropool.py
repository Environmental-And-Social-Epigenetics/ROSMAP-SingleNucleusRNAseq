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

males= ['2899847', '3713990', '6107196', '6804844', '7311370', '9841821', '15218541', '16068769', '18414513', '18659212', '21362537', '22789958', '22868024', '23004922', '24039289', '24680888', '31509843', '32705437', '33411712', '33501827', '34726040', '34962204', '36830117', '37178462', '37436329', '37527863', '39484737', '39989287', '43074402', '43485807', '44019405', '44299049', '44671043', '45566083', '46000440', '50101523', '50102114', '50102376', '50105301', '50106730', '50107583', '50107745', '50108912', '50109477', '50109927', '50301099', '50302266', '50302554', '50400259', '50400835', '50402693', '50405330', '50406057', '50410319', '51520126', '51815338', '52311825', '62404688', '64291829', '65292866', '65652206', '65736039', '66406040', '66754397', '67429065', '68015667', '68525196', '68539908', '68745332', '70816595', '70883578', '71626373', '72777797', '74284255', '75990666', '77143621', '77239958', '77891596', '78353027', '78452313', '82317494', '83984043', '89001223', '89546375', '89903942', '90267190', '90942860', '91018909', '91444029', '91707643', '92629514', '94092977', '95491648', '86408244']

females= ['3889845', '7253015', '7265221', '8109170', '8132197', '9650662', '14184286', '14498577', '17929065', '18920002', '20907534', '22776575', '23690880', '24747976', '26631069', '29933130', '31726180', '32383679', '35941263', '36492755', '37030589', '37065652', '41285665', '41773404', '45115248', '46246604', '46251007', '48023497', '50100806', '50100932', '50103129', '50105699', '50106280', '50106442', '50107619', '50107907', '50109639', '50197261', '50300084', '50300408', '50300822', '50301125', '50301387', '50301675', '50400385', '50401002', '50402567', '50402729', '50403446', '50404299', '50500136', '50500550', '51624179', '51826377', '52940446', '53772202', '54324848', '57978756', '60725338', '60747316', '60848460', '61029627', '61142759', '61344957', '61827429', '63188799', '64493027', '65001607', '66924745', '67185070', '69982533', '70153803', '70625336', '70917649', '72650337', '74753465', '75169847', '76733461', '76867532', '80790863', '81086436', '82684406', '83034844', '83173570', '84417209', '84642424', '85065193', '85171938', '86934089', '87410220', '87645445', '87779516', '89164957', '90780976', '92393245', '93815598', '94144536', '94430339', '94974890', '96095092', '98953007', '50102088']


variables = ['Exc','Inh','Oli','Ast','OPC','Mic']
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
		rootDirectory = "/om/scratch/Mon/mabdel03/SocialIsolation/Tsai/" 
		os.chdir(rootDirectory)
		folderPath="/om/scratch/Mon/mabdel03/SocialIsolation/Tsai/" #change path to folder
		db_paths = [folderPath+"hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
		folderPath+"hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"]
		dbs = [RankingDatabase(fname=p, name=os.path.basename(p)) for p in db_paths]
		f_tfs = folderPath+"hg.txt"
		tf_names = load_tf_names(f_tfs)
		CSVPatient = pd.read_csv('dataset_652_basic_03-23-2022.csv')
		adata = ad.read_h5ad(x.lower()+"2.h5ad")
		adata.obs["patient_id"]=adata.obs["patient_id"].astype("str")
		adata = adata[adata.obs['patient_id'].isin(males)]
		print(adata)
		print("read it in")
		metaDF = adata.obs.copy()
		metaDF['patient_id'] = metaDF['patient_id'].astype(str)
		CSVPatient['projid'] = CSVPatient['projid'].astype(str)
		metaDF = metaDF.merge(
		CSVPatient[['projid', 'social_isolation_avg','msex','age_death','pmi','niareagansc']],how='left',left_on='patient_id',right_on='projid')
		metaDF.drop(columns=['projid'], inplace=True)
		adata = adata[metaDF["msex"] == 1].copy()
		print("meta")
		adata = fixed_micropool(adata, patient_col="patient_id", pool_size=50)
		print("micropool")
		metaDF = adata.obs.copy()
		metaDF['patient_id'] = metaDF['patient_id'].astype(str)
		CSVPatient['projid'] = CSVPatient['projid'].astype(str)
		metaDF = metaDF.reset_index().merge(
		CSVPatient[['projid', 'social_isolation_avg','msex','age_death','pmi','niareagansc']],how='left',left_on='patient_id',right_on='projid')
		metaDF.drop(columns=['projid'], inplace=True)
		metaDF['age_death'] = ((metaDF['age_death'] - metaDF['age_death'].mean()) / metaDF['age_death'].std())
		metaDF['pmi'] = ((metaDF['pmi'] - metaDF['pmi'].mean()) / metaDF['pmi'].std())
		metaDF = metaDF.set_index('pool_id')
		metaDF.to_csv("metaDFMale"+x+".csv")
		print("pool id")
		# if adata.shape[0] > 0:
		# 	X_sparse = adata.X if sparse.issparse(adata.X) else sparse.csr_matrix(adata.X)
		# 	X_sparse = X_sparse.T
		# 	col_sums = X_sparse.sum(axis=0).A1 if hasattr(X_sparse, "A1") else np.array(X_sparse.sum(axis=0)).ravel()
		# 	col_sums[col_sums == 0] = 1
		# 	scaling_factors = 1e6 / col_sums
		# 	D = sparse.diags(scaling_factors)
		# 	X_sparseNew = X_sparse.dot(D)
		# 	df1 = pd.DataFrame.sparse.from_spmatrix(X_sparseNew,index=adata.var_names,columns=adata.obs_names)
		# 	df1.to_csv("carlesmale_matrix"+x+".tsv", sep="\t")

for x in variables:
	if __name__ == '__main__':
		#todo: put files into folder
		cellTypeToSubset =x
		rootDirectory = "/om/scratch/Mon/mabdel03/SocialIsolation/Tsai/" 
		os.chdir(rootDirectory)
		folderPath="/om/scratch/Mon/mabdel03/SocialIsolation/Tsai/" #change path to folder
		db_paths = [folderPath+"hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
		folderPath+"hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"]
		dbs = [RankingDatabase(fname=p, name=os.path.basename(p)) for p in db_paths]
		f_tfs = folderPath+"hg.txt"
		tf_names = load_tf_names(f_tfs)
		CSVPatient = pd.read_csv('dataset_652_basic_03-23-2022.csv')
		adata = ad.read_h5ad(x.lower()+"2.h5ad")
		adata.obs["patient_id"]=adata.obs["patient_id"].astype("str")
		adata = adata[adata.obs['patient_id'].isin(females)]
		print(adata)
		print("read it in")
		metaDF = adata.obs.copy()
		metaDF['patient_id'] = metaDF['patient_id'].astype(str)
		CSVPatient['projid'] = CSVPatient['projid'].astype(str)
		metaDF = metaDF.merge(
		CSVPatient[['projid', 'social_isolation_avg','msex','age_death','pmi','niareagansc']],how='left',left_on='patient_id',right_on='projid')
		metaDF.drop(columns=['projid'], inplace=True)
		adata = adata[metaDF["msex"] == 0].copy()
		print("meta")
		adata = fixed_micropool(adata, patient_col="patient_id", pool_size=50)
		print("micropool")
		metaDF = adata.obs.copy()
		metaDF['patient_id'] = metaDF['patient_id'].astype(str)
		CSVPatient['projid'] = CSVPatient['projid'].astype(str)
		metaDF = metaDF.reset_index().merge(
		CSVPatient[['projid', 'social_isolation_avg','msex','age_death','pmi','niareagansc']],how='left',left_on='patient_id',right_on='projid')
		metaDF.drop(columns=['projid'], inplace=True)
		metaDF['age_death'] = ((metaDF['age_death'] - metaDF['age_death'].mean()) / metaDF['age_death'].std())
		metaDF['pmi'] = ((metaDF['pmi'] - metaDF['pmi'].mean()) / metaDF['pmi'].std())
		metaDF = metaDF.set_index('pool_id')
		metaDF.to_csv("metaDFFemale"+x+".csv")
		print("pool id")
		# if adata.shape[0] > 0:
		# 	X_sparse = adata.X if sparse.issparse(adata.X) else sparse.csr_matrix(adata.X)
		# 	X_sparse = X_sparse.T
		# 	col_sums = X_sparse.sum(axis=0).A1 if hasattr(X_sparse, "A1") else np.array(X_sparse.sum(axis=0)).ravel()
		# 	col_sums[col_sums == 0] = 1
		# 	scaling_factors = 1e6 / col_sums
		# 	D = sparse.diags(scaling_factors)
		# 	X_sparseNew = X_sparse.dot(D)
		# 	df1 = pd.DataFrame.sparse.from_spmatrix(X_sparseNew,index=adata.var_names,columns=adata.obs_names)
		# 	df1.to_csv("carlesfemale_matrix"+x+".tsv", sep="\t")