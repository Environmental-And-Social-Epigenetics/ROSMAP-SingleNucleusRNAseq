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

variables=["OPC","Exc","Inh","Ast","Mic","Oli"]


def fixed_micropool(adata, patient_col="patient_id", pool_size=30):
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
		rootDirectory = os.path.join(os.environ["SOCISL_OUTPUT_ROOT"], "Tsai", "") 
		os.chdir(rootDirectory)
		folderPath=os.path.join(os.environ["SOCISL_OUTPUT_ROOT"], "Tsai", "") #change path to folder
		db_paths = [folderPath+"hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
		folderPath+"hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"]
		dbs = [RankingDatabase(fname=p, name=os.path.basename(p)) for p in db_paths]
		f_tfs = folderPath+"hg.txt"
		tf_names = load_tf_names(f_tfs)
		CSVPatient = pd.read_csv('dataset_652_basic_03-23-2022.csv')
		adata_filtered=ad.AnnData()
		if cellTypeToSubset=="Exc":
			adata_filtered = ad.read_h5ad("excAnno.h5ad")
		elif cellTypeToSubset=="Inh":
			adata_filtered = ad.read_h5ad("inhAnno.h5ad")
		elif cellTypeToSubset=="Mic":
			adata_filtered = ad.read_h5ad("micAnno.h5ad")
		elif cellTypeToSubset=="Oli":
			adata_filtered = ad.read_h5ad("oliAnno.h5ad")
		elif cellTypeToSubset=="OPC":
			adata_filtered = ad.read_h5ad("opcAnno.h5ad")
		elif cellTypeToSubset=="Ast":
			adata_filtered = ad.read_h5ad("astAnno.h5ad")
		adata = adata_filtered.copy()
		adata.raw = adata_filtered.raw.to_adata()
		metaDF=adata.obs
		metaDF['batch'] = metaDF['batch'].astype(str)
		CSVPatient['projid'] = CSVPatient['projid'].astype(str)
		metaDF = metaDF.merge(
		CSVPatient[['projid', 'social_isolation_avg','msex']],how='left',left_on='batch',right_on='projid')
		metaDF.drop(columns=['projid'], inplace=True)
		adata = adata[metaDF["msex"] == 1].copy()
		adata = fixed_micropool(adata, patient_col="patient_id", pool_size=100)
		metaDF=adata.obs
		metaDF['batch'] = metaDF['batch'].astype(str)
		CSVPatient['projid'] = CSVPatient['projid'].astype(str)
		metaDF = metaDF.reset_index().merge(
		CSVPatient[['projid', 'social_isolation_avg','msex']],how='left',left_on='batch',right_on='projid')
		metaDF.drop(columns=['projid'], inplace=True)
		metaDF = metaDF.set_index('pool_id')
		if adata.shape[0] > 0:
			X_sparse = adata.X if sparse.issparse(adata.X) else sparse.csr_matrix(adata.X)
			X_sparse = X_sparse.T
			col_sums = X_sparse.sum(axis=0).A1 if hasattr(X_sparse, "A1") else np.array(X_sparse.sum(axis=0)).ravel()
			col_sums[col_sums == 0] = 1
			scaling_factors = 1e6 / col_sums
			D = sparse.diags(scaling_factors)
			X_sparseNew = X_sparse.dot(D)
			df1 = pd.DataFrame.sparse.from_spmatrix(X_sparseNew,index=adata.var_names,columns=adata.obs_names)
			df1.to_csv("Nmale_matrix"+x+".tsv", sep="\t")
			X_sparse = X_sparse.T
			df = pd.DataFrame.sparse.from_spmatrix(X_sparse,index=adata.obs_names,columns=adata.var_names)
			# df = pd.DataFrame(X_sparse,index=adata.obs_names, columns=adata.var_names)
			adjacencies = grnboost2(expression_data=df, tf_names=tf_names)
			modules= list(modules_from_adjacencies(adjacencies, df))
			enrichedMot = prune2df(dbs, modules, folderPath+"motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl", num_workers=4)  # <-- safer to set single-threaded
			regulons = df2regulons(enrichedMot) 
			with open(rootDirectory+"male__regulons_Tsai"+x+".p", "wb") as f:
				pickle.dump(regulons, f)
			dfReg = pd.DataFrame (regulons, columns = ['regulon'])
			#List of enriched regulons
			dfReg.to_csv("male_regulonsFULL_Tsai"+x+".csv")
			# AUCell
			auc_mtx = aucell(df, regulons, num_workers=4)
			auc_mtx.to_csv(rootDirectory+"male_auc_mtxFULL_Tsai"+x+".csv");
			#Calculating RSS: score forchanges in expression of different regulons across groups/strength of regulon importance in these conditions
			# map=sns.clustermap(auc_mtx, figsize=(8,8))
			# plt.savefig(rootDirectory+'fem_figMtxFULL'+x+'.png', dpi=400)
			# print(rss_cellType.index)
			auc_scores = auc_mtx
			results = []
			for regulon in auc_scores.columns:
			# if True:
				if True:
					# regulon="SREBF2(+)"
					print("hi")
					metadata=metaDF
					print(auc_mtx.head)
					print(auc_mtx[regulon])
					data = pd.DataFrame({'AUCell': auc_mtx[regulon],'Condition': metadata['social_isolation_avg'],'Batch': metadata['batch']})
					print(data.head)
					plt.figure(figsize=(6,4))
					sns.regplot(x="Condition",y="AUCell",data=data,scatter_kws={"s": 10, "alpha": 0.4},line_kws={"color": "red"})
					plt.title(x+" AUCell vs Condition")
					plt.xlabel("Condition (continuous)")
					plt.ylabel("AUCell score")
					plt.tight_layout()
					plt.savefig(x+"_male_AUCell_vs_Condition.png", dpi=300)
					try:
						print("hii")
						data["Condition"] = pd.to_numeric(data["Condition"], errors="coerce")
						model = smf.ols("AUCell ~ Condition", data=data)
						fit = model.fit()
						term = next((k for k in fit.params.keys() if "Condition" in k), None)
						coef = fit.params[term]
						pval0 = fit.pvalues[term]
						print("summary:"+str(regulon))
						print(fit.summary())
						print(regulon)
						results.append({"regulon": regulon,"coef": coef,"pval": pval0})
						print(regulon + ":" + str(pval))
					except Exception as e:
						results.append({"error": str(e)})
				results_df = pd.DataFrame(results).dropna(subset=["coef", "pval"])
				_, pvals_fdr, _, _ = multipletests(results_df["pval"], alpha=0.05, method='fdr_bh')
				results_df["fdr"] = pvals_fdr
				results_df["log_fdr"] = -np.log10(results_df["fdr"])
				results_df["significance"] = "Not Significant"
				results_df.loc[(results_df["fdr"] < 0.05) & (results_df["coef"] > 0), "significance"] = "Up"
				results_df.loc[(results_df["fdr"] < 0.05) & (results_df["coef"] < 0), "significance"] = "Down"
				results_df.to_csv("male_pVals"+x+".csv")
				palette = {"Up": "red","Down": "blue","Not Significant": "gray"}
				plt.figure(figsize=(10, 6))
				sns.scatterplot(
					data=results_df, x="coef", y="log_fdr", hue="significance",
					palette=palette, edgecolor=None
				)
				for _, row in results_df[(results_df["significance"] == "Up") | (results_df["significance"] == "Down")].iterrows():
					plt.text(
						row["coef"], row["log_fdr"], row["regulon"],
						fontsize=8, ha='right', va='bottom'
					)
				plt.axhline(-np.log10(0.05), linestyle="--", color="black", linewidth=1)
				plt.axvline(0, linestyle="--", color="black", linewidth=1)
				plt.xlabel("effect size")
				plt.ylabel("-log10(FDR)")
				plt.title(f"{x} Social Isolation Transcription Factor Activity - Volcano")
				plt.legend(title="Significance", loc="upper right")
				plt.tight_layout()
				plt.savefig(x + "TsaiMaleVolcano.png") #volcano plot!
				plt.show()


for x in variables:
	if __name__ == '__main__':
		# variables=['Ex_L2_3']
		#todo: put files into folder
		cellTypeToSubset =x
		#defining paths to files
		rootDirectory = os.path.join(os.environ["SOCISL_OUTPUT_ROOT"], "Tsai", "") 
		os.chdir(rootDirectory)
		folderPath=os.path.join(os.environ["SOCISL_OUTPUT_ROOT"], "Tsai", "") #change path to folder
		db_paths = [folderPath+"hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
		folderPath+"hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"]
		dbs = [RankingDatabase(fname=p, name=os.path.basename(p)) for p in db_paths]
		f_tfs = folderPath+"hg.txt"
		tf_names = load_tf_names(f_tfs)
		CSVPatient = pd.read_csv('dataset_652_basic_03-23-2022.csv')
		adata_filtered=ad.AnnData()
		if cellTypeToSubset=="Exc":
			adata_filtered = ad.read_h5ad("excAnno.h5ad")
		elif cellTypeToSubset=="Inh":
			adata_filtered = ad.read_h5ad("inhAnno.h5ad")
		elif cellTypeToSubset=="Mic":
			adata_filtered = ad.read_h5ad("micAnno.h5ad")
		elif cellTypeToSubset=="Oli":
			adata_filtered = ad.read_h5ad("oliAnno.h5ad")
		elif cellTypeToSubset=="OPC":
			adata_filtered = ad.read_h5ad("opcAnno.h5ad")
		elif cellTypeToSubset=="Ast":
			adata_filtered = ad.read_h5ad("astAnno.h5ad")
		adata = adata_filtered.copy()
		adata.raw = adata_filtered.raw.to_adata()
		metaDF=adata.obs
		metaDF['batch'] = metaDF['batch'].astype(str)
		CSVPatient['projid'] = CSVPatient['projid'].astype(str)
		metaDF = metaDF.merge(
		CSVPatient[['projid', 'social_isolation_avg','msex']],how='left',left_on='batch',right_on='projid')
		metaDF.drop(columns=['projid'], inplace=True)
		adata = adata[metaDF["msex"] == 0].copy()
		adata = fixed_micropool(adata, patient_col="patient_id", pool_size=150)
		metaDF=adata.obs
		metaDF['batch'] = metaDF['batch'].astype(str)
		CSVPatient['projid'] = CSVPatient['projid'].astype(str)
		metaDF = metaDF.reset_index().merge(
		CSVPatient[['projid', 'social_isolation_avg','msex']],how='left',left_on='batch',right_on='projid')
		metaDF.drop(columns=['projid'], inplace=True)
		metaDF = metaDF.set_index('pool_id')
		if adata.shape[0] > 0:
			X_sparse = adata.X if sparse.issparse(adata.X) else sparse.csr_matrix(adata.X)
			X_sparse = X_sparse.T
			col_sums = X_sparse.sum(axis=0).A1 if hasattr(X_sparse, "A1") else np.array(X_sparse.sum(axis=0)).ravel()
			col_sums[col_sums == 0] = 1
			scaling_factors = 1e6 / col_sums
			D = sparse.diags(scaling_factors)
			X_sparseNew = X_sparse.dot(D)
			df1 = pd.DataFrame.sparse.from_spmatrix(X_sparseNew,index=adata.var_names,columns=adata.obs_names)
			df1.to_csv("Nfemale_matrix"+x+".tsv", sep="\t")
			X_sparse = X_sparse.T
			df = pd.DataFrame.sparse.from_spmatrix(X_sparse,index=adata.obs_names,columns=adata.var_names)
			# df = pd.DataFrame(X_sparse,index=adata.obs_names, columns=adata.var_names)
			adjacencies = grnboost2(expression_data=df, tf_names=tf_names)
			modules= list(modules_from_adjacencies(adjacencies, df))
			enrichedMot = prune2df(dbs, modules, folderPath+"motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl", num_workers=4)  # <-- safer to set single-threaded
			regulons = df2regulons(enrichedMot) 
			for r in regulons:
				if r.name.startswith("SREBF2"):
					print("Regulon:", r.name)
					print("Number of genes:", len(r.genes))
					print("Genes:", r.genes[:20])  # print first 20 genes
			with open(rootDirectory+"female__regulons_Tsai"+x+".p", "wb") as f:
				pickle.dump(regulons, f)
			dfReg = pd.DataFrame (regulons, columns = ['regulon'])
			#List of enriched regulons
			dfReg.to_csv("female_regulonsFULL_Tsai"+x+".csv")
			# AUCell
			auc_mtx = aucell(df, regulons, num_workers=4)
			auc_mtx.to_csv(rootDirectory+"female_auc_mtxFULL_Tsai"+x+".csv");
			auc_scores = auc_mtx
			results = []
			for regulon in auc_scores.columns:
			# if True:
				if True:
					# regulon="SREBF2(+)"
					print("hi")
					metadata=metaDF
					print(auc_mtx.head)
					print(auc_mtx[regulon])
					data = pd.DataFrame({'AUCell': auc_mtx[regulon],'Condition': metadata['social_isolation_avg'],'Batch': metadata['batch']})
					print(data.head)
					plt.figure(figsize=(6,4))
					sns.regplot(x="Condition",y="AUCell",data=data,scatter_kws={"s": 10, "alpha": 0.4},line_kws={"color": "red"})
					plt.title(x+" AUCell vs Condition")
					plt.xlabel("Condition (continuous)")
					plt.ylabel("AUCell score")
					plt.tight_layout()
					plt.savefig(x+"_female_AUCell_vs_Condition.png", dpi=300)
					try:
						print("hii")
						data["Condition"] = pd.to_numeric(data["Condition"], errors="coerce")
						model = smf.ols("AUCell ~ Condition", data=data)
						fit = model.fit()
						term = next((k for k in fit.params.keys() if "Condition" in k), None)
						coef = fit.params[term]
						pval0 = fit.pvalues[term]
						print("summary:"+str(regulon))
						print(fit.summary())
						print(regulon)
						results.append({"regulon": regulon,"coef": coef,"pval": pval0})
						print(regulon + ":" + str(pval))
					except Exception as e:
						results.append({"error": str(e)})
				results_df = pd.DataFrame(results).dropna(subset=["coef", "pval"])
				_, pvals_fdr, _, _ = multipletests(results_df["pval"], alpha=0.05, method='fdr_bh')
				results_df["fdr"] = pvals_fdr
				results_df["log_fdr"] = -np.log10(results_df["fdr"])
				results_df["significance"] = "Not Significant"
				results_df.loc[(results_df["fdr"] < 0.05) & (results_df["coef"] > 0), "significance"] = "Up"
				results_df.loc[(results_df["fdr"] < 0.05) & (results_df["coef"] < 0), "significance"] = "Down"
				results_df.to_csv("female_pVals"+x+".csv")
				palette = {"Up": "red","Down": "blue","Not Significant": "gray"}
				plt.figure(figsize=(10, 6))
				sns.scatterplot(
					data=results_df, x="coef", y="log_fdr", hue="significance",
					palette=palette, edgecolor=None
				)
				for _, row in results_df[(results_df["significance"] == "Up") | (results_df["significance"] == "Down")].iterrows():
					plt.text(
						row["coef"], row["log_fdr"], row["regulon"],
						fontsize=8, ha='right', va='bottom'
					)
				plt.axhline(-np.log10(0.05), linestyle="--", color="black", linewidth=1)
				plt.axvline(0, linestyle="--", color="black", linewidth=1)
				plt.xlabel("effect size")
				plt.ylabel("-log10(FDR)")
				plt.title(f"{x} Social Isolation Transcription Factor Activity - Volcano")
				plt.legend(title="Significance", loc="upper right")
				plt.tight_layout()
				plt.savefig(x + "TsaiFemaleVolcano.png") #volcano plot!
				plt.show()