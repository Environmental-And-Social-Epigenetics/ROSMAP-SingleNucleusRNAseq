import anndata as ad
import pandas as pd
file_path_target = 'totalAdataAnno3D.h5ad'
# variables = ['In_PV (Basket)', 'In_PV (Chandelier)']
#, 'Ex_L4_5''Ex_L5/6-CC', 'Ex_NRGN','Ast', 'Endo', 'Ex_L2_3', 'Ex_L4', 'Ex_L5', 'Ex_L5_6', 'In_Rosehip', 'In_SST', 'In_VIP', 'Mic', 'Oli', 'OPC'
from scipy.stats import ranksums
from scipy import sparse
from scipy.sparse import vstack
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
adataOG = ad.read_h5ad(file_path_target)

print(adataOG.shape)


import statsmodels.formula.api as smf
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
#
variables=["OPC","Oli","Exc","Inh","Ast","Mic"]
#


def fixed_micropool(adata, patient_col="patient_id", pool_size=30):
	#micropooling
	pooled_X = []
	pooled_obs = []
	obs_cols = [c for c in adata.obs.columns if c != patient_col]
	X=adata.raw.to_adata().X
	for patient in adata.obs[patient_col].unique():
		idx = np.where(adata.obs[patient_col] == patient)[0]
		n_cells = len(idx)
		n_pools = n_cells // pool_size
		for i in range(n_pools):
			pool_idx = idx[i*pool_size : (i+1)*pool_size]
			pooled_counts = X[pool_idx].sum(axis=0)
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
from statsmodels.stats.multitest import multipletests

#males
for x in variables:
	if __name__ == '__main__':
		#todo: put files into folder
		cellTypeToSubset =x
		#defining paths to files
		rootDirectory = os.path.join(os.environ.get("ACE_OUTPUT_ROOT", "."), "SCENIC", "DeJager", "") 
		os.chdir(rootDirectory)
		folderPath=os.environ.get("SCENIC_RANKING_DIR", ".") #change path to folder
		db_paths = [folderPath+"hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
		folderPath+"hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"]
		dbs = [RankingDatabase(fname=p, name=os.path.basename(p)) for p in db_paths]
		f_tfs = folderPath+"hg.txt"
		tf_names = load_tf_names(f_tfs)
		CSVPatient = pd.read_csv(os.path.join(os.environ.get('PHENOTYPE_DIR', '.'), 'TSAI_DEJAGER_all_patients_wACEscores.csv'))
		adata_filtered=adataOG
		print(adata_filtered.obs['cell_type'])
		#cell type subsetting
		if cellTypeToSubset=="Exc":
			adata_filtered = adataOG[adataOG.obs['cell_type'].isin(["Ex_L4_5", "Ex_L5_6", "Ex_L5_6_CC", "Ex_L2_3", "Ex_L4","Ex_L5"])]
		elif cellTypeToSubset=="Inh":
			adata_filtered = adataOG[adataOG.obs['cell_type'].isin(["In_VIP", "In_PV (Basket)", "In_PV (Chandelier)", "In_SST", "In_Rosehip"])]
		else:
			adata_filtered = adataOG[adataOG.obs['cell_type']==cellTypeToSubset]
		adata = adata_filtered.to_memory()
		sc.pp.filter_genes(adata, min_cells=200)
		print(adata)
		adata.raw = adata_filtered.raw.to_adata()[:, adata.var_names]
		metaDF=adata.obs #metadata
		metaDF['patient_id'] = metaDF['patient_id'].astype(str)
		CSVPatient['projid'] = CSVPatient['projid'].astype(str)
		metaDF = metaDF.merge(
		CSVPatient[['projid', 'ace_aggregate','msex']],how='left',left_on='patient_id',right_on='projid')
		metaDF.drop(columns=['projid'], inplace=True)
		adata = adata[metaDF["msex"] == 1] #subset for sex
		print(adata)
		adata = fixed_micropool(adata, patient_col="patient_id", pool_size=30) #micropooling
		metaDF=adata.obs #regenerate metadata
		metaDF['patient_id'] = metaDF['patient_id'].astype(str)
		CSVPatient['projid'] = CSVPatient['projid'].astype(str)
		metaDF = metaDF.reset_index().merge(
		CSVPatient[['projid', 'ace_aggregate','msex']],how='left',left_on='patient_id',right_on='projid')
		metaDF.drop(columns=['projid'], inplace=True)
		metaDF = metaDF.set_index('pool_id')
		metaDF.to_csv("male_metaDF"+x+".csv")
		if adata.shape[0] > 0:
			X_sparse = adata.X if sparse.issparse(adata.X) else sparse.csr_matrix(adata.X)
			print(adata)
			# col_sums = X_sparse.sum(axis=0).A1 if hasattr(X_sparse, "A1") else np.array(X_sparse.sum(axis=0)).ravel()
			# col_sums[col_sums == 0] = 1
			# scaling_factors = 1e6 / col_sums
			# D = sparse.diags(scaling_factors)
			# X_sparseNew = X_sparse.dot(D)
			# df1 = pd.DataFrame.sparse.from_spmatrix(X_sparseNew,index=adata.obs_names,columns=adata.var_names)
			# df1.to_csv("male_matrix"+x+".tsv", sep="\t")
			df = pd.DataFrame.sparse.from_spmatrix(X_sparse,index=adata.obs_names,columns=adata.var_names)
			#SCENIC steps
			adjacencies = grnboost2(expression_data=df, tf_names=tf_names) #GRNboost
			print("hi")
			modules= list(modules_from_adjacencies(adjacencies, df))
			enrichedMot = prune2df(dbs, modules, folderPath+"motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl", num_workers=1)  # <-- safer to set single-threaded
			regulons = df2regulons(enrichedMot) 
			print("hi")
			with open(rootDirectory+"male__regulons_Tsai"+x+".p", "wb") as f:
				pickle.dump(regulons, f)
			dfReg = pd.DataFrame (regulons, columns = ['regulon'])
			#List of enriched regulons
			dfReg.to_csv("male_regulonsFULL_Tsai"+x+".csv") #male regulons
			print("hi")
			# AUCell
			auc_mtx = aucell(df, regulons, num_workers=4)
			auc_mtx.to_csv(rootDirectory+"male_auc_mtxFULL_Tsai"+x+".csv");
			auc_scores = auc_mtx
			results = []
			if True:
				for regulon in auc_scores.columns:
					metadata=metaDF
					data = pd.DataFrame({'AUCell': auc_mtx[regulon],'Condition': metadata['ace_aggregate'],'Batch': metadata['batch']})
					plt.figure(figsize=(6,4))
					sns.regplot(x="Condition",y="AUCell",data=data,scatter_kws={"s": 10, "alpha": 0.4},line_kws={"color": "red"})
					plt.title(x+" AUCell vs Condition")
					plt.xlabel("Condition (continuous)")
					plt.ylabel("AUCell score")
					plt.tight_layout()
					plt.savefig(x+"_male_AUCell_vs_Condition.png", dpi=300)
					try:
						data["Condition"] = pd.to_numeric(data["Condition"], errors="coerce")
						model = smf.ols("AUCell ~ Condition", data=data)
						#statistical testing
						fit = model.fit()
						term = next((k for k in fit.params.keys() if "Condition" in k), None)
						coef = fit.params[term]
						pval0 = fit.pvalues[term]
						results.append({"regulon": regulon,"coef": coef,"pval": pval0})
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
				plt.savefig(x + "MaleVolcano.png") #volcano plot!
				plt.show()


variables=["Exc","Inh","Ast","Mic","OPC","Oli"]

for x in variables:
	if __name__ == '__main__':
		# variables=['Ex_L2_3']
		#todo: put files into folder
		cellTypeToSubset =x
		#defining paths to files
		rootDirectory = os.path.join(os.environ.get("ACE_OUTPUT_ROOT", "."), "SCENIC", "DeJager", "") 
		os.chdir(rootDirectory)
		folderPath=os.environ.get("SCENIC_RANKING_DIR", ".") #change path to folder
		db_paths = [folderPath+"hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
		folderPath+"hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"]
		dbs = [RankingDatabase(fname=p, name=os.path.basename(p)) for p in db_paths]
		f_tfs = folderPath+"hg.txt"
		tf_names = load_tf_names(f_tfs)
		CSVPatient = pd.read_csv(os.path.join(os.environ.get('PHENOTYPE_DIR', '.'), 'TSAI_DEJAGER_all_patients_wACEscores.csv'))
		adata_filtered=adataOG
		if cellTypeToSubset=="Exc":
			adata_filtered = adataOG[adataOG.obs['cell_type'].isin(["Ex_L4_5", "Ex_L5_6", "Ex_L5_6_CC", "Ex_L2_3", "Ex_L4","Ex_L5"])]
		elif cellTypeToSubset=="Inh":
			adata_filtered = adataOG[adataOG.obs['cell_type'].isin(["In_VIP", "In_PV (Basket)", "In_PV (Chandelier)", "In_SST", "In_Rosehip"])]
		else:
			adata_filtered = adataOG[adataOG.obs['cell_type']==cellTypeToSubset]
		adata = adata_filtered
		adata.raw = adata_filtered.raw.to_adata()[:, adata.var_names]
		metaDF=adata.obs
		metaDF['patient_id'] = metaDF['patient_id'].astype(str)
		CSVPatient['projid'] = CSVPatient['projid'].astype(str)
		metaDF = metaDF.merge(
		CSVPatient[['projid', 'ace_aggregate','msex']],how='left',left_on='patient_id',right_on='projid')
		metaDF.drop(columns=['projid'], inplace=True)
		adata = adata[metaDF["msex"] == 0]
		adata = fixed_micropool(adata, patient_col="patient_id", pool_size=39)
		metaDF=adata.obs
		metaDF['patient_id'] = metaDF['patient_id'].astype(str)
		CSVPatient['projid'] = CSVPatient['projid'].astype(str)
		metaDF = metaDF.reset_index().merge(
		CSVPatient[['projid', 'ace_aggregate','msex']],how='left',left_on='patient_id',right_on='projid')
		metaDF.drop(columns=['projid'], inplace=True)
		metaDF = metaDF.set_index('pool_id')
		metaDF.to_csv("female_metaDF"+x+".csv")
		if adata.shape[0] > 0:
			X_sparse = adata.X if sparse.issparse(adata.X) else sparse.csr_matrix(adata.X)
			# col_sums = X_sparse.sum(axis=0).A1 if hasattr(X_sparse, "A1") else np.array(X_sparse.sum(axis=0)).ravel()
			# col_sums[col_sums == 0] = 1
			# scaling_factors = 1e6 / col_sums
			# D = sparse.diags(scaling_factors)
			# X_sparseNew = X_sparse.dot(D)
			# df1 = pd.DataFrame.sparse.from_spmatrix(X_sparseNew,index=adata.obs_names,columns=adata.var_names)
			# df1.to_csv("female_matrix"+x+".tsv", sep="\t")
			df = pd.DataFrame.sparse.from_spmatrix(X_sparse,index=adata.obs_names,columns=adata.var_names)
			# df = pd.DataFrame(X_sparse,index=adata.obs_names, columns=adata.var_names)
			adjacencies = grnboost2(expression_data=df, tf_names=tf_names) #grnboost
			modules= list(modules_from_adjacencies(adjacencies, df))
			enrichedMot = prune2df(dbs, modules, folderPath+"motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl", num_workers=1)  # <-- safer to set single-threaded
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
			auc_mtx.to_csv(rootDirectory+"female_auc_mtxFULL_Tsai"+x+".csv")
			auc_scores = auc_mtx
			results = []
			if True:
				for regulon in auc_scores.columns:
					metadata=metaDF
					#statistical testing
					data = pd.DataFrame({'AUCell': auc_mtx[regulon],'Condition': metadata['ace_aggregate'],'Batch': metadata['batch']})
					plt.figure(figsize=(6,4))
					sns.regplot(x="Condition",y="AUCell",data=data,scatter_kws={"s": 10, "alpha": 0.4},line_kws={"color": "red"})
					plt.title(x+" AUCell vs Condition")
					plt.xlabel("Condition (continuous)")
					plt.ylabel("AUCell score")
					plt.tight_layout()
					plt.savefig(x+"_female_AUCell_vs_Condition.png", dpi=300)
					try:
						data["Condition"] = pd.to_numeric(data["Condition"], errors="coerce")
						model = smf.ols("AUCell ~ Condition", data=data)
						fit = model.fit()
						term = next((k for k in fit.params.keys() if "Condition" in k), None)
						coef = fit.params[term]
						pval0 = fit.pvalues[term]
						results.append({"regulon": regulon,"coef": coef,"pval": pval0}) 
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
				plt.savefig(x + "FemaleVolcano.png") #volcano plot!
				plt.show()
#run for tsai and dejager



# if __name__ == '__main__':
#	 for x in variables:
#		 print(x)
#		 results = []
#		 adata_filtered = adataOG[adataOG.obs['cell_type']==x]
#		 CSVPatient = pd.read_csv(os.path.join(os.environ.get('PHENOTYPE_DIR', '.'), 'TSAI_DEJAGER_all_patients_wACEscores.csv'))
#		 adata = adata_filtered.copy()
#		 # sc.pp.filter_genes(adata, min_cells=2000, inplace=True)
#		 adata.raw = adata_filtered.raw.to_adata()[:, adata.var_names]
#		 X_dense = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
#		 df = pd.DataFrame(X_dense,index=adata.obs_names, columns=adata.var_names)
#		 metaDF=adata.obs
#		 metaDF['barcode'] = adata.obs_names
#		 metaDF['batch'] = metaDF['batch'].astype(str)
#		 CSVPatient['projid'] = CSVPatient['projid'].astype(str)
#		 metaDF = metaDF.merge(CSVPatient[['projid', 'ace_aggregate','msex']],how='left',left_on='batch',right_on='projid')
#		 metaDF.drop(columns=['projid'], inplace=True)
#		 metaDF.index = metaDF['barcode']
#		 adata = adata[metaDF["msex"] == 1].copy()
#		 adata = fixed_micropool(adata, patient_col="patient_id", pool_size=30)
#		 metaDF=adata.obs
#		 metaDF['barcode'] = adata.obs_names
#		 metaDF['batch'] = metaDF['batch'].astype(str)
#		 CSVPatient['projid'] = CSVPatient['projid'].astype(str)
#		 metaDF = metaDF.merge(CSVPatient[['projid', 'ace_aggregate','msex']],how='left',left_on='batch',right_on='projid')
#		 metaDF.drop(columns=['projid'], inplace=True)
#		 metaDF.index = metaDF['barcode']
#		 print("hi")
#		 if os.path.exists("male__regulons_Tsai" + x + ".p"):
#			 print("hi")
#			 with open("male__regulons_Tsai" + x + ".p", "rb") as f:
#				 regulons = pickle.load(f)
#			 if __name__ == '__main__':
#				 print("hi")
#				 auc_scores = aucell(df, regulons, num_workers=2)
#				 metadata = metaDF.loc[auc_scores.index]
#				 for regulon in auc_scores.columns:
#					 if True:
#						 data = pd.DataFrame({
#							 'AUCell': auc_scores[regulon],
#							 'Condition': metadata['ace_aggregate'],
#							 'Batch': metadata['batch']
#						 })
#						 data["Condition"] = data["Condition"].astype(str)
#						 try:
#							 model = smf.mixedlm("AUCell ~ Condition", data, groups=data["Batch"])
#							 fit = model.fit()
#							 print(fit.params.keys())
#							 term = next((k for k in fit.params.keys() if "Condition" in k and "T." in k), None)
#							 # print(term)
#							 coef = fit.params[term]
#							 pval = fit.pvalues[term]
#							 results.append({
#								 "regulon": regulon,
#								 "coef": coef,
#								 "pval": pval
#							 })
#							 print(regulon + ":" + str(pval))
#						 except Exception as e:
#							 results.append({"error": str(e)})
#				 results_df = pd.DataFrame(results).dropna(subset=["coef", "pval"])
#				 print(results_df)
#				 results_df["log_pval"] = -np.log10(results_df["pval"])
#				 results_df["significant"] = (results_df["pval"] < 0.05) 
#				 plt.figure(figsize=(10, 6))
#				 sns.scatterplot(
#					 data=results_df, x="coef", y="log_pval", hue="significant",
#					 palette={True: "red", False: "gray"}, edgecolor=None
#				 )
#				 for _, row in results_df[results_df["significant"]].iterrows():
#					 plt.text(
#						 row["coef"], row["log_pval"], row["regulon"],
#						 fontsize=8, ha='right', va='bottom'
#					 )
#				 plt.axhline(-np.log10(0.05), linestyle="--", color="black", linewidth=1)
#				 plt.axvline(0.05, linestyle="--", color="black", linewidth=1)
#				 plt.axvline(-0.05, linestyle="--", color="black", linewidth=1)
#				 plt.xlabel("Effect size for AUCell Scoring")
#				 plt.ylabel("-log10(p-value)")
#				 plt.title(f"Volcano Plot: Resilience Relation to Regulon Activity ({x})")
#				 plt.legend(title="Significant", loc="upper right")
#				 plt.tight_layout()
#				 plt.savefig(x + "TsaiMaleVolcano.png")
#				 plt.show()

