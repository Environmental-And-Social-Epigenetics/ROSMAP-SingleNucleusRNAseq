import pandas as pd
import os
import numpy as np
import seaborn as sns
import scanpy as sc
from scipy.stats import median_abs_deviation

obj = sc.read_10x_h5('/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/Tsai/Preprocessing/Preprocessed_Counts/ACE/14184286/processed_feature_bc_matrix_filtered.h5')
obj.var_names_make_unique()
obj.var["mt"] = obj.var_names.str.startswith("MT-")
# ribosomal genes
obj.var["ribo"] = obj.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes.
obj.var["hb"] = obj.var_names.str.contains(("^HB[^(P)]"))

#Calculate the qc metrics
sc.pp.calculate_qc_metrics(obj, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)

# Function that takes as input number of MADS that is permissible and metric of interest (column from .obs)
# and returns whether it is an outlier

def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric] #pull column of values of interest (metric we are looking at)
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    ) # 'outliers' are are values who are plus or minus nmads (which we specify) from the median
    return outlier

#outlier status for non mito metrics --> use MAD for this
obj.obs["outlier"] = (
    is_outlier(obj, "log1p_total_counts", 5)
    | is_outlier(obj, "log1p_n_genes_by_counts", 5)
    | is_outlier(obj, "pct_counts_in_top_20_genes", 5)
)

#outlier status for mito metrics
obj.obs["mt_outlier"] = obj.obs["pct_counts_mt"] > 10 #classify mitochondrial outliers as cells with >10 percent mito

obj = obj[(~obj.obs.outlier) & (~obj.obs.mt_outlier)].copy() #subset out cells that are not outliers or mitochondrial outliers

obj.write_h5ad('/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/QC/Outliers/14184286/14184286_outlier_filtered.h5')
