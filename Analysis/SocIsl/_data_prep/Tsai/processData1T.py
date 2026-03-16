import pandas as pd
import anndata as ad
import os
import numpy as np
import seaborn as sns
import scanpy as sc
import torch
import tempfile
from scipy.stats import median_abs_deviation
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
#Figure settings --> can adjust as desired but these were the ones used in SC Best Practices
sc.settings.verbosity = 0
os.chdir('/net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/Preprocessed_Counts/')

sc.set_figure_params(figsize=(6, 6), frameon=False)
sns.set_theme()
torch.set_float32_matmul_precision("high")
save_dir = tempfile.TemporaryDirectory()
plot_dir = './qc_plots'
os.makedirs(plot_dir, exist_ok=True)
sc.settings.figdir = plot_dir

import os
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

# Replace 'your_directory_path' with the path to the directory you want to scan

# List all folders in the directory
# libraries = ['200908-B50-A', '200825-B48-B','191219-B9-B','200303-B14-B','200707-B30-B','200701-B28-B','200313-B22-A','200227-B12-A','190403-B4-B','201007-B57-B','200930-B55-B','200915-B53-A','200909-B51-B','200810-B46-B']
libraries = ['86408244', '50102088']

# libraries = [folder for folder in os.listdir(directory_path) if os.path.isdir(os.path.join(directory_path, folder))]

# libraries.remove("2518573")
# libraries.remove("82317494")
# libraries.remove("83034844")
# libraries.remove("94430339")
# libraries.remove("95919181")
# libraries.remove("52311825")
# libraries.remove("50101523")
# libraries.remove("65499271")

# library_adataist = []
# for library in libraries:
#     root = library
#     libraryPath = os.path.join(root, 'processed_feature_bc_matrix_filtered.h5')
#     if os.path.exists(libraryPath):
#         library_adata = sc.read_10x_h5(libraryPath)
#         # Make indices unique
#         library_adata.obs_names_make_unique()
#         library_adata.var_names_make_unique()
#         # library_adata.obs['patient_id'] = library_adata.obs['patient_id'].apply(lambda x: str(int(float(x))) if pd.notna(x) else '')
#         # print(library_adata.obs['patient_id'].head())
#         library_adataist.append(library_adata)
#     #how about tag patients now??

# # Concatenate after ensuring unique indices
# adata = ad.concat(library_adataist, join="outer")

# print("shape originally")
# print(adata.shape)


# libraries.remove("figures")
# libraries.remove("concatObjFinal.zarr")
# libraries.remove("OutputP2")
# print("Folders in the directory:")
# for folder in folders:
#     print(folder)

# libraries = ['200227-B12-B','200930-B55-A','200826-B49-A','200806-B44-B','200225-B10-A','190403-B4-B','190409-B5-A','200306-B16-A','200312-B20-A','200313-B23-B']
#library = '200313-B22-B'
#if True:

for library in libraries:
    print(library)
    root = library
    if os.path.exists(os.path.join("/om/scratch/Mon/mabdel03/SocIsl/Preprocessed_Counts/SocIsl/"+root, 'processed_feature_bc_matrix_filtered.h5')):
        libraryPath=os.path.join("/om/scratch/Mon/mabdel03/SocIsl/Preprocessed_Counts/SocIsl/"+root, 'processed_feature_bc_matrix_filtered.h5')
    elif os.path.exists(os.path.join("/net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/Tsai/Preprocessing/Preprocessed_Counts/ACE/"+root, 'processed_feature_bc_matrix_filtered.h5')):
        libraryPath=os.path.join("/net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/Tsai/Preprocessing/Preprocessed_Counts/ACE/"+root, 'processed_feature_bc_matrix_filtered.h5')
    elif os.path.exists(os.path.join("/net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/Tsai/Preprocessing/Preprocessed_Counts/Resilient/"+root, 'processed_feature_bc_matrix_filtered.h5')):
        libraryPath=os.path.join("/net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/Tsai/Preprocessing/Preprocessed_Counts/Resilient/"+root, 'processed_feature_bc_matrix_filtered.h5')
    if True:
        library_adata = sc.read_10x_h5(libraryPath)
        original_adata = library_adata
        if len(library_adata.obs_names) !=0:
            library_adata.var_names_make_unique()
            # mitochondrial genes
            library_adata.var["mt"] = library_adata.var_names.str.startswith("MT-")
            library_adata.var["ribo"] = library_adata.var_names.str.startswith(("RPS", "RPL"))
            library_adata.var["hb"] = library_adata.var_names.str.contains(("^HB[^(P)]"))
            # QC Metrics
            sc.pp.calculate_qc_metrics(library_adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)
            # counts_hist, bin_edges = np.histogram(library_adata.obs["total_counts"], bins=50)
            # peaks, _ = find_peaks(counts_hist)
            print(library_adata)
            def is_outlier_percentile(adata, metric: str, lower_percentile: float, upper_percentile: float):
                M = adata.obs[metric]
                lower_thresh = np.percentile(M, lower_percentile)
                upper_thresh = np.percentile(M, upper_percentile)
                return (M < lower_thresh) | (M > upper_thresh)
            def is_low_outlier(adata, metric: str, nmads: int):
                #low outlier as opposed to bidirectional!
                M = adata.obs[metric]
                outlier = M < np.median(M) - nmads * median_abs_deviation(M)
                return outlier
            def is_outlier(adata, metric: str, nmads: int):
                M = adata.obs[metric] #pull column of values of interest (metric we are looking at)
                outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
                    np.median(M) + nmads * median_abs_deviation(M) < M
                ) # 'outliers' are are values who are plus or minus nmads (which we specify) from the median
                return outlier
            print(library_adata.obs)
            sc.pl.violin(library_adata, ["total_counts", "n_genes_by_counts", "pct_counts_mt"], jitter=0.4, multi_panel=True)
            library_adata.obs["outlier"] = (
                is_outlier_percentile(library_adata, "log1p_total_counts", 4.5,96) |  # Low counts outliers
                is_outlier_percentile(library_adata, "log1p_n_genes_by_counts", 5,100)  # Low genes outliers
                # (library_adata.obs["log1p_total_counts"] > 11)  # Absolute maximum total counts
            )
            library_adata.obs["mt_outlier"] = library_adata.obs["pct_counts_mt"] > 10
            library_adata = library_adata[
                (~library_adata.obs["outlier"]) &
                (~library_adata.obs["mt_outlier"])
                # (~library_adata.obs["top_20_outlier"])
            ].copy()
            metrics = ["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt","pct_counts_in_top_20_genes"]
            filename = '/om/scratch/Mon/mabdel03/SocialIsolation/Tsai/'+library+'OutputP1.h5ad'
            library_adata.write_h5ad(filename)
            # filename2 = library+'OutputP1B.h5ad'
            # library_adata2.write_h5ad(filename2)
            # filename3 = library+'OutputP1C.h5ad'
            # library_adata3.write_h5ad(filename3)