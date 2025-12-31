import pandas as pd
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

# Replace 'your_directory_path' with the path to the directory you want to scan
directory_path = '/net/vast-storage/scratch/vast/lhtsai/mabdel03/files/ACE_Analysis/Data/DeJager/Preprocessed_Counts/'

# List all folders in the directory
libraries = [folder for folder in os.listdir(directory_path) if os.path.isdir(os.path.join(directory_path, folder))]

libraries.remove("qc_plots")
libraries.remove("figures")
libraries.remove("concatObjFinal.zarr")
libraries.remove("OutputP2")
# print("Folders in the directory:")
# for folder in folders:
#     print(folder)

# libraries = ['200227-B12-B','200930-B55-A','200826-B49-A','200806-B44-B','200225-B10-A','190403-B4-B','190409-B5-A','200306-B16-A','200312-B20-A','200313-B23-B']
for library in libraries:
    print(library)
    root = library
    library1 = os.path.join(root, 'processed_feature_bc_matrix_filtered.h5')
    if os.path.exists(library1):
        library_adata = sc.read_10x_h5(library1)
        original_adata = library_adata

        library_adata.var_names_make_unique()

        # mitochondrial genes
        library_adata.var["mt"] = library_adata.var_names.str.startswith("MT-")
        library_adata.var["ribo"] = library_adata.var_names.str.startswith(("RPS", "RPL"))
        library_adata.var["hb"] = library_adata.var_names.str.contains(("^HB[^(P)]"))

        # QC Metrics
        sc.pp.calculate_qc_metrics(library_adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)

        # counts_hist, bin_edges = np.histogram(library_adata.obs["total_counts"], bins=50)
        # peaks, _ = find_peaks(counts_hist)

        # # Identify the main peak
        # main_peak_idx = peaks[np.argmax(counts_hist[peaks])]
        # threshold = bin_edges[main_peak_idx] + 300 
        # print(f"Count Depth Threshold: {threshold}")


        # genes_hist, gene_bin_edges = np.histogram(library_adata.obs["n_genes_by_counts"], bins=50)
        # noise_peaks, _ = find_peaks(genes_hist, height=5)  # Low peak for noise
        # noise_peak_idx = noise_peaks[0]  # Assume first peak is noise

        # # Set threshold above noise peak
        # gene_threshold = gene_bin_edges[noise_peak_idx] + 100
        # print(f"Gene Count Threshold: {gene_threshold}")

        # from kneed import KneeLocator

        # # Example data: Replace with actual count-depth data
        # sorted_counts = np.sort(library_adata.obs["total_counts"])[::-1]

        # #maximum second derivative of total count plot
        # sorted_counts = np.array(sorted_counts) 
        # first_derivative = np.diff(sorted_counts)
        # second_derivative = np.diff(first_derivative)
        # elbow_point = np.argmax(second_derivative) + 1
        # elbow_threshold = sorted_counts[elbow_point+10]

        # print(f"Second Derivative Elbow Threshold: {elbow_threshold}")

        # low_quality_mask = library_adata.obs["pct_counts_mt"] > 0.09
        # high_quality_mask = ~low_quality_mask

        # #mask for filtration
        # count_threshold = np.percentile(library_adata.obs["total_counts"][high_quality_mask], 5)  # Lower Xth percentile of total counts of cells who have percent mito below 10% 
        # gene_threshold = np.percentile(library_adata.obs["n_genes_by_counts"][high_quality_mask], 5) # Lower Xth percentile of number of genes per cell of cells who have percent mito below 10%  

        # print(f"New Count Threshold: {count_threshold}")
        # print(f"New Gene Threshold: {gene_threshold}")

        # pct_threshold = 80
        # plt.figure(figsize=(10, 6))
        # plt.hist(library_adata.obs["total_counts"], bins=50, alpha=0.7, label="All Cells")
        # plt.axvline(count_threshold, color="red", linestyle="--", label=f"Threshold: {count_threshold}")
        # plt.title("Count Depth per Cell") 
        # plt.xlabel("Total Counts")
        # plt.ylabel("Frequency")
        # plt.legend()
        # plt.savefig("histogram_count_depth.png", dpi=300)
        # plt.close()

        # plt.figure(figsize=(10, 6))
        # plt.hist(library_adata.obs["n_genes_by_counts"], bins=50, alpha=0.7, label="All Cells")
        # plt.axvline(gene_threshold, color="red", linestyle="--", label=f"Threshold: {gene_threshold}")
        # plt.title("Number of Genes Detected per Cell")
        # plt.xlabel("Number of Genes")
        # plt.ylabel("Frequency")
        # plt.legend()
        # plt.savefig("histogram_genes_detected.png", dpi=300)
        # plt.close()

        # sorted_counts = np.sort(library_adata.obs["total_counts"])[::-1]
        # plt.figure(figsize=(10, 6))
        # plt.plot(sorted_counts, label="Count Depths")
        # plt.axhline(count_threshold, color="red", linestyle="--", label=f"Threshold: {count_threshold}")
        # plt.yscale("log")
        # plt.xscale("log")
        # plt.title("Count Depth Distribution (High to Low)")
        # plt.xlabel("Cell Rank")
        # plt.ylabel("Total Counts (log scale)")
        # plt.legend()
        # plt.savefig("count_depth_distribution.png", dpi=300)
        # plt.close()

        # plt.figure(figsize=(10, 6))
        # scatter = plt.scatter(
        #     library_adata.obs["total_counts"],
        #     library_adata.obs["n_genes_by_counts"],
        #     c=library_adata.obs["pct_counts_mt"],
        #     cmap="viridis",
        #     alpha=0.7
        # )
        # plt.axvline(count_threshold, color="red", linestyle="--", label=f"Count Threshold: {count_threshold}")
        # plt.axhline(gene_threshold, color="blue", linestyle="--", label=f"Gene Threshold: {gene_threshold}")
        # plt.colorbar(scatter, label="Fraction of Mitochondrial Reads")
        # plt.title("Number of Genes vs. Count Depth")
        # plt.xlabel("Total Counts")
        # plt.ylabel("Number of Genes")
        # plt.legend()
        # plt.savefig("genes_vs_counts_mito.png", dpi=300)
        # plt.close()

        # # QC Plots
        # print("LIBRARY METRICS")
        # p1 = sns.displot(library_adata.obs["total_counts"], bins=100, kde=False)
        # p1.savefig("total_counts_distribution.png", dpi=300)
        # plt.close() 

        # #violin plots!!
        # sc.pl.violin(
        #     library_adata, 
        #     "pct_counts_mt", 
        #     save="pct_counts_mt_violin.png"
        # )
        # sc.pl.violin(
        #     library_adata, 
        #     "log1p_total_counts", 
        #     save="log1p_total_counts.png"
        # )
        # sc.pl.violin(
        #     library_adata, 
        #     "n_genes_by_counts", 
        #     save="n_genes_by_counts.png"
        # )




        # def is_outlier_percentile(adata, metric: str, lower_percentile: float, upper_percentile: float):
        #     M = adata.obs[metric]
        #     lower_thresh = np.percentile(M, lower_percentile)
        #     upper_thresh = np.percentile(M, upper_percentile)
        #     return (M < lower_thresh) | (M > upper_thresh)

        def is_outlier(adata, metric: str, nmads: int):
            M = adata.obs[metric] #pull column of values of interest (metric we are looking at)
            outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
                np.median(M) + nmads * median_abs_deviation(M) < M
            ) # 'outliers' are are values who are plus or minus nmads (which we specify) from the median
            return outlier

        print(library_adata.obs['total_counts'].describe())

        print(library_adata.obs['n_genes_by_counts'].describe())

        print(library_adata.obs['pct_counts_in_top_20_genes'].describe())
        print(library_adata.obs['pct_counts_mt'].describe())

        # # Apply Outlier Detection

        # library_adata.obs["top_20_outlier"] = library_adata.obs['pct_counts_in_top_20_genes'] >= 40 #Top Gene Dominance

        
        # Create masks for cells meeting the thresholds

        # count_threshold = count_threshold
        # gene_threshold = gene_threshold
        # count_mask = library_adata.obs["total_counts"] > count_threshold
        # gene_mask = library_adata.obs["n_genes_by_counts"] > gene_threshold
        # pctCount_mask = library_adata.obs["pct_counts_in_top_20_genes"] <= pct_threshold

        #Set thresholds for everything

        # Combine masks (logical AND)
        # filter_mask = count_mask & gene_mask & pctCount_mask 

        # Filter the AnnData object
        # library_adata = library_adata[filter_mask].copy()

        # library_adata.obs["outlier"] = (
        #     is_outlier_percentile(library_adata, "log1p_total_counts",3,97)#6.7,90.9)
        #     | is_outlier_percentile(library_adata, "log1p_n_genes_by_counts",3,97)#6.7,91.7)
        #     | is_outlier_percentile(library_adata, "pct_counts_in_top_20_genes",3,97)#6.7,91.7)
           
        # )
        # # | is_outlier(library_adata, "pct_counts_in_top_20_genes",4)#0,90)

        # # Mitochondrial Outliers
        library_adata.obs["outlier"] = (
            is_outlier(library_adata, "log1p_total_counts",4)#6.7,90.9)
            | is_outlier(library_adata, "log1p_n_genes_by_counts",4)#6.7,91.7)
            | is_outlier(library_adata, "pct_counts_in_top_20_genes",4)#6.7,91.7)
        )
        # library_adata.obs["mt_outlier"] = library_adata.obs["pct_counts_mt"] > 10


        # library_adata.obs["outlier"] = (
        #     is_outlier(library_adata, "log1p_total_counts",4)#6.7,90.9)
        #     | is_outlier(library_adata, "log1p_n_genes_by_counts",4)#6.7,91.7)
        #     | is_outlier(library_adata, "pct_counts_in_top_20_genes",4)#6.7,91.7)
        # )
        library_adata.obs["mt_outlier"] = is_outlier(library_adata, "pct_counts_mt", 3) | (
            library_adata.obs["pct_counts_mt"] > 7.5
        )

        # # Additional Filters (Optional)
        # # Set thresholds dynamically based on visualization
        library_adata = library_adata[
            (~library_adata.obs["outlier"]) &
            (~library_adata.obs["mt_outlier"])
            # (~library_adata.obs["top_20_outlier"])
        ].copy()

        # library_adata = library_adata[
        #     (~library_adata.obs["mt_outlier"]) &
        #     (~library_adata.obs["top_20_outlier"])
        # ].copy()
        # min_cells = 5  #minimum cells a gene should be expressed in

        # #number of cells each gene is expressed in
        # gene_mask = np.sum(library_adata.X > 0, axis=0) >= min_cells

        # # filtering out genes that aren't expressed enough
        # library_adata = library_adata[:, gene_mask].copy()

        # min_genes = 200

        # #number of genes expressed in each cell 
        # n_genes_per_cell = np.sum(library_adata.X > 0, axis=1)

        # # filtering out cells that don't expressed enough
        # cell_mask = n_genes_per_cell >= min_genes

        # # filtering cells
        # library_adata = library_adata[cell_mask].copy()
        # Plot the calculated QC Metrics
        # print("LIBRARY METRICS")
        # p1 = sns.displot(library_adata.obs["total_counts"], bins=100, kde=False)
        # p2 = sc.pl.violin(library_adata, "pct_counts_mt")
        # p3 = sc.pl.scatter(library_adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

        # Function that takes as input number of MADS that is permissible and metric of interest (column from .obs)
        # and returns whether it is an outlier

        # def is_outlier_percentile(adata, metric: str, lower_percentile: float = 5, upper_percentile: float = 95):
        #     M = adata.obs[metric]
        #     lower_thresh = np.percentile(M, lower_percentile)
        #     upper_thresh = np.percentile(M, upper_percentile)
            
        #     outlier = (M < lower_thresh) | (M > upper_thresh)
        #     return outlier

        # # # Add outlier status to every cell in every subject
        # # #outlier status for non mito metrics --> use MAD for this
        # # #outlier status for non mito metrics --> use MAD for this
        # library_adata.obs["outlier"] = (
        #     is_outlier_percentile(library_adata, "log1p_total_counts")
        #     | is_outlier_percentile(library_adata, "log1p_n_genes_by_counts")
        #     | is_outlier_percentile(library_adata, "pct_counts_in_top_20_genes")
        # )

        # #outlier status for mito metrics
        # #sc best practices included: is_outlier(library_adata, "pct_counts_mt", 3) |
        # library_adata.obs["mt_outlier"] = library_adata.obs["pct_counts_mt"] > 10 #classify mitochondrial outliers as cells with >10 percent mito

        # # Filter the objects based on the constructed thresholds/outlier classifications

        # library_adata = library_adata[(~library_adata.obs.outlier) & (~library_adata.obs.mt_outlier)].copy() #subset out cells that are not outliers or mitochondrial outliers
        # p1 = sc.pl.scatter(library_adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt") # make plots

        metrics = ["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt","pct_counts_in_top_20_genes"]

        for metric in metrics:
            print(f"Metric: {metric}")
            print("Before Filtering:")
            print(original_adata.obs[metric].describe())
            print("\nAfter Filtering:")
            print(library_adata.obs[metric].describe())
            print("-" * 50)


        filename = library+'OutputP1.h5ad'
        library_adata.write_h5ad(filename)