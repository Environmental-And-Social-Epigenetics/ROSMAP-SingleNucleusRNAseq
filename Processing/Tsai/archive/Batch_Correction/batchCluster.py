#Imports
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import os
from functools import reduce
import seaborn as sns
import torch
from rich import print
from pathlib import Path
import matplotlib.pyplot as plt
from sklearn.metrics import silhouette_score
from scipy.stats import chi2_contingency
from gtfparse import read_gtf
import polars as pl


#Load in anndata object!
adata = ad.read_h5ad('/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/clearedObj2.h5ad')

#First, depth!
def calculate_transcriptome_size(gtf_file):
    #loading GTF file in
    gtf = pl.read_csv(
        gtf_file,
        separator="\t", 
        comment_char="#",  
        has_header=False,
        new_columns=[
            "seqname", "source", "feature", "start", "end", "score",
            "strand", "frame", "attribute"
        ]
    )

    #filtering for exonic regions of GTF file
    exons = gtf.filter(gtf["feature"] == "exon")

    #calculting exon lengths
    exons = exons.with_columns(
        (pl.col("end") - pl.col("start") + 1).alias("exon_length")
    )

    #summing up the lengths
    total_transcriptome_size = exons["exon_length"].sum()

    return total_transcriptome_size


def calculate_total_bases_sequenced(adata):
    #summing all counts in the matrix
    total_counts = adata.X.sum()
    return total_counts

def calculate_sequencing_depth(total_bases, transcriptome_size):
	#dividing to find sequencing depth!
    sequencing_depth = total_bases / transcriptome_size
    return sequencing_depth

#gtf input file
gtf_file = "gencode.v43.annotation.gtf" 

#find transcriptome size
transcriptome_size = calculate_transcriptome_size(gtf_file)
print(f"Transcriptome Size: {transcriptome_size / 1e6:.2f} Mb")

#calculating total bases sequenced
total_bases_sequenced = calculate_total_bases_sequenced(adata)
print(f"Total Bases Sequenced: {total_bases_sequenced:.2f} bp")

#using the above info to find sequencing depth
sequencing_depth = calculate_sequencing_depth(total_bases_sequenced, transcriptome_size)
print(f"Sequencing Depth: {sequencing_depth:.2f}X")

#Next, sparsity:)

import scipy.sparse as sp
#loading in counts matrix
if True:
	matrix=adata.X
	#assuming the matrix is sparse
	total_entries = matrix.shape[0] * matrix.shape[1]
	zero_entries = total_entries - matrix.nnz  #nnz corresponds to the number of nonzero entries
	# else: #using dense method if not sparse!
	# 	total_entries = matrix.size
	# 	zero_entries = (matrix == 0).sum()
	sparsity = zero_entries / total_entries 
	print(f"Sparsity: {sparsity:.2%}") #printing sparsity!

#BATCH effects!!

adata.obs['batch'] = adata.obs['batch'].astype('category')

import seaborn as sns

#with seaborn palette, generating dynamic one and assigning colors to batch variable
batch_categories = adata.obs['batch'].cat.categories
num_batches = len(batch_categories)
batch_colors = sns.color_palette('husl', num_batches) 
adata.uns['batch_colors'] = batch_colors

#contingency table - showing relation of leiden clusters and batch variable 
contingency_table = pd.crosstab(adata.obs['leiden_res0_2'], adata.obs['batch'])
print(contingency_table)

#chi squared statistic performed on this table
chi2, p, dof, expected = chi2_contingency(contingency_table)
print(f"Chi-Square statistic: {chi2}, p-value: {p}")

#calculating batch purity!
cluster_purity = contingency_table.div(contingency_table.sum(axis=1), axis=0).max(axis=1)
print(cluster_purity)
print("Mean Batch Purity Across Clusters:", cluster_purity.mean())

#a plot showing batch purity for each cluster:)
sc.pl.umap(adata, color=['batch', 'leiden_res0_2'],palette=adata.uns['batch_colors'], save='umapLeidenBatchPurity.png')


#Load in anndata object!
adata = ad.read_h5ad('/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/clearedObj2.h5ad')

#First, depth!
def calculate_transcriptome_size(gtf_file):
    #loading GTF file in
    gtf = pl.read_csv(
        gtf_file,
        separator="\t", 
        comment_char="#",  
        has_header=False,
        new_columns=[
            "seqname", "source", "feature", "start", "end", "score",
            "strand", "frame", "attribute"
        ]
    )

    #filtering for exonic regions of GTF file
    exons = gtf.filter(gtf["feature"] == "exon")

    #calculting exon lengths
    exons = exons.with_columns(
        (pl.col("end") - pl.col("start") + 1).alias("exon_length")
    )

    #summing up the lengths
    total_transcriptome_size = exons["exon_length"].sum()

    return total_transcriptome_size


def calculate_total_bases_sequenced(adata):
    #summing all counts in the matrix
    total_counts = adata.X.sum()
    return total_counts

def calculate_sequencing_depth(total_bases, transcriptome_size):
	#dividing to find sequencing depth!
    sequencing_depth = total_bases / transcriptome_size
    return sequencing_depth

#gtf input file
gtf_file = "gencode.v43.annotation.gtf" 

#find transcriptome size
transcriptome_size = calculate_transcriptome_size(gtf_file)
print(f"Transcriptome Size: {transcriptome_size / 1e6:.2f} Mb")

#calculating total bases sequenced
total_bases_sequenced = calculate_total_bases_sequenced(adata)
print(f"Total Bases Sequenced: {total_bases_sequenced:.2f} bp")

#using the above info to find sequencing depth
sequencing_depth = calculate_sequencing_depth(total_bases_sequenced, transcriptome_size)
print(f"Sequencing Depth: {sequencing_depth:.2f}X")

#Next, sparsity:)

#loading in counts matrix
if True:
	matrix=adata.X
	if issparse(matrix): #is the matrix sparse? (probably yes)
		total_entries = matrix.shape[0] * matrix.shape[1]
		zero_entries = total_entries - matrix.nnz  #nnz corresponds to the number of nonzero entries
	else: #using dense method if not sparse!
		total_entries = matrix.size
		zero_entries = (matrix == 0).sum()
	sparsity = zero_entries / total_entries 
	print(f"Sparsity: {sparsity:.2%}") #printing sparsity!

#BATCH effects!!

adata.obs['batch'] = adata.obs['batch'].astype('category')

import seaborn as sns

#with seaborn palette, generating dynamic one and assigning colors to batch variable
batch_categories = adata.obs['batch'].cat.categories
num_batches = len(batch_categories)
batch_colors = sns.color_palette('husl', num_batches) 
adata.uns['batch_colors'] = batch_colors

#contingency table - showing relation of leiden clusters and batch variable 
contingency_table = pd.crosstab(adata.obs['leiden_res0_2'], adata.obs['batch'])
print(contingency_table)

#chi squared statistic performed on this table
chi2, p, dof, expected = chi2_contingency(contingency_table)
print(f"Chi-Square statistic: {chi2}, p-value: {p}")

#calculating batch purity!
cluster_purity = contingency_table.div(contingency_table.sum(axis=1), axis=0).max(axis=1)
print(cluster_purity)
print("Mean Batch Purity Across Clusters:", cluster_purity.mean())

#a plot showing batch purity for each cluster:)
sc.pl.umap(adata, color=['batch', 'leiden_res0_2'],palette=adata.uns['batch_colors'], save='umapLeidenBatchPurity.png')