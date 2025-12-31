import pandas as pd
import anndata as ad
from scipy import sparse
from sklearn.preprocessing import normalize
import numpy as np
import pyreadr

import scanpy as sc
import h5py

# Open the file in read mode
# with h5py.File('/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/postAnnotation110224.h5ad', 'r') as f:
#     # List all top-level keys in the file
#     print(list(f.keys()))
#     # Optionally, inspect subkeys or specific attributes within 'obsm' or other problematic fields
#     print(list(f['obsm'].keys()))

adata = sc.read('/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/postAnnotation110224.h5ad')

# Define the list of patient labels to exclude
excluded_patients = [2518573, 82317494, 83034844, 94430339, 95919181, 52311825, 50101523, 65499271]

# Filter `adata` to exclude these patients before pseudobulking
# adata = adata_old[~adata_old.obs['batch'].isin(excluded_patients)].copy()

obs_to_keep=[]
# Use `adata_filtered` in the rest of the pseudobulking process
# For example:
adata_pb = aggregate_and_filter(adata, adata.obs["cell_type"].cat.categories[0], obs_to_keep=obs_to_keep)

for cell_type in adata.obs["cell_type"].cat.categories[1:]:
    adata_cell_type = aggregate_and_filter(adata, cell_type, obs_to_keep=obs_to_keep)
    adata_pb = adata_pb.concatenate(adata_cell_type)

# Proceed with the rest of your analysis using `adata_pb`

metadata = pd.read_csv('dataset_652_basic_12-23-2021.csv')

file_path_source = '/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/ACE/Final_Pipeline/Batch_Correction/postAnnotation110224.h5ad'
adata_source = ad.read_h5ad(file_path_source)

adata.obs['label'] = [metadata.loc[metadata['projid'] == int(projid_value), 'tot_adverse_exp'].values[0] for projid_value in adata_source.obs['batch']]

adata.obs['replicate'] = adata.obs['batch'][:]

import random  # Needed for shuffling indices

def aggregate_and_filter(
    adata,
    cell_identity,
    donor_key="sample",
    condition_key="label",
    cell_identity_key="cell_type",
    obs_to_keep=[],  # additional metadata to keep, e.g., gender, age, etc.
    replicates_per_patient=1,
):
    # Subset `adata` to the specified cell type
    adata_cell_pop = adata[adata.obs[cell_identity_key] == cell_identity].copy()

    # Check donors with cell counts below the threshold
    size_by_donor = adata_cell_pop.obs.groupby([donor_key]).size()
    donors_to_drop = [donor for donor in size_by_donor.index if size_by_donor[donor] <= NUM_OF_CELL_PER_DONOR]
    if donors_to_drop:
        print("Dropping the following samples:", donors_to_drop)

    # Prepare DataFrame for storing aggregated data
    df = pd.DataFrame(columns=[*adata_cell_pop.var_names, *obs_to_keep])

    adata_cell_pop.obs[donor_key] = adata_cell_pop.obs[donor_key].astype("category")
    for i, donor in enumerate(donors := adata_cell_pop.obs[donor_key].cat.categories):
        print(f"\tProcessing donor {i+1} out of {len(donors)}...", end="\r")
        
        if donor not in donors_to_drop:
            adata_donor = adata_cell_pop[adata_cell_pop.obs[donor_key] == donor]
            indices = list(adata_donor.obs_names)
            random.shuffle(indices)
            indices = np.array_split(np.array(indices), replicates_per_patient)
            
            for i, rep_idx in enumerate(indices):
                adata_replicate = adata_donor[rep_idx]
                agg_dict = {gene: "sum" for gene in adata_replicate.var_names}
                for obs in obs_to_keep:
                    agg_dict[obs] = "first"

                # Aggregate gene expression and metadata
                df_donor = pd.DataFrame(
                    adata_replicate.X.toarray() if sparse.issparse(adata_replicate.X) else adata_replicate.X,
                    index=adata_replicate.obs_names,
                    columns=adata_replicate.var_names,
                )
                df_donor = df_donor.join(adata_replicate.obs[obs_to_keep])
                df_donor = df_donor.groupby(donor_key).agg(agg_dict)
                df_donor[donor_key] = donor
                df.loc[f"donor_{donor}_{i}"] = df_donor.loc[donor]
    print("\n")

    # Create AnnData object from aggregated data
    adata_cell_pop = sc.AnnData(
        df[adata_cell_pop.var_names], obs=df.drop(columns=adata_cell_pop.var_names)
    )
    return adata_cell_pop


import numpy as np
import pandas as pd
# Initializing pseudobulk aggregation
cell_type = adata.obs["cell_type"].cat.categories[0]
adata_pb = aggregate_and_filter(adata, cell_type, obs_to_keep=obs_to_keep)

for i, cell_type in enumerate(adata.obs["cell_type"].cat.categories[1:]):
    adata_cell_type = aggregate_and_filter(adata, cell_type, obs_to_keep=obs_to_keep)
    adata_pb = adata_pb.concatenate(adata_cell_type)

# Proceed with normalization and analysis steps
adata_pb.layers['counts'] = adata_pb.X.copy()
sc.pp.normalize_total(adata_pb, target_sum=1e6)
sc.pp.log1p(adata_pb)
sc.pp.pca(adata_pb)

# Plotting `lib_size` and `log_lib_size` (optional)
lib_size = adata_pb.layers['counts'].sum(axis=1)
adata_pb.obs["lib_size"] = pd.to_numeric(lib_size)
adata_pb.obs["log_lib_size"] = np.log(adata_pb.obs["lib_size"])

for col in ['label', 'cell_type', 'replicate', 'sample', 'lib_size', 'log_lib_size']:
    try:
        sc.pl.pca(adata_pb, color=col, ncols=1, size=300)
        print(f"Column {col} plotted successfully.")
    except Exception as e:
        print(f"Error plotting column {col}: {e}")

# print(
#     f'Processing {cell_type} (1 out of {len(adata.obs["cell_type"].cat.categories)})...'
# )
# adata_pb = aggregate_and_filter(adata, cell_type, obs_to_keep=obs_to_keep)
# for i, cell_type in enumerate(adata.obs["cell_type"].cat.categories[1:]):
#     print(
#         f'Processing {cell_type} ({i+2} out of {len(adata.obs["cell_type"].cat.categories)})...'
#     )
#     adata_cell_type = aggregate_and_filter(adata, cell_type, obs_to_keep=obs_to_keep)
#     adata_pb = adata_pb.concatenate(adata_cell_type)

# adata_pb.layers['counts'] = adata_pb.X.copy()

# sc.pp.normalize_total(adata_pb, target_sum=1e6)
# sc.pp.log1p(adata_pb)
# sc.pp.pca(adata_pb)


# # Step 1: Ensure that counts are a numpy array
# counts_array = np.array(adata_pb.layers["counts"])

# # Step 2: Calculate the library size
# lib_size = np.sum(counts_array, axis=1)

# # Convert lib_size to a numeric type
# lib_size = pd.to_numeric(lib_size)

# # Step 3: Assign the library size to the 'lib_size' column in obs
# adata_pb.obs["lib_size"] = lib_size

# # Ensure that lib_size is now part of adata_pb.obs
# print("lib_size in obs:", adata_pb.obs["lib_size"])

# # Step 4: Calculate the log of the library size
# log_lib_size = np.log(lib_size)

# # Check the type and shape of log_lib_size
# print("log_lib_size type:", type(log_lib_size))
# print("log_lib_size shape:", log_lib_size.shape)

# # Step 5: Assign the log library size to the 'log_lib_size' column in obs
# adata_pb.obs["log_lib_size"] = log_lib_size

# # Ensure that log_lib_size is now part of adata_pb.obs
# print("log_lib_size in obs:", adata_pb.obs["log_lib_size"])

# for col in ['label', 'cell_type', 'replicate', 'sample', 'lib_size', 'log_lib_size']:
#     try:
#         sc.pl.pca(adata_pb, color=col, ncols=1, size=300)
#         print(f"Column {col} plotted successfully.")
#     except Exception as e:
#         print(f"Error plotting column {col}: {e}")

adata_pb.X = adata_pb.layers['counts'].copy()


import anndata
import anndata2ri
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, numpy2ri
from rpy2.robjects.conversion import localconverter

# Activate the automatic conversion between pandas and R data frames
pandas2ri.activate()
numpy2ri.activate()
anndata2ri.activate()

# Subset the AnnData object by cell type
adata_mono = adata_pb[adata_pb.obs["cell_type"] == "Oli"]
adata_mono.obs_names = [
    name.split("_")[2] + "_" + name.split("_")[3] for name in adata_mono.obs_names
]

# Extract data from the AnnData object
X_data = adata_mono.X.toarray() if hasattr(adata_mono.X, "toarray") else adata_mono.X
obs_data = adata_mono.obs
var_data = adata_mono.var

# Convert extracted data to pandas DataFrames
import pandas as pd
X_df = pd.DataFrame(X_data, index=obs_data.index, columns=var_data.index)
obs_df = pd.DataFrame(obs_data)
var_df = pd.DataFrame(var_data)

# Transfer data to the R environment
with localconverter(ro.default_converter + pandas2ri.converter + numpy2ri.converter):
    ro.globalenv["X_df"] = ro.conversion.py2rpy(X_df)
    ro.globalenv["obs_df"] = ro.conversion.py2rpy(obs_df)
    ro.globalenv["var_df"] = ro.conversion.py2rpy(var_df)

# Define the fit_model function in R
fit_model_code = """
fit_model <- function(adata_){
    library(edgeR)
    library(SingleCellExperiment)
    library(rlang)

    # Create edgeR object with counts and grouping factor
    y <- DGEList(assay(adata_, "counts"), group = colData(adata_)$label)

    # Filter genes with low counts
    print("Dimensions before subsetting:")
    print(dim(y))
    keep <- filterByExpr(y)
    y <- y[keep, , keep.lib.sizes=FALSE]
    print("Dimensions after subsetting:")
    print(dim(y))

    # Normalize counts
    y <- calcNormFactors(y)

    # Create vector for contrasts
    group <- paste0(colData(adata_)$label, ".", colData(adata_)$cell_type)
    replicate <- colData(adata_)$replicate

    # Design matrix with multiple donors
    design <- model.matrix(~ 0 + group + replicate)

    # Estimate dispersion
    y <- estimateDisp(y, design = design)

    # Fit model
    fit <- glmQLFit(y, design)

    return(list("fit"=fit, "design"=design, "y"=y))
}
"""
# Pass the function definition to R
ro.r(fit_model_code)

# Create the SingleCellExperiment object in R
create_sce_code = """
library(SingleCellExperiment)
sce <- SingleCellExperiment(assays=list(counts=as.matrix(X_df)))
colData(sce) <- obs_df
rowData(sce) <- var_df
"""
# Pass the R script to R
ro.r(create_sce_code)

# Call the fit_model function in R
fit_model_call_code = """
outs <- fit_model(sce)
"""
# Execute the R function
ro.r(fit_model_call_code)
