# import os

# # Base template for the Python script
# base_script_template ="""
# import pandas as pd
# import os
# import numpy as np
# import seaborn as sns
# import scanpy as sc
# from scipy.stats import median_abs_deviation
# #Figure settings --> can adjust as desired but these were the ones used in SC Best Practices
# sc.settings.verbosity = 0
# sc.settings.set_figure_params(
#     dpi=80,
#     facecolor="white",
#     frameon=False,
# )
# #Figure settings --> can adjust as desired but these were the ones used in SC Best Practices
# sc.settings.verbosity = 0
# sc.settings.set_figure_params(
#     dpi=80,
#     facecolor="white",
#     frameon=False,
# )
# #Get paths to each patient's cellbender output
# library = '{variable}'

# root = os.path.join(os.environ['DATA_ROOT'], 'Data/Tsai/Preprocessing/Preprocessed_Counts/ACE/')+library
# library1 = os.path.join(root, 'processed_feature_bc_matrix_filtered.h5')
# library1_adata = sc.read_10x_h5(library1)
# obj_list=[library1_adata]
# for obj in obj_list:
#     obj.var_names_make_unique()

# for obj in obj_list:
#     # mitochondrial genes
#     obj.var["mt"] = obj.var_names.str.startswith("MT-")
#     # ribosomal genes
#     obj.var["ribo"] = obj.var_names.str.startswith(("RPS", "RPL"))
#     # hemoglobin genes.
#     obj.var["hb"] = obj.var_names.str.contains(("^HB[^(P)]"))

#     #Calculate the qc metrics
#     sc.pp.calculate_qc_metrics(obj, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)

# # Plot the calculated QC Metrics
# for idx, obj in enumerate(obj_list):
#     p1 = sns.displot(obj.obs["total_counts"], bins=100, kde=False)
#     p2 = sc.pl.violin(obj, "pct_counts_mt")
#     p3 = sc.pl.scatter(obj, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

# # Function that takes as input number of MADS that is permissible and metric of interest (column from .obs)
# # and returns whether it is an outlier

# def is_outlier(adata, metric: str, nmads: int):
#     M = adata.obs[metric] #pull column of values of interest (metric we are looking at)
#     outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
#         np.median(M) + nmads * median_abs_deviation(M) < M
#     ) # 'outliers' are are values who are plus or minus nmads (which we specify) from the median
#     return outlier

# # Add outlier status to every cell in every subject
# for obj in obj_list:
#     #outlier status for non mito metrics --> use MAD for this
#     obj.obs["outlier"] = (
#         is_outlier(obj, "log1p_total_counts", 5)
#         | is_outlier(obj, "log1p_n_genes_by_counts", 5)
#         | is_outlier(obj, "pct_counts_in_top_20_genes", 5)
#     )

#     #outlier status for mito metrics
#     #sc best practices included: is_outlier(obj, "pct_counts_mt", 3) | 
#     obj.obs["mt_outlier"] = obj.obs["pct_counts_mt"] > 10 #classify mitochondrial outliers as cells with >10 percent mito

# # Filter the objects based on the constructed thresholds/outlier classifications
# for idx, obj in enumerate(obj_list):
#     print(f"Total number of cells: {obj.n_obs}") #number of cells that the subject initially has
#     obj = obj[(~obj.obs.outlier) & (~obj.obs.mt_outlier)].copy() #subset out cells that are not outliers or mitochondrial outliers
#     print(f"Number of cells after filtering of low quality cells: {obj.n_obs}")
#     p1 = sc.pl.scatter(obj, "total_counts", "n_genes_by_counts", color="pct_counts_mt") # make plots

# for ix, obj in enumerate(obj_list):
#     filename = 'OutputP1.h5ad'
#     obj.write_h5ad(os.path.join(root, filename))
# """

# def generate_scripts(variables, output_dir):
#     os.makedirs(output_dir, exist_ok=True)

#     for variable in variables:
#         script_content = base_script_template.format(variable=variable)
#         script_filename = os.path.join(output_dir, f'script_{variable}.py')

#         with open(script_filename, 'w') as script_file:
#             script_file.write(script_content)

# # List of variables to iterate over
# libraries = os.listdir(os.path.join(os.environ['DATA_ROOT'], 'Data/DeJager/Preprocessed_Counts/'))

# # Directory to save the generated scripts
# output_dir = 'generated_scripts'

# # Generate the scripts
# generate_scripts(libraries, output_dir)

# sc.settings.verbosity = 0
# sc.settings.set_figure_params(
#     dpi=80,
#     facecolor="white",
#     frameon=False,
# )
# #Get paths to each patient's cellbender output
# library = 


#FILTER ALL BAMs script

import os


base_script_templateBAM ="""#!/bin/bash

#SBATCH -n 45                    # Number of cores requested
#SBATCH -t 12:00:00                # Runtime in hours
#SBATCH --mem=800G              # GB memory needed (memory PER CORE)
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-user=nkhera@college.harvard.edu
#SBATCH --mail-type=ALL

source "${CONDA_INIT_SCRIPT}"

echo "Hello, 1!"

cd "${SOCISL_OUTPUT_ROOT}/Tsai"

echo "Hello, 2!"

activate_env "${NEBULA_ENV}"

echo "Hello, 3!"

Rscript nebula_{variable}.Rscript

sbatch nebulaEnrich{variable}.sh

"""


base_script_templateBAMDejager ="""#!/bin/bash

#SBATCH -n 45                    # Number of cores requested
#SBATCH -t 12:00:00                # Runtime in hours
#SBATCH --mem=800G              # GB memory needed (memory PER CORE)
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-user=nkhera@college.harvard.edu
#SBATCH --mail-type=ALL

source "${CONDA_INIT_SCRIPT}"

echo "Hello, 1!"

cd "${SOCISL_OUTPUT_ROOT}/DeJager"

echo "Hello, 2!"

activate_env "${NEBULA_ENV}"

echo "Hello, 3!"

Rscript nebula_{variable}D.Rscript

sbatch nebulaEnrich{variable}.sh

"""


base_script_templateBAMM ="""#!/bin/bash

#SBATCH -n 45                    # Number of cores requested
#SBATCH -t 12:00:00                # Runtime in hours
#SBATCH --mem=800G              # GB memory needed (memory PER CORE)
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-user=nkhera@college.harvard.edu
#SBATCH --mail-type=ALL

source "${CONDA_INIT_SCRIPT}"

echo "Hello, 1!"

cd "${SOCISL_OUTPUT_ROOT}/Tsai"

echo "Hello, 2!"

activate_env "${NEBULA_ENV}"

echo "Hello, 3!"

Rscript male_nebula_{variable}.Rscript
sbatch male_nebulaEnrich{variable}.sh

"""


base_script_templateNebulaRVYay ="""#!/bin/bash
library(edgeR)
library(zellkonverter)
library(SingleCellExperiment)
library(magrittr)
library(data.table) 
library(scran)
library(igraph)
setwd(Sys.getenv("SOCISL_OUTPUT_ROOT"))

celltypes <- c("Ex_L2_3","Ex_L4","Ex_L4_5","Ex_L5_6","Ex_L5","Oli","Ast","Endo","Mic","In_VIP","In_SST","In_PV (Basket)","In_PV (Chandelier)","In_Rosehip","OPC")
sce <- readH5AD("totalAdataAnno012125.h5ad",raw = TRUE)  
ct <- "{variable}"
#loading in SCE object

  #subsetting for cell type
  scek <- sce[, colData(sce)$cell_type == ct]

  #barcode column! 
  colData(scek)$barcode <- colnames(scek) 

  #subsetting patients based on cognitive status
  library(dplyr)

  data <- read.csv("dataset_652_basic_03-23-2022.csv")

  desiredPatientsM <- c(50105301, 10518782, 74284255, 10202345, 15113169, 50101659, 11157783, 10253148, 3713990, 10490993, 
          50106730, 50104134, 11444465, 50405330, 10394182, 50402693, 50405042, 18414513, 44299049, 10101589, 
          10277308, 10502798, 11327005)

  desiredPatientsF <- c(21151608, 20339740, 21180847, 21408652, 20282974, 20248206, 20153010, 20109020, 21406768, 20208992, 
         20634274, 20297403, 50301125, 50401002, 20254902, 50108886, 20970441, 7265221, 66924745, 
         20344143, 60725338, 50402729, 31908032, 69982533, 50301675, 20380417, 20195344, 20929774, 50104008, 
         20504017, 92393245, 32383679, 20254588, 31726180)

  #482428, 18920002

  desired_patientsFR <- data %>%
    filter(projid %in% desiredPatientsF & cogdx == 1)

  desired_patientsFN <- data %>%
    filter(projid %in% desiredPatientsF & cogdx == 4)

  desired_patientsMR <- data %>%
    filter(projid %in% desiredPatientsM & cogdx == 1)

  desired_patientsMN <- data %>%
    filter(projid %in% desiredPatientsM & cogdx == 4)

  colData(scek)$batch <- as.character(colData(scek)$batch)

  desired_patientsFR[[1]] <- as.character(desired_patientsFR[[1]])
  desired_patientsFN[[1]] <- as.character(desired_patientsFN[[1]])
  desired_patientsMR[[1]] <- as.character(desired_patientsMR[[1]])
  desired_patientsMN[[1]] <- as.character(desired_patientsMN[[1]])

  #building female single cell experiment object
  sceFR <- scek[, colData(scek)$batch %in% desired_patientsFR[[1]]]
  sceFN <- scek[, colData(scek)$batch %in% desired_patientsFN[[1]]]


  colData(sceFR)$group <- rep(1, ncol(sceFR))
  colData(sceFN)$group <- rep(0, ncol(sceFN))
  sceF <- cbind(sceFR, sceFN)
  colData(sceF)$projid <- as.character(colData(sceF)$batch)

  #building male single cell experiment object

  sceMR <- scek[, colData(scek)$batch %in% desired_patientsMR[[1]]]
  sceMN <- scek[, colData(scek)$batch %in% desired_patientsMN[[1]]]
  colData(sceMR)$group <- rep(1, ncol(sceMR))
  colData(sceMN)$group <- rep(0, ncol(sceMN))
  data$projid <- as.character(data$projid)

  sceM <- cbind(sceMR, sceMN)
  colData(sceM)$projid <- as.character(colData(sceM)$batch)

  original_colnames <- colnames(sceF)

  meta_df <- as.data.frame(colData(sceF))
  meta_df <- dplyr::left_join(meta_df, data, by = "projid")

  colData(sceF) <- DataFrame(meta_df)
  colnames(sceF) <- original_colnames

  original_colnames <- colnames(sceM)

  meta_df <- as.data.frame(colData(sceM))
  meta_df <- dplyr::left_join(meta_df, data, by = "projid")

  colData(sceM) <- DataFrame(meta_df)
  colnames(sceM) <- original_colnames

  colData(sceF)$age_death <- as.numeric(colData(sceF)$age_death)
  colData(sceF)$pmi <- as.numeric(colData(sceF)$pmi)
  colData(sceF)$educ <- as.numeric(colData(sceF)$educ)
  colData(sceF)$group <- as.numeric(colData(sceF)$group)
  colData(sceF)$msex <- factor(colData(sceF)$msex)
  colData(sceF)$race <- factor(colData(sceF)$race)
  colData(sceF)$apoe <- factor(colData(sceF)$apoe_genotype)
  colData(sceM)$age_death <- as.numeric(colData(sceM)$age_death)
  colData(sceM)$pmi <- as.numeric(colData(sceM)$pmi)
  colData(sceM)$educ <- as.numeric(colData(sceM)$educ)
  colData(sceM)$group <- as.numeric(colData(sceM)$group)
  colData(sceM)$msex <- factor(colData(sceM)$msex)
  colData(sceM)$race <- factor(colData(sceM)$race)
  colData(sceM)$apoe <- factor(colData(sceM)$apoe_genotype)

  # making sure batch and group are factors!
  colData(sceF)$batch <- factor(colData(sceF)$batch)
  colData(sceM)$batch <- factor(colData(sceM)$batch)
  colData(sceF)$group <- factor(colData(sceF)$group, levels = c(0, 1))
  colData(sceM)$group <- factor(colData(sceM)$group, levels = c(0, 1))

  library(nebula)

  #labelling the genes in the raw experiment
  raw_exp <- altExp(sceF, "raw")
  raw_exp <- as(raw_exp, "SingleCellExperiment")
  var_names <- read.csv("var_names.csv", stringsAsFactors = FALSE)
  rownames(raw_exp) <- var_names$gene

  #subsetting for the filtered genes

  common_genes <- intersect(rownames(raw_exp), rownames(sceF))

  assayNames(raw_exp)[assayNames(raw_exp) == "X"] <- "counts"

  #assigning batch and group info
  colData(raw_exp)$batch <- colData(sceF)$batch[match(colnames(raw_exp), rownames(colData(sceF)))]
  colData(raw_exp)$group <- colData(sceF)$group[match(colnames(raw_exp), rownames(colData(sceF)))]
  colData(raw_exp)$pmi <- colData(sceF)$pmi[match(colnames(raw_exp), rownames(colData(sceF)))]
  colData(raw_exp)$educ <- colData(sceF)$educ[match(colnames(raw_exp), rownames(colData(sceF)))]
  colData(raw_exp)$apoe <- colData(sceF)$apoe[match(colnames(raw_exp), rownames(colData(sceF)))]
  colData(raw_exp)$race <- colData(sceF)$race[match(colnames(raw_exp), rownames(colData(sceF)))]
  colData(raw_exp)$age_death <- colData(sceF)$age_death[match(colnames(raw_exp), rownames(colData(sceF)))]
  colData(raw_exp)$age_death <- colData(raw_exp)$age_death / 10
  colData(raw_exp)$pmi <- colData(raw_exp)$pmi / 10

  apoe_odds_map <- c(
    "22" = 0.56,
    "23" = 0.56,
    "24" = 2.64,
    "33" = 1,
    "34" = 3.63,
    "44" = 14.49
  )

  colData(raw_exp)$apoe <- apoe_odds_map[as.character(colData(raw_exp)$apoe)]

  group_metadata <- data.frame(batch = colData(sceF)$batch, group = colData(sceF)$group)

  
  sceE <- raw_exp

  #RUV-III-NB code
  library(variancePartition)
  library(BiocParallel)

  library(NewWave)

  counts_matrix <- counts(sceE)
  keep_feature <- rowSums(counts(sceE) > 0) > ncol(sceE) * 0.1
  sceE <- sceE[keep_feature, ]
  raw_exp <- raw_exp[keep_feature, ]

  colData(raw_exp)$idName <- colData(raw_exp)$batch


  fluidigm_zinb <- newWave(Y=sceE, K = 2, X = "~ group", verbose=TRUE)
  W <- reducedDim(fluidigm_zinb)


  #accounting for ruv factors
  colData(raw_exp) <- cbind(colData(raw_exp), W)


  ruv_formula <- paste("~ group + ", paste(colnames(W), collapse=" + "))
  # model matrix!

  colData(raw_exp)$group <- factor(colData(raw_exp)$group, levels = c("0", "1"))
  design <- model.matrix(as.formula(ruv_formula), data = colData(raw_exp))

  #accounting for ruv factors

  count_df <- as.matrix(assay(raw_exp, "counts"))

  # DGElist + normalization
  dge <- DGEList(counts = count_df)
  dge <- calcNormFactors(dge)  # Normalize counts

  counts_matrix <- counts(raw_exp)
  keep_feature <- rowSums(counts_matrix > 0) > 0
  raw_exp <- raw_exp[keep_feature, ]
  colData(raw_exp)$idName <- colData(raw_exp)$batch
  seuratdata <- scToNeb(obj = raw_exp, assay = "counts", id = "idName", pred = c("group","W1","W2"))

  mat <- seuratdata$count
  nonzero_cols <- Matrix::colSums(mat) > 0
  mat <- mat[, nonzero_cols]
  seuratdata$count <- seuratdata$count[, nonzero_cols]
  seuratdata$pred <- seuratdata$pred[nonzero_cols, ]
  seuratdata$id   <- seuratdata$id[nonzero_cols]

  pred_df <- as.data.frame(seuratdata$pred)
  pred_df[] <- lapply(pred_df, function(x) {
    if (is.character(x)) {
      if (all(grepl("^-?[0-9.]+$", x))) {
        as.numeric(x)
      } else {
        as.factor(x)
      }
    } else {
      x
    }
  })

 
  pred_df$group <- as.character(pred_df$group)
  pred_df$group <- factor(pred_df$group)
  pred_df$group <- relevel(pred_df$group, ref = "0")

  seuratdata$pred <- pred_df

  dge <- DGEList(counts = seuratdata$count)
  dge <- calcNormFactors(dge)
  offset <- dge$samples$lib.size * dge$samples$norm.factors
  stopifnot(length(offset) == ncol(seuratdata$count))
  offset <- pmax(offset, 1e-6)
  seuratdata$offset <- offset


  re2 = nebula(
    count = seuratdata$count, 
    id = seuratdata$id, 
    pred = design, 
    offset = seuratdata$offset,
    output_re=TRUE,
    method ="HL"
  )

  # saving results
  re2_name <- paste0("Fnebula_analysis",ct,"RVNew")
  save(re2, file = paste0("./", re2_name, ".rda"))

  markers <- c("HES4","PDE10A","RPH3A","ST6GAL2","UST")
  print("Nebula Tsai")
  print("fem")
  print(re2$summary[markers,])

"""


base_script_templateNebulaRVYayM ="""#!/bin/bash
library(edgeR)
library(zellkonverter)
library(SingleCellExperiment)
library(magrittr)
library(data.table) 
library(scran)
library(igraph)
setwd(Sys.getenv("SOCISL_OUTPUT_ROOT"))

celltypes <- c("Ex_L2_3","Ex_L4","Ex_L4_5","Ex_L5_6","Ex_L5","Oli","Ast","Endo","Mic","In_VIP","In_SST","In_PV (Basket)","In_PV (Chandelier)","In_Rosehip","OPC")
sce <- readH5AD("totalAdataAnno012125.h5ad",raw = TRUE)  
ct <- "{variable}"
#loading in SCE object

  #subsetting for cell type
  scek <- sce[, colData(sce)$cell_type == ct]

  #barcode column! 
  colData(scek)$barcode <- colnames(scek) 

  #subsetting patients based on cognitive status
  library(dplyr)

  data <- read.csv("dataset_652_basic_03-23-2022.csv")

  desiredPatientsM <- c(50105301, 10518782, 74284255, 10202345, 15113169, 50101659, 11157783, 10253148, 3713990, 10490993, 
          50106730, 50104134, 11444465, 50405330, 10394182, 50402693, 50405042, 18414513, 44299049, 10101589, 
          10277308, 10502798, 11327005)

  desiredPatientsF <- c(21151608, 20339740, 21180847, 21408652, 20282974, 20248206, 20153010, 20109020, 21406768, 20208992, 
         20634274, 20297403, 50301125, 50401002, 20254902, 50108886, 20970441, 7265221, 66924745, 
         20344143, 60725338, 50402729, 31908032, 69982533, 50301675, 20380417, 20195344, 20929774, 50104008, 
         20504017, 92393245, 32383679, 20254588, 31726180)

  #482428, 18920002

  desired_patientsFR <- data %>%
    filter(projid %in% desiredPatientsF & cogdx == 1)

  desired_patientsFN <- data %>%
    filter(projid %in% desiredPatientsF & cogdx == 4)

  desired_patientsMR <- data %>%
    filter(projid %in% desiredPatientsM & cogdx == 1)

  desired_patientsMN <- data %>%
    filter(projid %in% desiredPatientsM & cogdx == 4)

  colData(scek)$batch <- as.character(colData(scek)$batch)

  desired_patientsFR[[1]] <- as.character(desired_patientsFR[[1]])
  desired_patientsFN[[1]] <- as.character(desired_patientsFN[[1]])
  desired_patientsMR[[1]] <- as.character(desired_patientsMR[[1]])
  desired_patientsMN[[1]] <- as.character(desired_patientsMN[[1]])

  #building female single cell experiment object
  sceFR <- scek[, colData(scek)$batch %in% desired_patientsFR[[1]]]
  sceFN <- scek[, colData(scek)$batch %in% desired_patientsFN[[1]]]


  colData(sceFR)$group <- rep(1, ncol(sceFR))
  colData(sceFN)$group <- rep(0, ncol(sceFN))
  sceF <- cbind(sceFR, sceFN)
  colData(sceF)$projid <- as.character(colData(sceF)$batch)

  #building male single cell experiment object

  sceMR <- scek[, colData(scek)$batch %in% desired_patientsMR[[1]]]
  sceMN <- scek[, colData(scek)$batch %in% desired_patientsMN[[1]]]
  colData(sceMR)$group <- rep(1, ncol(sceMR))
  colData(sceMN)$group <- rep(0, ncol(sceMN))
  data$projid <- as.character(data$projid)

  sceM <- cbind(sceMR, sceMN)
  colData(sceM)$projid <- as.character(colData(sceM)$batch)

  original_colnames <- colnames(sceF)

  meta_df <- as.data.frame(colData(sceF))
  meta_df <- dplyr::left_join(meta_df, data, by = "projid")

  colData(sceF) <- DataFrame(meta_df)
  colnames(sceF) <- original_colnames

  original_colnames <- colnames(sceM)

  meta_df <- as.data.frame(colData(sceM))
  meta_df <- dplyr::left_join(meta_df, data, by = "projid")

  colData(sceM) <- DataFrame(meta_df)
  colnames(sceM) <- original_colnames

  colData(sceF)$age_death <- as.numeric(colData(sceF)$age_death)
  colData(sceF)$pmi <- as.numeric(colData(sceF)$pmi)
  colData(sceF)$educ <- as.numeric(colData(sceF)$educ)
  colData(sceF)$group <- as.numeric(colData(sceF)$group)
  colData(sceF)$msex <- factor(colData(sceF)$msex)
  colData(sceF)$race <- factor(colData(sceF)$race)
  colData(sceF)$apoe <- factor(colData(sceF)$apoe_genotype)
  colData(sceM)$age_death <- as.numeric(colData(sceM)$age_death)
  colData(sceM)$pmi <- as.numeric(colData(sceM)$pmi)
  colData(sceM)$educ <- as.numeric(colData(sceM)$educ)
  colData(sceM)$group <- as.numeric(colData(sceM)$group)
  colData(sceM)$msex <- factor(colData(sceM)$msex)
  colData(sceM)$race <- factor(colData(sceM)$race)
  colData(sceM)$apoe <- factor(colData(sceM)$apoe_genotype)

  # making sure batch and group are factors!
  colData(sceF)$batch <- factor(colData(sceF)$batch)
  colData(sceM)$batch <- factor(colData(sceM)$batch)
  colData(sceF)$group <- factor(colData(sceF)$group, levels = c(0, 1))
  colData(sceM)$group <- factor(colData(sceM)$group, levels = c(0, 1))

  library(nebula)

  #labelling the genes in the raw experiment
  raw_exp <- altExp(sceM, "raw")
  raw_exp <- as(raw_exp, "SingleCellExperiment")
  var_names <- read.csv("var_names.csv", stringsAsFactors = FALSE)
  rownames(raw_exp) <- var_names$gene

  #subsetting for the filtered genes

  common_genes <- intersect(rownames(raw_exp), rownames(sceM))

  assayNames(raw_exp)[assayNames(raw_exp) == "X"] <- "counts"

  #assigning batch and group info
  colData(raw_exp)$batch <- colData(sceM)$batch[match(colnames(raw_exp), rownames(colData(sceM)))]
  colData(raw_exp)$group <- colData(sceM)$group[match(colnames(raw_exp), rownames(colData(sceM)))]
  colData(raw_exp)$pmi <- colData(sceM)$pmi[match(colnames(raw_exp), rownames(colData(sceM)))]
  colData(raw_exp)$educ <- colData(sceM)$educ[match(colnames(raw_exp), rownames(colData(sceM)))]
  colData(raw_exp)$apoe <- colData(sceM)$apoe[match(colnames(raw_exp), rownames(colData(sceM)))]
  colData(raw_exp)$race <- colData(sceM)$race[match(colnames(raw_exp), rownames(colData(sceM)))]
  colData(raw_exp)$age_death <- colData(sceM)$age_death[match(colnames(raw_exp), rownames(colData(sceM)))]
  colData(raw_exp)$age_death <- colData(raw_exp)$age_death / 10
  colData(raw_exp)$pmi <- colData(raw_exp)$pmi / 10

  apoe_odds_map <- c(
    "22" = 0.56,
    "23" = 0.56,
    "24" = 2.64,
    "33" = 1,
    "34" = 3.63,
    "44" = 14.49
  )

  colData(raw_exp)$apoe <- apoe_odds_map[as.character(colData(raw_exp)$apoe)]

  group_metadata <- data.frame(batch = colData(sceM)$batch, group = colData(sceM)$group)

  
  sceE <- raw_exp

  #RUV-III-NB code
  library(variancePartition)
  library(BiocParallel)

  library(NewWave)

  counts_matrix <- counts(sceE)
  keep_feature <- rowSums(counts(sceE) > 0) > ncol(sceE) * 0.1
  sceE <- sceE[keep_feature, ]
  raw_exp <- raw_exp[keep_feature, ]

  colData(raw_exp)$idName <- colData(raw_exp)$batch


  fluidigm_zinb <- newWave(Y=sceE, K = 2, X = "~ group", verbose=TRUE)
  W <- reducedDim(fluidigm_zinb)


  #accounting for ruv factors
  colData(raw_exp) <- cbind(colData(raw_exp), W)


  ruv_formula <- paste("~ group + ", paste(colnames(W), collapse=" + "))
  # model matrix!

  colData(raw_exp)$group <- factor(colData(raw_exp)$group, levels = c("0", "1"))
  design <- model.matrix(as.formula(ruv_formula), data = colData(raw_exp))

  #accounting for ruv factors

  count_df <- as.matrix(assay(raw_exp, "counts"))

  # DGElist + normalization
  dge <- DGEList(counts = count_df)
  dge <- calcNormFactors(dge)  # Normalize counts

  counts_matrix <- counts(raw_exp)
  keep_feature <- rowSums(counts_matrix > 0) > 0
  raw_exp <- raw_exp[keep_feature, ]
  colData(raw_exp)$idName <- colData(raw_exp)$batch
  seuratdata <- scToNeb(obj = raw_exp, assay = "counts", id = "idName", pred = c("group","W1","W2"))

  mat <- seuratdata$count
  nonzero_cols <- Matrix::colSums(mat) > 0
  mat <- mat[, nonzero_cols]
  seuratdata$count <- seuratdata$count[, nonzero_cols]
  seuratdata$pred <- seuratdata$pred[nonzero_cols, ]
  seuratdata$id   <- seuratdata$id[nonzero_cols]

  pred_df <- as.data.frame(seuratdata$pred)
  pred_df[] <- lapply(pred_df, function(x) {
    if (is.character(x)) {
      if (all(grepl("^-?[0-9.]+$", x))) {
        as.numeric(x)
      } else {
        as.factor(x)
      }
    } else {
      x
    }
  })

 
  pred_df$group <- as.character(pred_df$group)
  pred_df$group <- factor(pred_df$group)
  pred_df$group <- relevel(pred_df$group, ref = "0")

  seuratdata$pred <- pred_df

  dge <- DGEList(counts = seuratdata$count)
  dge <- calcNormFactors(dge)
  offset <- dge$samples$lib.size * dge$samples$norm.factors
  stopifnot(length(offset) == ncol(seuratdata$count))
  offset <- pmax(offset, 1e-6)
  seuratdata$offset <- offset

  design <- model.matrix(~ group + W1 + W2, data = as.data.frame(seuratdata$pred))

  re2 = nebula(
    count = seuratdata$count, 
    id = seuratdata$id, 
    pred = design, 
    offset = seuratdata$offset,
    output_re=TRUE,
    method ="HL"
  )

  # saving results
  re2_name <- paste0("Mnebula_analysis",ct,"RVNew")
  save(re2, file = paste0("./", re2_name, ".rda"))

  markers <- c("HES4","PDE10A","RPH3A","ST6GAL2","UST")
  print("Nebula Tsai")
  print("fem")
  print(re2$summary[markers,])

"""



base_script_templateNebulaRVYayD ="""

library(edgeR)
library(zellkonverter)
library(SingleCellExperiment)
library(magrittr)
library(data.table) 
library(scran)
library(igraph)

setwd(Sys.getenv("SOCISL_OUTPUT_ROOT"))

celltypes <- c("Ex_L2_3","Ex_L4","Ex_L4_5","Ex_L5_6","Ex_L5","Oli","Ast","Endo","Mic","In_VIP","In_SST","In_PV (Basket)","In_PV (Chandelier)","In_Rosehip","OPC")
sce <- readH5AD(file.path(Sys.getenv("SOCISL_OUTPUT_ROOT"), "totalAdataAnno040825.h5ad",raw = TRUE)  
ct <- "{variable}"

  #loading in SCE object

  #subsetting for cell type
  scek <- sce[, colData(sce)$cell_type == ct]

  #barcode column! 
  colData(scek)$barcode <- colnames(scek) 

  #subsetting patients based on cognitive status
  library(dplyr)

  data <- read.csv("dataset_652_basic_03-23-2022.csv")

  desiredPatientsM <- c(10205244, 15218541, 15176592, 11165535, 10246987, 11353057, 10277308,10394182, 50405330, 10202345, 11200645, 10490993, 50101659, 70883578, 50405042, 10218339, 10101327, 10101589, 10248033, 10684501, 11157783, 11326252, 11327005,
  11606935, 18414513, 21362537, 22789958, 43485807, 44019405, 66754397,90267190, 50402693)

  desiredPatientsF <- c(21000180, 20252720, 20933324, 21176459, 20173942, 21293107, 98953007, 20358955, 20500815, 21000630, 81852640, 20242958, 20152393, 54122640, 14641578, 81874628, 20898476, 89614402, 60725338, 32383679, 20983400, 21157370, 20153984, 20236838,21411459, 63874408, 50106578, 50108886, 20221273, 20254588, 21140119, 61827429,75675336, 69432088)

  colData(scek)$batch <- as.character(colData(scek)$patient_id)
  colData(scek)$projid <- as.character(colData(scek)$patient_id)
  colData(scek)$projid <- as.factor(colData(scek)$projid)
  data$projid <- as.character(data$projid)

  desired_patientsFR <- data %>%
    filter(projid %in% desiredPatientsF & cogdx == 1)

  desired_patientsFN <- data %>%
    filter(projid %in% desiredPatientsF & cogdx == 4)

  desired_patientsMR <- data %>%
    filter(projid %in% desiredPatientsM & cogdx == 1)

  desired_patientsMN <- data %>%
    filter(projid %in% desiredPatientsM & cogdx == 4)

  colData(scek)$batch <- as.character(colData(scek)$batch)

  desired_patientsFR[[1]] <- as.character(desired_patientsFR[[1]])
  desired_patientsFN[[1]] <- as.character(desired_patientsFN[[1]])
  desired_patientsMR[[1]] <- as.character(desired_patientsMR[[1]])
  desired_patientsMN[[1]] <- as.character(desired_patientsMN[[1]])

  #building female single cell experiment object
  sceFR <- scek[, colData(scek)$batch %in% desired_patientsFR[[1]]]
  sceFN <- scek[, colData(scek)$batch %in% desired_patientsFN[[1]]]


  colData(sceFR)$group <- rep(1, ncol(sceFR))
  colData(sceFN)$group <- rep(0, ncol(sceFN))
  sceF <- cbind(sceFR, sceFN)
  colData(sceF)$projid <- as.character(colData(sceF)$batch)

  #building male single cell experiment object

  sceMR <- scek[, colData(scek)$batch %in% desired_patientsMR[[1]]]
  sceMN <- scek[, colData(scek)$batch %in% desired_patientsMN[[1]]]
  colData(sceMR)$group <- rep(1, ncol(sceMR))
  colData(sceMN)$group <- rep(0, ncol(sceMN))
  data$projid <- as.character(data$projid)

  sceM <- cbind(sceMR, sceMN)
  colData(sceM)$projid <- as.character(colData(sceM)$batch)

  original_colnames <- colnames(sceF)

  meta_df <- as.data.frame(colData(sceF))
  meta_df <- dplyr::left_join(meta_df, data, by = "projid")

  colData(sceF) <- DataFrame(meta_df)
  colnames(sceF) <- original_colnames

  original_colnames <- colnames(sceM)

  meta_df <- as.data.frame(colData(sceM))
  meta_df <- dplyr::left_join(meta_df, data, by = "projid")

  colData(sceM) <- DataFrame(meta_df)
  colnames(sceM) <- original_colnames

  colData(sceF)$age_death <- as.numeric(colData(sceF)$age_death)
  colData(sceF)$pmi <- as.numeric(colData(sceF)$pmi)
  colData(sceF)$educ <- as.numeric(colData(sceF)$educ)
  colData(sceF)$group <- as.numeric(colData(sceF)$group)
  colData(sceF)$msex <- factor(colData(sceF)$msex)
  colData(sceF)$race <- factor(colData(sceF)$race)
  colData(sceF)$apoe <- factor(colData(sceF)$apoe_genotype)
  colData(sceM)$age_death <- as.numeric(colData(sceM)$age_death)
  colData(sceM)$pmi <- as.numeric(colData(sceM)$pmi)
  colData(sceM)$educ <- as.numeric(colData(sceM)$educ)
  colData(sceM)$group <- as.numeric(colData(sceM)$group)
  colData(sceM)$msex <- factor(colData(sceM)$msex)
  colData(sceM)$race <- factor(colData(sceM)$race)
  colData(sceM)$apoe <- factor(colData(sceM)$apoe_genotype)

  # making sure batch and group are factors!
  colData(sceF)$batch <- factor(colData(sceF)$batch)
  colData(sceM)$batch <- factor(colData(sceM)$batch)
  colData(sceF)$group <- factor(colData(sceF)$group, levels = c(0, 1))
  colData(sceM)$group <- factor(colData(sceM)$group, levels = c(0, 1))

  library(nebula)

  #labelling the genes in the raw experiment
  raw_exp <- altExp(sceF, "raw")
  raw_exp <- as(raw_exp, "SingleCellExperiment")
  var_names <- read.csv("var_names.csv", stringsAsFactors = FALSE)
  rownames(raw_exp) <- var_names$gene

  #subsetting for the filtered genes

  common_genes <- intersect(rownames(raw_exp), rownames(sceF))

  assayNames(raw_exp)[assayNames(raw_exp) == "X"] <- "counts"

  #assigning batch and group info
  colData(raw_exp)$batch <- colData(sceF)$batch[match(colnames(raw_exp), rownames(colData(sceF)))]
  colData(raw_exp)$group <- colData(sceF)$group[match(colnames(raw_exp), rownames(colData(sceF)))]
  colData(raw_exp)$pmi <- colData(sceF)$pmi[match(colnames(raw_exp), rownames(colData(sceF)))]
  colData(raw_exp)$educ <- colData(sceF)$educ[match(colnames(raw_exp), rownames(colData(sceF)))]
  colData(raw_exp)$apoe <- colData(sceF)$apoe[match(colnames(raw_exp), rownames(colData(sceF)))]
  colData(raw_exp)$race <- colData(sceF)$race[match(colnames(raw_exp), rownames(colData(sceM)))]
  colData(raw_exp)$age_death <- colData(sceF)$age_death[match(colnames(raw_exp), rownames(colData(sceF)))]
  colData(raw_exp)$age_death <- colData(raw_exp)$age_death / 10
  colData(raw_exp)$pmi <- colData(raw_exp)$pmi / 10

  apoe_odds_map <- c(
    "22" = 0.56,
    "23" = 0.56,
    "24" = 2.64,
    "33" = 1,
    "34" = 3.63,
    "44" = 14.49
  )

  colData(raw_exp)$apoe <- apoe_odds_map[as.character(colData(raw_exp)$apoe)]

  group_metadata <- data.frame(batch = colData(sceF)$batch, group = colData(sceF)$group)

  
  sceE <- raw_exp

  #RUV-III-NB code
  library(variancePartition)
  library(BiocParallel)

  library(NewWave)

  counts_matrix <- counts(sceE)
  keep_feature <- rowSums(counts(sceE) > 0) > ncol(sceE) * 0.25
  sceE <- sceE[keep_feature, ]
  raw_exp <- raw_exp[keep_feature, ]

  colData(raw_exp)$idName <- colData(raw_exp)$batch

  fluidigm_zinb <- newWave(Y=sceE, K = 3, X = "~ group", verbose=TRUE)
  W <- reducedDim(fluidigm_zinb)


  #accounting for ruv factors
  colData(raw_exp) <- cbind(colData(raw_exp), W)


  ruv_formula <- paste("~ group + ", paste(colnames(W), collapse=" + "))
  # model matrix!

  colData(raw_exp)$group <- factor(colData(raw_exp)$group, levels = c("0", "1"))
  design <- model.matrix(as.formula(ruv_formula), data = colData(raw_exp))

  #accounting for ruv factors

  count_df <- as.matrix(assay(raw_exp, "counts"))

  # DGElist + normalization
  dge <- DGEList(counts = count_df)
  dge <- calcNormFactors(dge)  # Normalize counts

  counts_matrix <- counts(raw_exp)
  keep_feature <- rowSums(counts_matrix > 0) > 0
  raw_exp <- raw_exp[keep_feature, ]
  colData(raw_exp)$idName <- colData(raw_exp)$batch
  seuratdata <- scToNeb(obj = raw_exp, assay = "counts", id = "idName", pred = c("group","W1","W2","W3"))

  mat <- seuratdata$count
  nonzero_cols <- Matrix::colSums(mat) > 0
  mat <- mat[, nonzero_cols]
  seuratdata$count <- seuratdata$count[, nonzero_cols]
  seuratdata$pred <- seuratdata$pred[nonzero_cols, ]
  seuratdata$id   <- seuratdata$id[nonzero_cols]

  pred_df <- as.data.frame(seuratdata$pred)
  pred_df[] <- lapply(pred_df, function(x) {
    if (is.character(x)) {
      if (all(grepl("^-?[0-9.]+$", x))) {
        as.numeric(x)
      } else {
        as.factor(x)
      }
    } else {
      x
    }
  })

  
 
  pred_df$group <- as.character(pred_df$group)
  pred_df$group <- factor(pred_df$group)
  pred_df$group <- relevel(pred_df$group, ref = "0")

  seuratdata$pred <- pred_df

  dge <- DGEList(counts = seuratdata$count)
  dge <- calcNormFactors(dge)
  offset <- dge$samples$lib.size * dge$samples$norm.factors
  stopifnot(length(offset) == ncol(seuratdata$count))
  offset <- pmax(offset, 1e-6)
  seuratdata$offset <- offset

 
  seuratdata$pred$group <- as.character(seuratdata$pred$group)
  seuratdata$pred$group <- factor(seuratdata$pred$group)
  seuratdata$pred$group <- relevel(seuratdata$pred$group, ref = "0")

  design <- model.matrix(~ group + W1 + W2 + W3, data = as.data.frame(seuratdata$pred))

  seuratdata=group_cell(count=seuratdata$count,id=seuratdata$id,pred=design)

  re2 = nebula(
    count = seuratdata$count, 
    id = seuratdata$id, 
    pred = design, 
    offset = seuratdata$offset,
    output_re=TRUE,
    method ="HL"
  )

  # saving results
  re2_name <- paste0("Fdejnebula_analysis",ct,"RVN")
  save(re2, file = paste0("./", re2_name, ".rda"))

  markers <- c("HES4","PDE10A","RPH3A","ST6GAL2","UST")
  print("Nebula Tsai")
  print("fem")
  print(re2$summary[markers,])

"""

base_scriptTemplateMDejN = """


library(edgeR)
library(zellkonverter)
library(SingleCellExperiment)
library(magrittr)
library(data.table) 
library(scran)
library(igraph)

setwd(Sys.getenv("SOCISL_OUTPUT_ROOT"))

celltypes <- c("Ex_L2_3","Ex_L4","Ex_L4_5","Ex_L5_6","Ex_L5","Oli","Ast","Endo","Mic","In_VIP","In_SST","In_PV (Basket)","In_PV (Chandelier)","In_Rosehip","OPC")
sce <- readH5AD(file.path(Sys.getenv("SOCISL_OUTPUT_ROOT"), "totalAdataAnno040825.h5ad",raw = TRUE)  
ct <- "{variable}"

  #loading in SCE object

  #subsetting for cell type
  scek <- sce[, colData(sce)$cell_type == ct]

  #barcode column! 
  colData(scek)$barcode <- colnames(scek) 

  #subsetting patients based on cognitive status
  library(dplyr)

  data <- read.csv("dataset_652_basic_03-23-2022.csv")

  desiredPatientsM <- c(10205244, 15218541, 15176592, 11165535, 10246987, 11353057, 10277308,10394182, 50405330, 10202345, 11200645, 10490993, 50101659, 70883578, 50405042, 10218339, 10101327, 10101589, 10248033, 10684501, 11157783, 11326252, 11327005,
  11606935, 18414513, 21362537, 22789958, 43485807, 44019405, 66754397,90267190, 50402693)

  desiredPatientsF <- c(21000180, 20252720, 20933324, 21176459, 20173942, 21293107, 98953007, 20358955, 20500815, 21000630, 81852640, 20242958, 20152393, 54122640, 14641578, 81874628, 20898476, 89614402, 60725338, 32383679, 20983400, 21157370, 20153984, 20236838,21411459, 63874408, 50106578, 50108886, 20221273, 20254588, 21140119, 61827429,75675336, 69432088)

  colData(scek)$batch <- as.character(colData(scek)$patient_id)
  colData(scek)$projid <- as.character(colData(scek)$patient_id)
  colData(scek)$projid <- as.factor(colData(scek)$projid)
  data$projid <- as.character(data$projid)

  desired_patientsFR <- data %>%
    filter(projid %in% desiredPatientsF & cogdx == 1)

  desired_patientsFN <- data %>%
    filter(projid %in% desiredPatientsF & cogdx == 4)

  desired_patientsMR <- data %>%
    filter(projid %in% desiredPatientsM & cogdx == 1)

  desired_patientsMN <- data %>%
    filter(projid %in% desiredPatientsM & cogdx == 4)

  colData(scek)$batch <- as.character(colData(scek)$batch)

  desired_patientsFR[[1]] <- as.character(desired_patientsFR[[1]])
  desired_patientsFN[[1]] <- as.character(desired_patientsFN[[1]])
  desired_patientsMR[[1]] <- as.character(desired_patientsMR[[1]])
  desired_patientsMN[[1]] <- as.character(desired_patientsMN[[1]])

  #building female single cell experiment object
  sceFR <- scek[, colData(scek)$batch %in% desired_patientsFR[[1]]]
  sceFN <- scek[, colData(scek)$batch %in% desired_patientsFN[[1]]]


  colData(sceFR)$group <- rep(1, ncol(sceFR))
  colData(sceFN)$group <- rep(0, ncol(sceFN))
  sceF <- cbind(sceFR, sceFN)
  colData(sceF)$projid <- as.character(colData(sceF)$batch)

  #building male single cell experiment object

  sceMR <- scek[, colData(scek)$batch %in% desired_patientsMR[[1]]]
  sceMN <- scek[, colData(scek)$batch %in% desired_patientsMN[[1]]]
  colData(sceMR)$group <- rep(1, ncol(sceMR))
  colData(sceMN)$group <- rep(0, ncol(sceMN))
  data$projid <- as.character(data$projid)

  sceM <- cbind(sceMR, sceMN)
  colData(sceM)$projid <- as.character(colData(sceM)$batch)

  original_colnames <- colnames(sceF)

  meta_df <- as.data.frame(colData(sceF))
  meta_df <- dplyr::left_join(meta_df, data, by = "projid")

  colData(sceF) <- DataFrame(meta_df)
  colnames(sceF) <- original_colnames

  original_colnames <- colnames(sceM)

  meta_df <- as.data.frame(colData(sceM))
  meta_df <- dplyr::left_join(meta_df, data, by = "projid")

  colData(sceM) <- DataFrame(meta_df)
  colnames(sceM) <- original_colnames

  colData(sceF)$age_death <- as.numeric(colData(sceF)$age_death)
  colData(sceF)$pmi <- as.numeric(colData(sceF)$pmi)
  colData(sceF)$educ <- as.numeric(colData(sceF)$educ)
  colData(sceF)$group <- as.numeric(colData(sceF)$group)
  colData(sceF)$msex <- factor(colData(sceF)$msex)
  colData(sceF)$race <- factor(colData(sceF)$race)
  colData(sceF)$apoe <- factor(colData(sceF)$apoe_genotype)
  colData(sceM)$age_death <- as.numeric(colData(sceM)$age_death)
  colData(sceM)$pmi <- as.numeric(colData(sceM)$pmi)
  colData(sceM)$educ <- as.numeric(colData(sceM)$educ)
  colData(sceM)$group <- as.numeric(colData(sceM)$group)
  colData(sceM)$msex <- factor(colData(sceM)$msex)
  colData(sceM)$race <- factor(colData(sceM)$race)
  colData(sceM)$apoe <- factor(colData(sceM)$apoe_genotype)

  # making sure batch and group are factors!
  colData(sceF)$batch <- factor(colData(sceF)$batch)
  colData(sceM)$batch <- factor(colData(sceM)$batch)
  colData(sceF)$group <- factor(colData(sceF)$group, levels = c(0, 1))
  colData(sceM)$group <- factor(colData(sceM)$group, levels = c(0, 1))

  library(nebula)

  #labelling the genes in the raw experiment
  raw_exp <- altExp(sceM, "raw")
  raw_exp <- as(raw_exp, "SingleCellExperiment")
  var_names <- read.csv("var_names.csv", stringsAsFactors = FALSE)
  rownames(raw_exp) <- var_names$gene

  #subsetting for the filtered genes

  common_genes <- intersect(rownames(raw_exp), rownames(sceM))

  assayNames(raw_exp)[assayNames(raw_exp) == "X"] <- "counts"

  #assigning batch and group info
  colData(raw_exp)$batch <- colData(sceM)$batch[match(colnames(raw_exp), rownames(colData(sceM)))]
  colData(raw_exp)$group <- colData(sceM)$group[match(colnames(raw_exp), rownames(colData(sceM)))]
  colData(raw_exp)$pmi <- colData(sceM)$pmi[match(colnames(raw_exp), rownames(colData(sceM)))]
  colData(raw_exp)$educ <- colData(sceM)$educ[match(colnames(raw_exp), rownames(colData(sceM)))]
  colData(raw_exp)$apoe <- colData(sceM)$apoe[match(colnames(raw_exp), rownames(colData(sceM)))]
  colData(raw_exp)$race <- colData(sceM)$race[match(colnames(raw_exp), rownames(colData(sceM)))]
  colData(raw_exp)$age_death <- colData(sceM)$age_death[match(colnames(raw_exp), rownames(colData(sceM)))]
  colData(raw_exp)$age_death <- colData(raw_exp)$age_death / 10
  colData(raw_exp)$pmi <- colData(raw_exp)$pmi / 10

  apoe_odds_map <- c(
    "22" = 0.56,
    "23" = 0.56,
    "24" = 2.64,
    "33" = 1,
    "34" = 3.63,
    "44" = 14.49
  )

  colData(raw_exp)$apoe <- apoe_odds_map[as.character(colData(raw_exp)$apoe)]

  group_metadata <- data.frame(batch = colData(sceM)$batch, group = colData(sceM)$group)

  
  sceE <- raw_exp

  #RUV-III-NB code
  library(variancePartition)
  library(BiocParallel)

  library(NewWave)

  counts_matrix <- counts(sceE)
  keep_feature <- rowSums(counts(sceE) > 0) > ncol(sceE) * 0.1
  sceE <- sceE[keep_feature, ]
  raw_exp <- raw_exp[keep_feature, ]

  colData(raw_exp)$idName <- colData(raw_exp)$batch

  fluidigm_zinb <- newWave(Y=sceE, K = 2, X = "~ group", verbose=TRUE)
  W <- reducedDim(fluidigm_zinb)


  #accounting for ruv factors
  colData(raw_exp) <- cbind(colData(raw_exp), W)


  ruv_formula <- paste("~ group + ", paste(colnames(W), collapse=" + "))
  # model matrix!

  colData(raw_exp)$group <- factor(colData(raw_exp)$group, levels = c("0", "1"))
  design <- model.matrix(as.formula(ruv_formula), data = colData(raw_exp))

  #accounting for ruv factors

  count_df <- as.matrix(assay(raw_exp, "counts"))

  # DGElist + normalization
  dge <- DGEList(counts = count_df)
  dge <- calcNormFactors(dge)  # Normalize counts

  counts_matrix <- counts(raw_exp)
  keep_feature <- rowSums(counts_matrix > 0) > 0
  raw_exp <- raw_exp[keep_feature, ]
  colData(raw_exp)$idName <- colData(raw_exp)$batch
  print(colData(raw_exp))
  seuratdata <- scToNeb(obj = raw_exp, assay = "counts", id = "idName", pred = c("group","W1","W2"))

  mat <- seuratdata$count
  nonzero_cols <- Matrix::colSums(mat) > 0
  mat <- mat[, nonzero_cols]
  seuratdata$count <- seuratdata$count[, nonzero_cols]
  seuratdata$pred <- seuratdata$pred[nonzero_cols, ]
  seuratdata$id   <- seuratdata$id[nonzero_cols]

  pred_df <- as.data.frame(seuratdata$pred)
  pred_df[] <- lapply(pred_df, function(x) {
    if (is.character(x)) {
      if (all(grepl("^-?[0-9.]+$", x))) {
        as.numeric(x)
      } else {
        as.factor(x)
      }
    } else {
      x
    }
  })

  
 
  pred_df$group <- as.character(pred_df$group)
  pred_df$group <- factor(pred_df$group)
  pred_df$group <- relevel(pred_df$group, ref = "0")

  seuratdata$pred <- pred_df

  dge <- DGEList(counts = seuratdata$count)
  dge <- calcNormFactors(dge)
  offset <- dge$samples$lib.size * dge$samples$norm.factors
  stopifnot(length(offset) == ncol(seuratdata$count))
  offset <- pmax(offset, 1e-6)
  seuratdata$offset <- offset

  design <- model.matrix(~ group + W1 + W2, data = as.data.frame(seuratdata$pred))

  seuratdata=group_cell(count=seuratdata$count,id=seuratdata$id,pred=design)

  pred_df <- as.data.frame(seuratdata$pred)
 
  pred_df$group <- as.character(pred_df$group)
  pred_df$group <- factor(pred_df$group)
  pred_df$group <- relevel(pred_df$group, ref = "0")

  design2 <- model.matrix(~ group1 + W1 + W2, data = as.data.frame(pred_df))

  re2 = nebula(
    count = seuratdata$count, 
    id = seuratdata$id, 
    pred = design2, 
    offset = seuratdata$offset,
    output_re=TRUE,
    method ="HL"
  )

  # saving results
  re2_name <- paste0("dejnebula_analysis",ct,"RVNew")
  save(re2, file = paste0("./", re2_name, ".rda"))

  markers <- c("HES4","PDE10A","RPH3A","ST6GAL2","UST")
  print("Nebula Tsai")
  print("male")
  print(re2$summary[markers,])

  """



base_scriptTemplateFDejN = """


library(edgeR)
library(zellkonverter)
library(SingleCellExperiment)
library(magrittr)
library(data.table) 
library(scran)
library(igraph)

setwd(Sys.getenv("SOCISL_OUTPUT_ROOT"))

celltypes <- c("Ex_L2_3","Ex_L4","Ex_L4_5","Ex_L5_6","Ex_L5","Oli","Ast","Endo","Mic","In_VIP","In_SST","In_PV (Basket)","In_PV (Chandelier)","In_Rosehip","OPC")
sce <- readH5AD(file.path(Sys.getenv("SOCISL_OUTPUT_ROOT"), "totalAdataAnno040825.h5ad",raw = TRUE)  
ct <- "{variable}"


  #loading in SCE object

  #subsetting for cell type
  scek <- sce[, colData(sce)$cell_type == ct]

  #barcode column! 
  colData(scek)$barcode <- colnames(scek) 

  #subsetting patients based on cognitive status
  library(dplyr)

  data <- read.csv("dataset_652_basic_03-23-2022.csv")

  desiredPatientsM <- c(10205244, 15218541, 15176592, 11165535, 10246987, 11353057, 10277308,10394182, 50405330, 10202345, 11200645, 10490993, 50101659, 70883578, 50405042, 10218339, 10101327, 10101589, 10248033, 10684501, 11157783, 11326252, 11327005,
  11606935, 18414513, 21362537, 22789958, 43485807, 44019405, 66754397,90267190, 50402693)

  desiredPatientsF <- c(21000180, 20252720, 20933324, 21176459, 20173942, 21293107, 98953007, 20358955, 20500815, 21000630, 81852640, 20242958, 20152393, 54122640, 14641578, 81874628, 20898476, 89614402, 60725338, 32383679, 20983400, 21157370, 20153984, 20236838,21411459, 63874408, 50106578, 50108886, 20221273, 20254588, 21140119, 61827429,75675336, 69432088)

  colData(scek)$batch <- as.character(colData(scek)$patient_id)
  colData(scek)$projid <- as.character(colData(scek)$patient_id)
  colData(scek)$projid <- as.factor(colData(scek)$projid)
  data$projid <- as.character(data$projid)

  desired_patientsFR <- data %>%
    filter(projid %in% desiredPatientsF & cogdx == 1)

  desired_patientsFN <- data %>%
    filter(projid %in% desiredPatientsF & cogdx == 4)

  desired_patientsMR <- data %>%
    filter(projid %in% desiredPatientsM & cogdx == 1)

  desired_patientsMN <- data %>%
    filter(projid %in% desiredPatientsM & cogdx == 4)

  colData(scek)$batch <- as.character(colData(scek)$batch)

  desired_patientsFR[[1]] <- as.character(desired_patientsFR[[1]])
  desired_patientsFN[[1]] <- as.character(desired_patientsFN[[1]])
  desired_patientsMR[[1]] <- as.character(desired_patientsMR[[1]])
  desired_patientsMN[[1]] <- as.character(desired_patientsMN[[1]])

  #building female single cell experiment object
  sceFR <- scek[, colData(scek)$batch %in% desired_patientsFR[[1]]]
  sceFN <- scek[, colData(scek)$batch %in% desired_patientsFN[[1]]]


  colData(sceFR)$group <- rep(1, ncol(sceFR))
  colData(sceFN)$group <- rep(0, ncol(sceFN))
  sceF <- cbind(sceFR, sceFN)
  colData(sceF)$projid <- as.character(colData(sceF)$batch)

  #building male single cell experiment object

  sceMR <- scek[, colData(scek)$batch %in% desired_patientsMR[[1]]]
  sceMN <- scek[, colData(scek)$batch %in% desired_patientsMN[[1]]]
  colData(sceMR)$group <- rep(1, ncol(sceMR))
  colData(sceMN)$group <- rep(0, ncol(sceMN))
  data$projid <- as.character(data$projid)

  sceM <- cbind(sceMR, sceMN)
  colData(sceM)$projid <- as.character(colData(sceM)$batch)

  original_colnames <- colnames(sceF)

  meta_df <- as.data.frame(colData(sceF))
  meta_df <- dplyr::left_join(meta_df, data, by = "projid")

  colData(sceF) <- DataFrame(meta_df)
  colnames(sceF) <- original_colnames

  original_colnames <- colnames(sceM)

  meta_df <- as.data.frame(colData(sceM))
  meta_df <- dplyr::left_join(meta_df, data, by = "projid")

  colData(sceM) <- DataFrame(meta_df)
  colnames(sceM) <- original_colnames

  colData(sceF)$age_death <- as.numeric(colData(sceF)$age_death)
  colData(sceF)$pmi <- as.numeric(colData(sceF)$pmi)
  colData(sceF)$educ <- as.numeric(colData(sceF)$educ)
  colData(sceF)$group <- as.numeric(colData(sceF)$group)
  colData(sceF)$msex <- factor(colData(sceF)$msex)
  colData(sceF)$race <- factor(colData(sceF)$race)
  colData(sceF)$apoe <- factor(colData(sceF)$apoe_genotype)
  colData(sceM)$age_death <- as.numeric(colData(sceM)$age_death)
  colData(sceM)$pmi <- as.numeric(colData(sceM)$pmi)
  colData(sceM)$educ <- as.numeric(colData(sceM)$educ)
  colData(sceM)$group <- as.numeric(colData(sceM)$group)
  colData(sceM)$msex <- factor(colData(sceM)$msex)
  colData(sceM)$race <- factor(colData(sceM)$race)
  colData(sceM)$apoe <- factor(colData(sceM)$apoe_genotype)

  # making sure batch and group are factors!
  colData(sceF)$batch <- factor(colData(sceF)$batch)
  colData(sceM)$batch <- factor(colData(sceM)$batch)
  colData(sceF)$group <- factor(colData(sceF)$group, levels = c(0, 1))
  colData(sceM)$group <- factor(colData(sceM)$group, levels = c(0, 1))

  library(nebula)

  #labelling the genes in the raw experiment
  raw_exp <- altExp(sceF, "raw")
  raw_exp <- as(raw_exp, "SingleCellExperiment")
  var_names <- read.csv("var_names.csv", stringsAsFactors = FALSE)
  rownames(raw_exp) <- var_names$gene

  #subsetting for the filtered genes

  common_genes <- intersect(rownames(raw_exp), rownames(sceF))

  assayNames(raw_exp)[assayNames(raw_exp) == "X"] <- "counts"

  #assigning batch and group info
  colData(raw_exp)$batch <- colData(sceF)$batch[match(colnames(raw_exp), rownames(colData(sceF)))]
  colData(raw_exp)$group <- colData(sceF)$group[match(colnames(raw_exp), rownames(colData(sceF)))]
  colData(raw_exp)$pmi <- colData(sceF)$pmi[match(colnames(raw_exp), rownames(colData(sceF)))]
  colData(raw_exp)$educ <- colData(sceF)$educ[match(colnames(raw_exp), rownames(colData(sceF)))]
  colData(raw_exp)$apoe <- colData(sceF)$apoe[match(colnames(raw_exp), rownames(colData(sceF)))]
  colData(raw_exp)$race <- colData(sceF)$race[match(colnames(raw_exp), rownames(colData(sceF)))]
  colData(raw_exp)$age_death <- colData(sceF)$age_death[match(colnames(raw_exp), rownames(colData(sceF)))]
  colData(raw_exp)$age_death <- colData(raw_exp)$age_death / 10
  colData(raw_exp)$pmi <- colData(raw_exp)$pmi / 10

  apoe_odds_map <- c(
    "22" = 0.56,
    "23" = 0.56,
    "24" = 2.64,
    "33" = 1,
    "34" = 3.63,
    "44" = 14.49
  )

  colData(raw_exp)$apoe <- apoe_odds_map[as.character(colData(raw_exp)$apoe)]

  group_metadata <- data.frame(batch = colData(sceF)$batch, group = colData(sceF)$group)

  
  sceE <- raw_exp

  #RUV-III-NB code
  library(variancePartition)
  library(BiocParallel)

  library(NewWave)

  counts_matrix <- counts(sceE)

  keep_feature <- rowSums(counts(sceE) > 0) > ncol(sceE) * 0.1
  sceE <- sceE[keep_feature, ]
  raw_exp <- raw_exp[keep_feature, ]

  colData(raw_exp)$idName <- colData(raw_exp)$batch

  fluidigm_zinb <- newWave(Y=sceE, K = 2, X = "~ group", verbose=TRUE)
  W <- reducedDim(fluidigm_zinb)


  #accounting for ruv factors
  colData(raw_exp) <- cbind(colData(raw_exp), W)


  ruv_formula <- paste("~ group + ", paste(colnames(W), collapse=" + "))
  # model matrix!

  colData(raw_exp)$group <- factor(colData(raw_exp)$group, levels = c("0", "1"))
  design <- model.matrix(as.formula(ruv_formula), data = colData(raw_exp))

  #accounting for ruv factors

  count_df <- as.matrix(assay(raw_exp, "counts"))

  # DGElist + normalization
  dge <- DGEList(counts = count_df)
  dge <- calcNormFactors(dge)  # Normalize counts

  counts_matrix <- counts(raw_exp)
  keep_feature <- rowSums(counts_matrix > 0) > 0
  raw_exp <- raw_exp[keep_feature, ]
  colData(raw_exp)$idName <- colData(raw_exp)$batch
  print(colData(raw_exp))
  seuratdata <- scToNeb(obj = raw_exp, assay = "counts", id = "idName", pred = c("group","W1","W2"))

  mat <- seuratdata$count
  nonzero_cols <- Matrix::colSums(mat) > 0
  mat <- mat[, nonzero_cols]
  seuratdata$count <- seuratdata$count[, nonzero_cols]
  seuratdata$pred <- seuratdata$pred[nonzero_cols, ]
  seuratdata$id   <- seuratdata$id[nonzero_cols]

  pred_df <- as.data.frame(seuratdata$pred)
  pred_df[] <- lapply(pred_df, function(x) {
    if (is.character(x)) {
      if (all(grepl("^-?[0-9.]+$", x))) {
        as.numeric(x)
      } else {
        as.factor(x)
      }
    } else {
      x
    }
  })

 
  pred_df$group <- as.character(pred_df$group)
  pred_df$group <- factor(pred_df$group)
  pred_df$group <- relevel(pred_df$group, ref = "0")

  seuratdata$pred <- pred_df

  dge <- DGEList(counts = seuratdata$count)
  dge <- calcNormFactors(dge)
  offset <- dge$samples$lib.size * dge$samples$norm.factors
  stopifnot(length(offset) == ncol(seuratdata$count))
  offset <- pmax(offset, 1e-6)
  seuratdata$offset <- offset

  design <- model.matrix(~ group + W1 + W2, data = as.data.frame(seuratdata$pred))

  seuratdata=group_cell(count=seuratdata$count,id=seuratdata$id,pred=design)

  pred_df <- as.data.frame(seuratdata$pred)
 
  pred_df$group <- as.character(pred_df$group)
  pred_df$group <- factor(pred_df$group)
  pred_df$group <- relevel(pred_df$group, ref = "0")

  design2 <- model.matrix(~ group1 + W1 + W2, data = as.data.frame(pred_df))

  re2 = nebula(
    count = seuratdata$count, 
    id = seuratdata$id, 
    pred = design2, 
    offset = seuratdata$offset,
    output_re=TRUE,
    method ="HL"
  )


  # saving results
  re2_name <- paste0("Fdejnebula_analysis",ct,"RVNew")
  save(re2, file = paste0("./", re2_name, ".rda"))

  markers <- c("HES4","PDE10A","RPH3A","ST6GAL2","UST")
  print("Nebula Tsai")
  print("fem")
  print(re2$summary[markers,])

  """

base_script_templateBAMDejagerM ="""#!/bin/bash

#SBATCH -n 45                    # Number of cores requested
#SBATCH -t 12:00:00                # Runtime in hours
#SBATCH --mem=800G              # GB memory needed (memory PER CORE)
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-user=nkhera@college.harvard.edu
#SBATCH --mail-type=ALL

source "${CONDA_INIT_SCRIPT}"

echo "Hello, 1!"

cd "${SOCISL_OUTPUT_ROOT}/DeJager"

echo "Hello, 2!"

activate_env "${NEBULA_ENV}"

echo "Hello, 3!"

Rscript male_nebula_{variable}D.Rscript

sbatch male_nebulaEnrich{variable}.sh

"""

base_script_templateNebula ="""library(edgeR)
library(zellkonverter)
library(SingleCellExperiment)
library(magrittr)
library(data.table) 
library(ruvIIInb)
library(scran)
library(igraph)

#loading in SCE object
sce <- readH5AD("totalAdataAnno012125.h5ad",raw = TRUE)  

#subsetting for cell type
sce <- sce[, colData(sce)$cell_type == "{variable}"]

#barcode column! 
colData(sce)$barcode <- colnames(sce) 

#subsetting patients based on cognitive status
library(dplyr)

data <- read.csv("dataset_652_basic_03-23-2022.csv")

desiredPatientsM <- c(50105301, 10518782, 74284255, 10202345, 15113169, 50101659, 11157783, 10253148, 3713990, 10490993, 
        50106730, 50104134, 11444465, 50405330, 10394182, 50402693, 50405042, 18414513, 44299049, 10101589, 
        10277308, 10502798, 11327005)

desiredPatientsF <- c(21151608, 20339740, 21180847, 21408652, 20282974, 20248206, 20153010, 20109020, 21406768, 20208992, 
       20634274, 20297403, 50301125, 50401002, 20254902, 50108886, 20970441, 7265221, 66924745, 
       20344143, 60725338, 50402729, 31908032, 69982533, 50301675, 20380417, 20195344, 20929774, 50104008, 
       20504017, 92393245, 32383679, 20254588, 31726180)

#482428, 18920002

desired_patientsFR <- data %>%
  filter(projid %in% desiredPatientsF & cogdx == 1)

desired_patientsFN <- data %>%
  filter(projid %in% desiredPatientsF & cogdx == 4)

desired_patientsMR <- data %>%
  filter(projid %in% desiredPatientsM & cogdx == 1)

desired_patientsMN <- data %>%
  filter(projid %in% desiredPatientsM & cogdx == 4)

colData(sce)$batch <- as.character(colData(sce)$batch)

desired_patientsFR[[1]] <- as.character(desired_patientsFR[[1]])
desired_patientsFN[[1]] <- as.character(desired_patientsFN[[1]])
desired_patientsMR[[1]] <- as.character(desired_patientsMR[[1]])
desired_patientsMN[[1]] <- as.character(desired_patientsMN[[1]])

#building female single cell experiment object
sceFR <- sce[, colData(sce)$batch %in% desired_patientsFR[[1]]]
sceFN <- sce[, colData(sce)$batch %in% desired_patientsFN[[1]]]

colData(sceFR)$group <- rep(1, ncol(sceFR))
colData(sceFN)$group <- rep(0, ncol(sceFN))
sceF <- cbind(sceFR, sceFN)
colData(sceF)$projid <- as.character(colData(sceF)$batch)

#building male single cell experiment object

sceMR <- sce[, colData(sce)$batch %in% desired_patientsMR[[1]]]
sceMN <- sce[, colData(sce)$batch %in% desired_patientsMN[[1]]]
colData(sceMR)$group <- rep(1, ncol(sceMR))
colData(sceMN)$group <- rep(0, ncol(sceMN))

sceM <- cbind(sceMR, sceMN)
colData(sceM)$projid <- as.character(colData(sceM)$batch)
original_colnames <- colnames(sceF)

meta_df <- as.data.frame(colData(sceF))
meta_df <- merge(meta_df, data, by = "projid", all.x = TRUE)

meta_df <- meta_df[match(colnames(sceF), meta_df$barcode), ]

colData(sceF) <- DataFrame(meta_df)
colnames(sceF) <- original_colnames

original_colnames <- colnames(sceM)

meta_df <- as.data.frame(colData(sceM))
meta_df <- merge(meta_df, data, by = "projid", all.x = TRUE)

meta_df <- meta_df[match(colnames(sceM), meta_df$barcode), ]

colData(sceM) <- DataFrame(meta_df)
colnames(sceM) <- original_colnames


colData(sceF)$age_death <- as.numeric(colData(sceF)$age_death)
colData(sceF)$pmi <- as.numeric(colData(sceF)$pmi)
colData(sceF)$educ <- as.numeric(colData(sceF)$educ)
colData(sceF)$msex <- factor(colData(sceF)$msex)
colData(sceF)$race <- factor(colData(sceF)$race)
colData(sceF)$apoe <- factor(colData(sceF)$apoe_genotype)
colData(sceM)$age_death <- as.numeric(colData(sceM)$age_death)
colData(sceM)$pmi <- as.numeric(colData(sceM)$pmi)
colData(sceM)$educ <- as.numeric(colData(sceM)$educ)
colData(sceM)$msex <- factor(colData(sceM)$msex)
colData(sceM)$race <- factor(colData(sceM)$race)
colData(sceM)$apoe <- factor(colData(sceM)$apoe_genotype)

# making sure batch and group are factors!
colData(sceF)$batch <- factor(colData(sceF)$batch)
colData(sceM)$batch <- factor(colData(sceM)$batch)
colData(sceF)$group <- factor(colData(sceF)$group, levels = c("0", "1"))
colData(sceM)$group <- factor(colData(sceM)$group, levels = c("0", "1"))
library(nebula)

#labelling the genes in the raw experiment
raw_exp <- altExp(sceF, "raw")
raw_exp <- as(raw_exp, "SingleCellExperiment")
var_names <- read.csv("var_names.csv", stringsAsFactors = FALSE)
rownames(raw_exp) <- var_names$gene

#subsetting for the filtered genes

common_genes <- intersect(rownames(raw_exp), rownames(sceF))
raw_exp <- raw_exp[common_genes, ]

assayNames(raw_exp)[assayNames(raw_exp) == "X"] <- "counts"


#assigning batch and group info
colData(raw_exp)$batch <- colData(sceF)$batch[match(colnames(raw_exp), rownames(colData(sceF)))]
colData(raw_exp)$group <- colData(sceF)$group[match(colnames(raw_exp), rownames(colData(sceF)))]
colData(raw_exp)$pmi <- colData(sceF)$pmi[match(colnames(raw_exp), rownames(colData(sceF)))]
colData(raw_exp)$educ <- colData(sceF)$educ[match(colnames(raw_exp), rownames(colData(sceF)))]
colData(raw_exp)$apoe <- colData(sceF)$apoe[match(colnames(raw_exp), rownames(colData(sceF)))]
colData(raw_exp)$race <- colData(sceF)$race[match(colnames(raw_exp), rownames(colData(sceF)))]
colData(raw_exp)$age_death <- colData(sceF)$age_death[match(colnames(raw_exp), rownames(colData(sceF)))]

apoe_odds_map <- c(
  "22" = 0.56,
  "23" = 0.56,
  "24" = 2.64,
  "33" = 1,
  "34" = 3.63,
  "44" = 14.49
)

colData(raw_exp)$apoe <- apoe_odds_map[as.character(colData(raw_exp)$apoe)]

group_metadata <- data.frame(batch = colData(sceF)$batch, group = colData(sceF)$group)


sceE <- raw_exp

#RUV-III-NB code
library(variancePartition)
library(BiocParallel)

library(NewWave)

counts_matrix <- counts(sceE)
keep_feature <- rowSums(counts_matrix > 0) > 0
sceE <- sceE[keep_feature, ]
raw_exp <- raw_exp[keep_feature, ]

colData(raw_exp)$idName <- colData(raw_exp)$batch

fluidigm_zinb <- newWave(Y=sceE, K = 6, X = "~ group", verbose=TRUE)
W <- reducedDim(fluidigm_zinb)


#accounting for ruv factors
colData(raw_exp) <- cbind(colData(raw_exp), W)


ruv_formula <- paste("~ group + apoe + ", paste(colnames(W), collapse=" + "))
# model matrix!

colData(raw_exp)$group <- factor(colData(raw_exp)$group, levels = c("0", "1"))
design <- model.matrix(as.formula(ruv_formula), data = colData(raw_exp))

#accounting for ruv factors

count_df <- as.matrix(assay(raw_exp, "counts"))

# DGElist + normalization
dge <- DGEList(counts = count_df)
dge <- calcNormFactors(dge)  # Normalize counts
offset <- dge$samples$lib.size * dge$samples$norm.factors

dge_log <- edgeR::cpm(dge, log = TRUE, prior.count = 3)  # Log-transformed CPM

expr_mat <- dge_log

meta <- as.data.frame(colData(raw_exp))

meta$group <- factor(meta$group)
meta$batch <- factor(meta$batch)
form <- ~ batch + W1 + W2 + W3 + W4 + W5 + W6 

varPart <- fitExtractVarPartModel(expr_mat, form, meta)
attr(varPart, "errors")
summary_stats <- sort(colMeans(varPart))
print(round(summary_stats, 3))

# seurat input for nebula - creation
seuratdata <- scToNeb(obj = raw_exp, assay = "logcounts", id = "idName", pred = c("group","apoe","W1","W2","W3","W4","W5","W6"))

#creating offset based on past method
seuratdata$offset <- offset

seuratdata$offset <- pmax(seuratdata$offset, 1e-6)  # avoiding log(0) issues


#running nebula
re2 = nebula(
  count = seuratdata$count, 
  id = seuratdata$id, 
  pred = design, 
  offset = seuratdata$offset
)

# saving results
re2_name <- paste0("nebula_analysis{variable}")
save(re2, file = paste0("./", re2_name, ".rda"))
"""


base_script_templateNebulaM ="""library(edgeR)
library(zellkonverter)
library(SingleCellExperiment)
library(magrittr)
library(data.table) 
library(ruvIIInb)
library(scran)
library(igraph)

#loading in SCE object
sce <- readH5AD("totalAdataAnno012125.h5ad",raw = TRUE)  

#subsetting for cell type
sce <- sce[, colData(sce)$cell_type == "{variable}"]

#barcode column! 
colData(sce)$barcode <- colnames(sce) 

#subsetting patients based on cognitive status
library(dplyr)

data <- read.csv("dataset_652_basic_03-23-2022.csv")

desiredPatientsM <- c(50105301, 10518782, 74284255, 10202345, 15113169, 50101659, 11157783, 10253148, 3713990, 10490993, 
        50106730, 50104134, 11444465, 50405330, 10394182, 50402693, 50405042, 18414513, 44299049, 10101589, 
        10277308, 10502798, 11327005)

desiredPatientsF <- c(21151608, 20339740, 21180847, 21408652, 20282974, 20248206, 20153010, 20109020, 21406768, 20208992, 
       20634274, 20297403, 50301125, 50401002, 20254902, 50108886, 20970441, 7265221, 66924745, 
       20344143, 60725338, 50402729, 31908032, 69982533, 50301675, 20380417, 20195344, 20929774, 50104008, 
       20504017, 92393245, 32383679, 20254588, 31726180)

#482428, 18920002

desired_patientsFR <- data %>%
  filter(projid %in% desiredPatientsF & cogdx == 1)

desired_patientsFN <- data %>%
  filter(projid %in% desiredPatientsF & cogdx == 4)

desired_patientsMR <- data %>%
  filter(projid %in% desiredPatientsM & cogdx == 1)

desired_patientsMN <- data %>%
  filter(projid %in% desiredPatientsM & cogdx == 4)

colData(sce)$batch <- as.character(colData(sce)$batch)

desired_patientsFR[[1]] <- as.character(desired_patientsFR[[1]])
desired_patientsFN[[1]] <- as.character(desired_patientsFN[[1]])
desired_patientsMR[[1]] <- as.character(desired_patientsMR[[1]])
desired_patientsMN[[1]] <- as.character(desired_patientsMN[[1]])

#building female single cell experiment object
sceFR <- sce[, colData(sce)$batch %in% desired_patientsFR[[1]]]
sceFN <- sce[, colData(sce)$batch %in% desired_patientsFN[[1]]]

colData(sceFR)$group <- rep(1, ncol(sceFR))
colData(sceFN)$group <- rep(0, ncol(sceFN))
sceF <- cbind(sceFR, sceFN)
colData(sceF)$projid <- as.character(colData(sceF)$batch)

#building male single cell experiment object

sceMR <- sce[, colData(sce)$batch %in% desired_patientsMR[[1]]]
sceMN <- sce[, colData(sce)$batch %in% desired_patientsMN[[1]]]
colData(sceMR)$group <- rep(1, ncol(sceMR))
colData(sceMN)$group <- rep(0, ncol(sceMN))

sceM <- cbind(sceMR, sceMN)
colData(sceM)$projid <- as.character(colData(sceM)$batch)
original_colnames <- colnames(sceF)

meta_df <- as.data.frame(colData(sceF))
meta_df <- merge(meta_df, data, by = "projid", all.x = TRUE)

meta_df <- meta_df[match(colnames(sceF), meta_df$barcode), ]

colData(sceF) <- DataFrame(meta_df)
colnames(sceF) <- original_colnames

original_colnames <- colnames(sceM)

meta_df <- as.data.frame(colData(sceM))
meta_df <- merge(meta_df, data, by = "projid", all.x = TRUE)

meta_df <- meta_df[match(colnames(sceM), meta_df$barcode), ]

colData(sceM) <- DataFrame(meta_df)
colnames(sceM) <- original_colnames


colData(sceF)$age_death <- as.numeric(colData(sceF)$age_death)
colData(sceF)$pmi <- as.numeric(colData(sceF)$pmi)
colData(sceF)$educ <- as.numeric(colData(sceF)$educ)
colData(sceF)$msex <- factor(colData(sceF)$msex)
colData(sceF)$race <- factor(colData(sceF)$race)
colData(sceF)$apoe <- factor(colData(sceF)$apoe_genotype)
colData(sceM)$age_death <- as.numeric(colData(sceM)$age_death)
colData(sceM)$pmi <- as.numeric(colData(sceM)$pmi)
colData(sceM)$educ <- as.numeric(colData(sceM)$educ)
colData(sceM)$msex <- factor(colData(sceM)$msex)
colData(sceM)$race <- factor(colData(sceM)$race)
colData(sceM)$apoe <- factor(colData(sceM)$apoe_genotype)

# making sure batch and group are factors!
colData(sceF)$batch <- factor(colData(sceF)$batch)
colData(sceM)$batch <- factor(colData(sceM)$batch)
colData(sceF)$group <- factor(colData(sceF)$group, levels = c("0", "1"))
colData(sceM)$group <- factor(colData(sceM)$group, levels = c("0", "1"))
library(nebula)

#labelling the genes in the raw experiment
raw_exp <- altExp(sceM, "raw")
raw_exp <- as(raw_exp, "SingleCellExperiment")
var_names <- read.csv("var_names.csv", stringsAsFactors = FALSE)
rownames(raw_exp) <- var_names$gene

#subsetting for the filtered genes

common_genes <- intersect(rownames(raw_exp), rownames(sceF))
raw_exp <- raw_exp[common_genes, ]

assayNames(raw_exp)[assayNames(raw_exp) == "X"] <- "counts"


#assigning batch and group info
colData(raw_exp)$batch <- colData(sceM)$batch[match(colnames(raw_exp), rownames(colData(sceM)))]
colData(raw_exp)$group <- colData(sceM)$group[match(colnames(raw_exp), rownames(colData(sceM)))]
colData(raw_exp)$pmi <- colData(sceM)$pmi[match(colnames(raw_exp), rownames(colData(sceM)))]
colData(raw_exp)$educ <- colData(sceM)$educ[match(colnames(raw_exp), rownames(colData(sceM)))]
colData(raw_exp)$apoe <- colData(sceM)$apoe[match(colnames(raw_exp), rownames(colData(sceM)))]
colData(raw_exp)$race <- colData(sceM)$race[match(colnames(raw_exp), rownames(colData(sceM)))]
colData(raw_exp)$age_death <- colData(sceM)$age_death[match(colnames(raw_exp), rownames(colData(sceM)))]

apoe_odds_map <- c(
  "22" = 0.56,
  "23" = 0.56,
  "24" = 2.64,
  "33" = 1,
  "34" = 3.63,
  "44" = 14.49
)

colData(raw_exp)$apoe <- apoe_odds_map[as.character(colData(raw_exp)$apoe)]

group_metadata <- data.frame(batch = colData(sceM)$batch, group = colData(sceM)$group)


sceE <- raw_exp

#RUV-III-NB code
library(variancePartition)
library(BiocParallel)

library(NewWave)

counts_matrix <- counts(sceE)
keep_feature <- rowSums(counts_matrix > 0) > 0
sceE <- sceE[keep_feature, ]
raw_exp <- raw_exp[keep_feature, ]

colData(raw_exp)$idName <- colData(raw_exp)$batch

fluidigm_zinb <- newWave(Y=sceE, K = 6, X = "~ group", verbose=TRUE)
W <- reducedDim(fluidigm_zinb)


#accounting for ruv factors
colData(raw_exp) <- cbind(colData(raw_exp), W)


ruv_formula <- paste("~ group + apoe + ", paste(colnames(W), collapse=" + "))
# model matrix!

colData(raw_exp)$group <- factor(colData(raw_exp)$group, levels = c("0", "1"))
design <- model.matrix(as.formula(ruv_formula), data = colData(raw_exp))

#accounting for ruv factors

count_df <- as.matrix(assay(raw_exp, "counts"))

# DGElist + normalization
dge <- DGEList(counts = count_df)
dge <- calcNormFactors(dge)  # Normalize counts
offset <- dge$samples$lib.size * dge$samples$norm.factors

dge_log <- edgeR::cpm(dge, log = TRUE, prior.count = 3)  # Log-transformed CPM

expr_mat <- dge_log

meta <- as.data.frame(colData(raw_exp))

meta$group <- factor(meta$group)
meta$batch <- factor(meta$batch)
form <- ~ batch + W1 + W2 + W3 + W4 + W5 + W6 

varPart <- fitExtractVarPartModel(expr_mat, form, meta)
attr(varPart, "errors")
summary_stats <- sort(colMeans(varPart))
print(round(summary_stats, 3))

# seurat input for nebula - creation
seuratdata <- scToNeb(obj = raw_exp, assay = "logcounts", id = "idName", pred = c("group","apoe","W1","W2","W3","W4","W5","W6"))

#creating offset based on past method
seuratdata$offset <- offset

seuratdata$offset <- pmax(seuratdata$offset, 1e-6)  # avoiding log(0) issues


#running nebula
re2 = nebula(
  count = seuratdata$count, 
  id = seuratdata$id, 
  pred = design, 
  offset = seuratdata$offset
)

# saving results
re2_name <- paste0("male_nebula_analysis{variable}")
save(re2, file = paste0("./", re2_name, ".rda"))"""


base_script_templateNebulaDejager ="""
library(edgeR)
library(zellkonverter)
library(SingleCellExperiment)
library(magrittr)
library(data.table) 
library(ruvIIInb)
library(scran)
library(igraph)

#loading in SCE object
sce <- readH5AD(file.path(Sys.getenv("SOCISL_OUTPUT_ROOT"), "totalAdataAnno040825.h5ad",raw = TRUE)  

#subsetting for cell type
sce <- sce[, colData(sce)$cell_type == "{variable}"]

#barcode column! 
colData(sce)$barcode <- colnames(sce) 

#subsetting patients based on cognitive status
library(dplyr)

data <- read.csv("dataset_652_basic_03-23-2022.csv")

desiredPatientsM <- c(10205244, 15218541, 15176592, 11165535, 10246987, 11353057, 10277308, 10394182, 50405330, 10202345, 11200645, 10490993, 50101659, 70883578, 50405042, 10218339,46547648, 43485807, 44019405, 50402693, 50108912, 82317494, 11157783, 10248033,11327005, 46000440, 50106730, 21362537, 90267190, 10684501, 66754397, 10101327)

desiredPatientsF <- c(21000180, 20252720, 20933324, 21176459, 20173942, 21293107, 98953007, 20358955, 20500815, 21000630, 81852640, 20242958, 20152393, 54122640, 14641578, 81874628, 20898476, 89614402, 60725338, 32383679, 20983400, 21157370, 20153984, 20236838,21411459, 63874408, 50106578, 50108886, 20221273, 20254588, 21140119, 61827429,75675336, 69432088)

desired_patientsFR <- data %>%
  filter(projid %in% desiredPatientsF & cogdx == 1)

desired_patientsFN <- data %>%
  filter(projid %in% desiredPatientsF & cogdx == 4)

desired_patientsMR <- data %>%
  filter(projid %in% desiredPatientsM & cogdx == 1)

desired_patientsMN <- data %>%
  filter(projid %in% desiredPatientsM & cogdx == 4)

colData(sce)$batch <- as.character(colData(sce)$patient_id)
colData(sce)$batch[colData(sce)$batch == "2"] <- "60725338"

desired_patientsFR[[1]] <- as.character(desired_patientsFR[[1]])
desired_patientsFN[[1]] <- as.character(desired_patientsFN[[1]])
desired_patientsMR[[1]] <- as.character(desired_patientsMR[[1]])
desired_patientsMN[[1]] <- as.character(desired_patientsMN[[1]])

#building female single cell experiment object
sceFR <- sce[, colData(sce)$batch %in% desired_patientsFR[[1]]]
sceFN <- sce[, colData(sce)$batch %in% desired_patientsFN[[1]]]

colData(sceFR)$group <- rep(1, ncol(sceFR))
colData(sceFN)$group <- rep(0, ncol(sceFN))
sceF <- cbind(sceFR, sceFN)
colData(sceF)$projid <- as.character(colData(sceF)$patient_id)

#building male single cell experiment object

sceMR <- sce[, colData(sce)$batch %in% desired_patientsMR[[1]]]
sceMN <- sce[, colData(sce)$batch %in% desired_patientsMN[[1]]]
colData(sceMR)$group <- rep(1, ncol(sceMR))
colData(sceMN)$group <- rep(0, ncol(sceMN))

sceM <- cbind(sceMR, sceMN)
colData(sceM)$projid <- as.character(colData(sceM)$batch)
original_colnames <- colnames(sceF)

meta_df <- as.data.frame(colData(sceF))
meta_df <- merge(meta_df, data, by = "projid", all.x = TRUE)

meta_df <- meta_df[match(colnames(sceF), meta_df$barcode), ]

colData(sceF) <- DataFrame(meta_df)
colnames(sceF) <- original_colnames

original_colnames <- colnames(sceM)

meta_df <- as.data.frame(colData(sceM))
meta_df <- merge(meta_df, data, by = "projid", all.x = TRUE)
meta_df <- meta_df[match(colnames(sceM), meta_df$barcode), ]

colData(sceM) <- DataFrame(meta_df)
colnames(sceM) <- original_colnames


colData(sceF)$age_death <- as.numeric(colData(sceF)$age_death)
colData(sceF)$pmi <- as.numeric(colData(sceF)$pmi)
colData(sceF)$educ <- as.numeric(colData(sceF)$educ)
colData(sceF)$msex <- factor(colData(sceF)$msex)
colData(sceF)$race <- factor(colData(sceF)$race)
colData(sceF)$apoe <- factor(colData(sceF)$apoe_genotype)
colData(sceM)$age_death <- as.numeric(colData(sceM)$age_death)
colData(sceM)$pmi <- as.numeric(colData(sceM)$pmi)
colData(sceM)$educ <- as.numeric(colData(sceM)$educ)
colData(sceM)$msex <- factor(colData(sceM)$msex)
colData(sceM)$race <- factor(colData(sceM)$race)
colData(sceM)$apoe <- factor(colData(sceM)$apoe_genotype)


# making sure batch and group are factors!
colData(sceF)$batch <- factor(colData(sceF)$batch)
colData(sceM)$batch <- factor(colData(sceM)$batch)
colData(sceF)$group <- factor(colData(sceF)$group, levels = c("0", "1"))
colData(sceM)$group <- factor(colData(sceM)$group, levels = c("0", "1"))


library(nebula)

#labelling the genes in the raw experiment
raw_exp <- altExp(sceF, "raw")
raw_exp <- as(raw_exp, "SingleCellExperiment")
var_names <- read.csv("var_names.csv", stringsAsFactors = FALSE)
rownames(raw_exp) <- var_names$gene

#subsetting for the filtered genes

common_genes <- intersect(rownames(raw_exp), rownames(sceF))
raw_exp <- raw_exp[common_genes, ]

assayNames(raw_exp)[assayNames(raw_exp) == "X"] <- "counts"


#assigning batch and group info
colData(raw_exp)$batch <- colData(sceF)$batch[match(colnames(raw_exp), rownames(colData(sceF)))]
colData(raw_exp)$group <- colData(sceF)$group[match(colnames(raw_exp), rownames(colData(sceF)))]
colData(raw_exp)$pmi <- colData(sceF)$pmi[match(colnames(raw_exp), rownames(colData(sceF)))]
colData(raw_exp)$educ <- colData(sceF)$educ[match(colnames(raw_exp), rownames(colData(sceF)))]
colData(raw_exp)$apoe <- colData(sceF)$apoe[match(colnames(raw_exp), rownames(colData(sceF)))]
colData(raw_exp)$race <- colData(sceF)$race[match(colnames(raw_exp), rownames(colData(sceF)))]
colData(raw_exp)$age_death <- colData(sceF)$age_death[match(colnames(raw_exp), rownames(colData(sceF)))]

# valid_cells <- complete.cases(colData(raw_exp)[, c("pmi", "educ")])
# raw_exp <- raw_exp[, valid_cells]

apoe_odds_map <- c(
  "22" = 0.56,
  "23" = 0.56,
  "24" = 2.64,
  "33" = 1,
  "34" = 3.63,
  "44" = 14.49
)

colData(raw_exp)$apoe <- apoe_odds_map[as.character(colData(raw_exp)$apoe)]

group_metadata <- data.frame(batch = colData(sceF)$batch, group = colData(sceF)$group)


# seurat input for nebula - creation
sceE <- raw_exp

library(DelayedArray)

library(variancePartition)
library(BiocParallel)

library(NewWave)

counts_matrix <- counts(sceE)
keep_feature <- rowSums(counts_matrix > 0) > 0
sceE <- sceE[keep_feature, ]
raw_exp <- raw_exp[keep_feature, ]

colData(raw_exp)$idName <- colData(raw_exp)$batch


fluidigm_zinb <- newWave(Y=sceE, K = 6, X = "~ group", verbose=TRUE)
W <- reducedDim(fluidigm_zinb)

#accounting for ruv factors
colData(raw_exp) <- cbind(colData(raw_exp), W)

ruv_formula <- paste("~ group + apoe + ", paste(colnames(W), collapse=" + "))

# model matrix!
colData(raw_exp)$group <- factor(colData(raw_exp)$group, levels = c("0", "1"))
design <- model.matrix(as.formula(ruv_formula), data = colData(raw_exp))

count_df <- as.matrix(assay(raw_exp, "counts"))

# DGElist + normalization
dge <- DGEList(counts = count_df)
dge <- calcNormFactors(dge)  # Normalize counts
offset <- dge$samples$lib.size * dge$samples$norm.factors

dge_log <- edgeR::cpm(dge, log = TRUE, prior.count = 3)  # Log-transformed CPM

expr_mat <- dge_log

meta <- as.data.frame(colData(raw_exp))

meta$group <- factor(meta$group)
meta$batch <- factor(meta$batch)
form <- ~ batch + W1 + W2 + W3 + W4 + W5 + W6 

varPart <- fitExtractVarPartModel(expr_mat, form, meta)
attr(varPart, "errors")
summary_stats <- sort(colMeans(varPart))

# seurat input for nebula - creation
seuratdata <- scToNeb(obj = raw_exp, assay = "logcounts", id = "idName", pred = c("group","apoe","W1","W2","W3","W4","W5","W6"))

#creating offset based on past method
seuratdata$offset <- offset
seuratdata$offset <- pmax(seuratdata$offset, 1e-6)  # avoiding log(0) issues

seuratdata=group_cell(count=seuratdata$count,id=seuratdata$id,pred=design)

#running nebula
re2 = nebula(
  count = seuratdata$count, 
  id = seuratdata$id, 
  pred = design, 
  offset = seuratdata$offset
)

# saving results
re2_name <- paste0("nebula_analysis{variable}")
save(re2, file = paste0("./", re2_name, ".rda"))
"""



base_script_templateNebulaDejagerM ="""
library(edgeR)
library(zellkonverter)
library(SingleCellExperiment)
library(magrittr)
library(data.table) 
library(ruvIIInb)
library(scran)
library(igraph)

#loading in SCE object
sce <- readH5AD(file.path(Sys.getenv("SOCISL_OUTPUT_ROOT"), "totalAdataAnno040825.h5ad",raw = TRUE)  

#subsetting for cell type
sce <- sce[, colData(sce)$cell_type == "{variable}"]

#barcode column! 
colData(sce)$barcode <- colnames(sce) 

#subsetting patients based on cognitive status
library(dplyr)

data <- read.csv("dataset_652_basic_03-23-2022.csv")

desiredPatientsM <- c(10205244, 15218541, 15176592, 11165535, 10246987, 11353057, 10277308, 10394182, 50405330, 10202345, 11200645, 10490993, 50101659, 70883578, 50405042, 10218339,46547648, 43485807, 44019405, 50402693, 50108912, 82317494, 11157783, 10248033,11327005, 46000440, 50106730, 21362537, 90267190, 10684501, 66754397, 10101327)

desiredPatientsF <- c(21000180, 20252720, 20933324, 21176459, 20173942, 21293107, 98953007, 20358955, 20500815, 21000630, 81852640, 20242958, 20152393, 54122640, 14641578, 81874628, 20898476, 89614402, 60725338, 32383679, 20983400, 21157370, 20153984, 20236838,21411459, 63874408, 50106578, 50108886, 20221273, 20254588, 21140119, 61827429,75675336, 69432088)

desired_patientsFR <- data %>%
  filter(projid %in% desiredPatientsF & cogdx == 1)

desired_patientsFN <- data %>%
  filter(projid %in% desiredPatientsF & cogdx == 4)

desired_patientsMR <- data %>%
  filter(projid %in% desiredPatientsM & cogdx == 1)

desired_patientsMN <- data %>%
  filter(projid %in% desiredPatientsM & cogdx == 4)

colData(sce)$batch <- as.character(colData(sce)$patient_id)
colData(sce)$batch[colData(sce)$batch == "2"] <- "60725338"

desired_patientsFR[[1]] <- as.character(desired_patientsFR[[1]])
desired_patientsFN[[1]] <- as.character(desired_patientsFN[[1]])
desired_patientsMR[[1]] <- as.character(desired_patientsMR[[1]])
desired_patientsMN[[1]] <- as.character(desired_patientsMN[[1]])

#building female single cell experiment object
sceFR <- sce[, colData(sce)$batch %in% desired_patientsFR[[1]]]
sceFN <- sce[, colData(sce)$batch %in% desired_patientsFN[[1]]]

colData(sceFR)$group <- rep(1, ncol(sceFR))
colData(sceFN)$group <- rep(0, ncol(sceFN))
sceF <- cbind(sceFR, sceFN)
colData(sceF)$projid <- as.character(colData(sceF)$patient_id)

#building male single cell experiment object

sceMR <- sce[, colData(sce)$batch %in% desired_patientsMR[[1]]]
sceMN <- sce[, colData(sce)$batch %in% desired_patientsMN[[1]]]
colData(sceMR)$group <- rep(1, ncol(sceMR))
colData(sceMN)$group <- rep(0, ncol(sceMN))

sceM <- cbind(sceMR, sceMN)
colData(sceM)$projid <- as.character(colData(sceM)$batch)
original_colnames <- colnames(sceF)

meta_df <- as.data.frame(colData(sceF))
meta_df <- merge(meta_df, data, by = "projid", all.x = TRUE)

meta_df <- meta_df[match(colnames(sceF), meta_df$barcode), ]

colData(sceF) <- DataFrame(meta_df)
colnames(sceF) <- original_colnames

original_colnames <- colnames(sceM)

meta_df <- as.data.frame(colData(sceM))
meta_df <- merge(meta_df, data, by = "projid", all.x = TRUE)
meta_df <- meta_df[match(colnames(sceM), meta_df$barcode), ]

colData(sceM) <- DataFrame(meta_df)
colnames(sceM) <- original_colnames


colData(sceF)$age_death <- as.numeric(colData(sceF)$age_death)
colData(sceF)$pmi <- as.numeric(colData(sceF)$pmi)
colData(sceF)$educ <- as.numeric(colData(sceF)$educ)
colData(sceF)$msex <- factor(colData(sceF)$msex)
colData(sceF)$race <- factor(colData(sceF)$race)
colData(sceF)$apoe <- factor(colData(sceF)$apoe_genotype)
colData(sceM)$age_death <- as.numeric(colData(sceM)$age_death)
colData(sceM)$pmi <- as.numeric(colData(sceM)$pmi)
colData(sceM)$educ <- as.numeric(colData(sceM)$educ)
colData(sceM)$msex <- factor(colData(sceM)$msex)
colData(sceM)$race <- factor(colData(sceM)$race)
colData(sceM)$apoe <- factor(colData(sceM)$apoe_genotype)


# making sure batch and group are factors!
colData(sceF)$batch <- factor(colData(sceF)$batch)
colData(sceM)$batch <- factor(colData(sceM)$batch)
colData(sceF)$group <- factor(colData(sceF)$group, levels = c("0", "1"))
colData(sceM)$group <- factor(colData(sceM)$group, levels = c("0", "1"))


library(nebula)

#labelling the genes in the raw experiment
raw_exp <- altExp(sceF, "raw")
raw_exp <- as(raw_exp, "SingleCellExperiment")
var_names <- read.csv("var_names.csv", stringsAsFactors = FALSE)
rownames(raw_exp) <- var_names$gene

#subsetting for the filtered genes

common_genes <- intersect(rownames(raw_exp), rownames(sceF))
raw_exp <- raw_exp[common_genes, ]

assayNames(raw_exp)[assayNames(raw_exp) == "X"] <- "counts"


#assigning batch and group info
colData(raw_exp)$batch <- colData(sceM)$batch[match(colnames(raw_exp), rownames(colData(sceM)))]
colData(raw_exp)$group <- colData(sceM)$group[match(colnames(raw_exp), rownames(colData(sceM)))]
colData(raw_exp)$pmi <- colData(sceM)$pmi[match(colnames(raw_exp), rownames(colData(sceM)))]
colData(raw_exp)$educ <- colData(sceM)$educ[match(colnames(raw_exp), rownames(colData(sceM)))]
colData(raw_exp)$apoe <- colData(sceM)$apoe[match(colnames(raw_exp), rownames(colData(sceM)))]
colData(raw_exp)$race <- colData(sceM)$race[match(colnames(raw_exp), rownames(colData(sceM)))]
colData(raw_exp)$age_death <- colData(sceM)$age_death[match(colnames(raw_exp), rownames(colData(sceM)))]

# valid_cells <- complete.cases(colData(raw_exp)[, c("pmi", "educ")])
# raw_exp <- raw_exp[, valid_cells]

apoe_odds_map <- c(
  "22" = 0.56,
  "23" = 0.56,
  "24" = 2.64,
  "33" = 1,
  "34" = 3.63,
  "44" = 14.49
)

colData(raw_exp)$apoe <- apoe_odds_map[as.character(colData(raw_exp)$apoe)]

group_metadata <- data.frame(batch = colData(sceM)$batch, group = colData(sceM)$group)


# seurat input for nebula - creation
sceE <- raw_exp

library(DelayedArray)

library(variancePartition)
library(BiocParallel)

library(NewWave)

counts_matrix <- counts(sceE)
keep_feature <- rowSums(counts_matrix > 0) > 0
sceE <- sceE[keep_feature, ]
raw_exp <- raw_exp[keep_feature, ]

colData(raw_exp)$idName <- colData(raw_exp)$batch

fluidigm_zinb <- newWave(Y=sceE, K = 6, X = "~ group", verbose=TRUE)
W <- reducedDim(fluidigm_zinb)

#accounting for ruv factors
colData(raw_exp) <- cbind(colData(raw_exp), W)

ruv_formula <- paste("~ group + apoe + ", paste(colnames(W), collapse=" + "))

# model matrix!
colData(raw_exp)$group <- factor(colData(raw_exp)$group, levels = c("0", "1"))
design <- model.matrix(as.formula(ruv_formula), data = colData(raw_exp))

count_df <- as.matrix(assay(raw_exp, "counts"))

# DGElist + normalization
dge <- DGEList(counts = count_df)
dge <- calcNormFactors(dge)  # Normalize counts
offset <- dge$samples$lib.size * dge$samples$norm.factors

dge_log <- edgeR::cpm(dge, log = TRUE, prior.count = 3)  # Log-transformed CPM

expr_mat <- dge_log

meta <- as.data.frame(colData(raw_exp))

meta$group <- factor(meta$group)
meta$batch <- factor(meta$batch)
form <- ~ batch + W1 + W2 + W3 + W4 + W5 + W6 

varPart <- fitExtractVarPartModel(expr_mat, form, meta)

attr(varPart, "errors")
summary_stats <- sort(colMeans(varPart))

# seurat input for nebula - creation
seuratdata <- scToNeb(obj = raw_exp, assay = "logcounts", id = "idName", pred = c("group","apoe","W1","W2","W3","W4","W5","W6"))

#creating offset based on past method
seuratdata$offset <- offset
seuratdata$offset <- pmax(seuratdata$offset, 1e-6)  # avoiding log(0) issues

seuratdata=group_cell(count=seuratdata$count,id=seuratdata$id,pred=design)

#running nebula
re2 = nebula(
  count = seuratdata$count, 
  id = seuratdata$id, 
  pred = design, 
  offset = seuratdata$offset
)

# saving results
re2_name <- paste0("male_nebula_analysis{variable}")
save(re2, file = paste0("./", re2_name, ".rda"))
"""


base_script_templateNebulaEnrichment ="""

library(ggplot2)
library(dplyr)
library(readr)
library(tibble)
require(WebGestaltR)
load("Fnebula_analysis{variable}RVNew.rda")
diff.analysis <- re2
#VERSION TO GENERATE DATAFRAME. 
master.webgestaltR <- vector(mode = "list")
GSEA.list <- c("geneontology_Biological_Process_noRedundant","geneontology_Cellular_Component_noRedundant","geneontology_Molecular_Function_noRedundant","pathway_KEGG","pathway_Panther","pathway_Reactome","pathway_Wikipathway","network_Transcription_Factor_target")
diff.analysis$summary$rank <- sign(diff.analysis$summary$logFC_group1)*-log10(diff.analysis$summary$p_group1)
gsea.rank <- diff.analysis$summary[,c("gene","rank")]
gsea.rank <- gsea.rank[order(-gsea.rank$rank),]
for (j in 1:length(GSEA.list)){{
    listname <- paste0("{variable}","_",GSEA.list[j])
    master.webgestaltR[[listname]] <- WebGestaltR(enrichMethod = "GSEA", 
        organism = "hsapiens",
        minNum = 10,
        fdrThr = 0.2,
        enrichDatabase= GSEA.list[j],  
        enrichDatabaseType="genesymbol",
        interestGene =  gsea.rank,
        interestGeneType="genesymbol",
        sigMethod = "fdr",
        is.output = FALSE)
        saveRDS(master.webgestaltR,"overallWebGestaltRResult{variable}F.rds")
    }}
saveRDS(master.webgestaltR,"overallWebGestaltRResult{variable}F.rds")

"""


base_script_templateNebulaEnrichmentM ="""

library(ggplot2)
library(dplyr)
library(readr)
library(tibble)
require(WebGestaltR)
load("Mnebula_analysis{variable}RVNew.rda")
diff.analysis <- re2
#VERSION TO GENERATE DATAFRAME. 
master.webgestaltR <- vector(mode = "list")
GSEA.list <- c("geneontology_Biological_Process_noRedundant","geneontology_Cellular_Component_noRedundant","geneontology_Molecular_Function_noRedundant","pathway_KEGG","pathway_Panther","pathway_Reactome","pathway_Wikipathway","network_Transcription_Factor_target")
diff.analysis$summary$rank <- sign(diff.analysis$summary$logFC_group1)*-log10(diff.analysis$summary$p_group1)
gsea.rank <- diff.analysis$summary[,c("gene","rank")]
gsea.rank <- gsea.rank[order(-gsea.rank$rank),]
for (j in 1:length(GSEA.list)){{
    listname <- paste0("{variable}","_",GSEA.list[j])
    master.webgestaltR[[listname]] <- WebGestaltR(enrichMethod = "GSEA", 
        organism = "hsapiens",
        minNum = 10,
        fdrThr = 0.2,
        enrichDatabase= GSEA.list[j],  
        enrichDatabaseType="genesymbol",
        interestGene =  gsea.rank,
        interestGeneType="genesymbol",
        sigMethod = "fdr",
        is.output = FALSE)
        saveRDS(master.webgestaltR,"overallWebGestaltRResult{variable}M.rds")
    }}
saveRDS(master.webgestaltR,"overallWebGestaltRResult{variable}M.rds")


"""



base_script_templateNebulaEnrichmentD ="""

library(ggplot2)
library(dplyr)
library(readr)
library(tibble)
require(WebGestaltR)
load("Fdejnebula_analysis{variable}RVNew.rda")
diff.analysis <- re2
#VERSION TO GENERATE DATAFRAME. 
master.webgestaltR <- vector(mode = "list")
GSEA.list <- c("geneontology_Biological_Process_noRedundant","geneontology_Cellular_Component_noRedundant","geneontology_Molecular_Function_noRedundant","pathway_KEGG","pathway_Panther","pathway_Reactome","pathway_Wikipathway","network_Transcription_Factor_target")
diff.analysis$summary$rank <- sign(diff.analysis$summary$logFC_group1)*-log10(diff.analysis$summary$p_group1)
gsea.rank <- diff.analysis$summary[,c("gene","rank")]
gsea.rank <- gsea.rank[order(-gsea.rank$rank),]
for (j in 1:length(GSEA.list)){{
    listname <- paste0("{variable}","_",GSEA.list[j])
    master.webgestaltR[[listname]] <- WebGestaltR(enrichMethod = "GSEA", 
        organism = "hsapiens",
        minNum = 10,
        fdrThr = 0.2,
        enrichDatabase= GSEA.list[j],  
        enrichDatabaseType="genesymbol",
        interestGene =  gsea.rank,
        interestGeneType="genesymbol",
        sigMethod = "fdr",
        is.output = FALSE)
        saveRDS(master.webgestaltR,"overallWebGestaltRResult{variable}FDej.rds")
    }}
saveRDS(master.webgestaltR,"overallWebGestaltRResult{variable}FDej.rds")

"""


base_script_templateNebulaEnrichmentMD ="""

library(ggplot2)
library(dplyr)
library(readr)
library(tibble)
require(WebGestaltR)
load("dejnebula_analysis{variable}RVNew.rda")
diff.analysis <- re2
#VERSION TO GENERATE DATAFRAME. 
master.webgestaltR <- vector(mode = "list")
GSEA.list <- c("geneontology_Biological_Process_noRedundant","geneontology_Cellular_Component_noRedundant","geneontology_Molecular_Function_noRedundant","pathway_KEGG","pathway_Panther","pathway_Reactome","pathway_Wikipathway","network_Transcription_Factor_target")
diff.analysis$summary$rank <- sign(diff.analysis$summary$logFC_group1)*-log10(diff.analysis$summary$p_group1)
gsea.rank <- diff.analysis$summary[,c("gene","rank")]
gsea.rank <- gsea.rank[order(-gsea.rank$rank),]
for (j in 1:length(GSEA.list)){{
    listname <- paste0("{variable}","_",GSEA.list[j])
    master.webgestaltR[[listname]] <- WebGestaltR(enrichMethod = "GSEA", 
        organism = "hsapiens",
        minNum = 10,
        fdrThr = 0.2,
        enrichDatabase= GSEA.list[j],  
        enrichDatabaseType="genesymbol",
        interestGene =  gsea.rank,
        interestGeneType="genesymbol",
        sigMethod = "fdr",
        is.output = FALSE)
        saveRDS(master.webgestaltR,"overallWebGestaltRResult{variable}MDej.rds")
    }}
saveRDS(master.webgestaltR,"overallWebGestaltRResult{variable}MDej.rds")


"""


base_script_templateNebulaRRHO2 ="""


library(RRHO)
library(RRHO2)
load(file.path(Sys.getenv("SOCISL_OUTPUT_ROOT"), "Fnebula_analysis{variable}RVNew.rda")
df1 <- re2
load(file.path(Sys.getenv("SOCISL_OUTPUT_ROOT"), "Fdejnebula_analysis{variable}RVNew.rda")
df2 <- re2

df1$summary$signed_stat <- df1$summary$logFC_group1 / df1$summary$se_group1
genes1 <- setNames(df1$summary$signed_stat, df1$summary$gene)
genes1 <- sort(genes1, decreasing = TRUE)
df2$summary$signed_stat <- df2$summary$logFC_group1 / df2$summary$se_group1
genes2 <- setNames(df2$summary$signed_stat, df2$summary$gene)
genes2 <- sort(genes2, decreasing = TRUE)

common_genes <- intersect(names(genes1), names(genes2))
genes1_sub <- genes1[common_genes]
genes2_sub <- genes2[common_genes]

str(genes2)
genes1_df <- data.frame(Gene = names(genes1_sub), logFC = genes1_sub)
genes2_df <- data.frame(Gene = names(genes2_sub), logFC = genes2_sub)

RRHO_obj <-  RRHO2_initialize(genes1_df, genes2_df, labels = c("tsai", "dejager"), log10.ind=FALSE)

outp <- cor.test(genes1_sub, genes2_sub, method = "spearman")

print(outp)

png("RRHO_heatmap{variable}.png", width=800, height=800)
RRHO2_heatmap(RRHO_obj)

lab <- sprintf("Spearman ρ = %.2f\np = %.2e",
               outp$estimate, outp$p.value)

mtext(lab,
      side = 3,        # top
      adj = 1,         # right aligned
      line = -1,       # move inside plot
      cex = 0.9)

dev.off()

load(file.path(Sys.getenv("SOCISL_OUTPUT_ROOT"), "Mnebula_analysis{variable}RVNew.rda")
df1M <- re2
load(file.path(Sys.getenv("SOCISL_OUTPUT_ROOT"), "dejnebula_analysis{variable}RVNew.rda")
df2M <- re2

df1M$summary$signed_stat <- df1M$summary$logFC_group1 / df1M$summary$se_group1
genes1M <- setNames(df1M$summary$signed_stat, df1M$summary$gene)
genes1M <- sort(genes1M, decreasing = TRUE)
df2M$summary$signed_stat <- df2M$summary$logFC_group1 / df2M$summary$se_group1
genes2M <- setNames(df2M$summary$signed_stat, df2M$summary$gene)
genes2M <- sort(genes2M, decreasing = TRUE)

common_genes <- intersect(names(genes1M), names(genes2M))
genes1_subM <- genes1M[common_genes]
genes2_subM <- genes2M[common_genes]

genes1_dfM <- data.frame(Gene = names(genes1_subM), logFC = genes1_subM)
genes2_dfM <- data.frame(Gene = names(genes2_subM), logFC = genes2_subM)

RRHO_objM <-  RRHO2_initialize(genes1_dfM, genes2_dfM, labels = c("tsai", "dejager"), log10.ind=FALSE)
outp <- cor.test(genes1_subM, genes2_subM, method = "spearman")

png("RRHO_heatmap{variable}M.png", width=800, height=800)
RRHO2_heatmap(RRHO_objM)

lab <- sprintf("Spearman ρ = %.2f\np = %.2e",
               outp$estimate, outp$p.value)

mtext(lab,
      side = 3,        # top
      adj = 1,         # right aligned
      line = -1,       # move inside plot
      cex = 0.9)

dev.off()

"""


base_script_templateNebulaScenicDejager ="""
#Ast,'Endo', 'Ex_L2_3',
import anndata as ad
import pandas as pd
# variables = ['In_PV (Basket)', 'Ex_L4_5', 'Ex_NRGN','Ast', 'Endo', 'Ex_L2_3', 'Ex_L4', 'Ex_L5', 'Ex_L5_6', 'In_Rosehip', 'In_SST', 'In_VIP', 'Mic', 'Oli', 'OPC'
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

file_path_target = 'totalAdataAnno040825.h5ad'
adataOG = ad.read_h5ad(file_path_target)
variables = adataOG.obs['cell_type'].unique()

def fixed_micropool(adata, patient_col="patient_id", pool_size=30):
    pooled_X = []
    pooled_obs = []
    obs_cols = [c for c in adata.obs.columns]
    print(obs_cols)
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
            pooled_info = {{
                "batch": patient,
                "pool_id": patient+"_pool"+str(i)
            }}
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

x="{variable}"
cellTypeToSubset=x

if True:
    if __name__ == '__main__':
        import os
        import glob
        import pickle
        import numpy as np
        # import ACTIONet as an 
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
        #todo: put files into folder
        rootDirectory = os.path.join(os.environ["SOCISL_OUTPUT_ROOT"], "DeJager", "") 
        os.chdir(rootDirectory)
        folderPath=os.path.join(os.environ["SOCISL_OUTPUT_ROOT"], "") #change path to folder
        db_paths = [folderPath+"hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
        folderPath+"hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"]
        dbs = [RankingDatabase(fname=p, name=os.path.basename(p)) for p in db_paths]
        f_tfs = folderPath+"hg.txt"
        tf_names = load_tf_names(f_tfs)
        CSVPatient = pd.read_csv('dataset_652_basic_03-23-2022.csv')
        adata = adataOG[adataOG.obs['cell_type']==cellTypeToSubset]
        metaDF=adata.obs
        metaDF['patient_id'] = metaDF['patient_id'].astype(str)
        CSVPatient['projid'] = CSVPatient['projid'].astype(str)
        metaDF = metaDF.merge(
        CSVPatient[['projid', 'cogdx','msex']],how='left',left_on='patient_id',right_on='projid')
        metaDF.drop(columns=['projid'], inplace=True)
        adata = adata[metaDF["msex"] == 0].copy()
        adata = fixed_micropool(adata, patient_col="patient_id", pool_size=50)
        metaDF=adata.obs
        metaDF['patient_id'] = metaDF['patient_id'].astype(str)
        CSVPatient['projid'] = CSVPatient['projid'].astype(str)
        metaDF = metaDF.reset_index().merge(
        CSVPatient[['projid', 'cogdx','msex']],how='left',left_on='patient_id',right_on='projid')
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
            with open(rootDirectory+"female__regulons_Dejager"+x+".p", "wb") as f:
                pickle.dump(regulons, f)
            dfReg = pd.DataFrame (regulons, columns = ['regulon'])
            #List of enriched regulons
            dfReg.to_csv("female_regulonsFULL_Dejager"+x+".csv")
            # AUCell
            auc_mtx = aucell(df, regulons, num_workers=4)
            auc_mtx.to_csv(rootDirectory+"female_auc_mtxFULL_Dejager"+x+".csv");
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
                    data = pd.DataFrame({{'AUCell': auc_mtx[regulon],'Condition': metadata['cogdx'],'Batch': metadata['patient_id']}})
                    print(data.head)
                    plt.figure(figsize=(6,4))
                    sns.regplot(x="Condition",y="AUCell",data=data,scatter_kws={{"s": 10, "alpha": 0.4}},line_kws={{"color": "red"}})
                    plt.title(x+" AUCell vs Condition")
                    plt.xlabel("Condition (continuous)")
                    plt.ylabel("AUCell score")
                    plt.tight_layout()
                    plt.savefig(x+"_female_AUCell_vs_ConditionDejager.png", dpi=300)
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
                        results.append({{"regulon": regulon,"coef": coef,"pval": pval0}})
                        print(regulon + ":" + str(pval))
                    except Exception as e:
                        results.append({{"error": str(e)}})
                results_df = pd.DataFrame(results).dropna(subset=["coef", "pval"])
                _, pvals_fdr, _, _ = multipletests(results_df["pval"], alpha=0.05, method='fdr_bh')
                results_df["fdr"] = pvals_fdr
                results_df["log_fdr"] = -np.log10(results_df["fdr"])
                results_df["significance"] = "Not Significant"
                results_df.loc[(results_df["fdr"] < 0.05) & (results_df["coef"] > 0), "significance"] = "Up"
                results_df.loc[(results_df["fdr"] < 0.05) & (results_df["coef"] < 0), "significance"] = "Down"
                results_df.to_csv("Dejagerfemale_pVals"+x+".csv")
                palette = {{"Up": "red","Down": "blue","Not Significant": "gray"}}
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
                plt.title(f"{{x}} Social Isolation Transcription Factor Activity - Volcano")
                plt.legend(title="Significance", loc="upper right")
                plt.tight_layout()
                plt.savefig(x + "DejagerFemaleVolcano.png") #volcano plot!
                plt.show()

if True:
    if __name__ == '__main__':
        import os
        import glob
        import pickle
        import numpy as np
        # import ACTIONet as an 
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
        # patients = ["50105301", "10518782", "74284255", "10202345", "15113169", "50101659", "11157783",
        #     "10253148", "3713990", "10490993", "50106730", "50104134", "11444465", "50405330",
        #     "10394182", "50402693", "50405042", "18414513", "44299049", "10101589", "10277308",
        #     "10502798", "11327005"]
        #todo: put files into folder
        rootDirectory = os.path.join(os.environ["SOCISL_OUTPUT_ROOT"], "DeJager", "") 
        os.chdir(rootDirectory)
        folderPath=os.path.join(os.environ["SOCISL_OUTPUT_ROOT"], "") #change path to folder
        db_paths = [folderPath+"hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
        folderPath+"hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"]
        dbs = [RankingDatabase(fname=p, name=os.path.basename(p)) for p in db_paths]
        f_tfs = folderPath+"hg.txt"
        tf_names = load_tf_names(f_tfs)
        CSVPatient = pd.read_csv('dataset_652_basic_03-23-2022.csv')
        adata = adataOG[adataOG.obs['cell_type']==cellTypeToSubset]
        metaDF=adata.obs
        metaDF['patient_id'] = metaDF['patient_id'].astype(str)
        CSVPatient['projid'] = CSVPatient['projid'].astype(str)
        metaDF = metaDF.merge(
        CSVPatient[['projid', 'cogdx','msex']],how='left',left_on='patient_id',right_on='projid')
        metaDF.drop(columns=['projid'], inplace=True)
        adata = adata[metaDF["msex"] == 1].copy()
        adata = fixed_micropool(adata, patient_col="patient_id", pool_size=50)
        metaDF=adata.obs
        metaDF['patient_id'] = metaDF['patient_id'].astype(str)
        CSVPatient['projid'] = CSVPatient['projid'].astype(str)
        metaDF = metaDF.reset_index().merge(
        CSVPatient[['projid', 'cogdx','msex']],how='left',left_on='patient_id',right_on='projid')
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
            with open(rootDirectory+"male__regulons_Dejager"+x+".p", "wb") as f:
                pickle.dump(regulons, f)
            dfReg = pd.DataFrame (regulons, columns = ['regulon'])
            #List of enriched regulons
            dfReg.to_csv("male_regulonsFULL_Dejager"+x+".csv")
            # AUCell
            auc_mtx = aucell(df, regulons, num_workers=4)
            auc_mtx.to_csv(rootDirectory+"male_auc_mtxFULL_Dejager"+x+".csv");
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
                    data = pd.DataFrame({{'AUCell': auc_mtx[regulon],'Condition': metadata['cogdx'],'Batch': metadata['patient_id']}})
                    print(data.head)
                    plt.figure(figsize=(6,4))
                    sns.regplot(x="Condition",y="AUCell",data=data,scatter_kws={{"s": 10, "alpha": 0.4}},line_kws={{"color": "red"}})
                    plt.title(x+" AUCell vs Condition")
                    plt.xlabel("Condition (continuous)")
                    plt.ylabel("AUCell score")
                    plt.tight_layout()
                    plt.savefig(x+"_male_AUCell_vs_ConditionDejager.png", dpi=300)
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
                        results.append({{"regulon": regulon,"coef": coef,"pval": pval0}})
                        print(regulon + ":" + str(pval))
                    except Exception as e:
                        results.append({{"error": str(e)}})
                results_df = pd.DataFrame(results).dropna(subset=["coef", "pval"])
                _, pvals_fdr, _, _ = multipletests(results_df["pval"], alpha=0.05, method='fdr_bh')
                results_df["fdr"] = pvals_fdr
                results_df["log_fdr"] = -np.log10(results_df["fdr"])
                results_df["significance"] = "Not Significant"
                results_df.loc[(results_df["fdr"] < 0.05) & (results_df["coef"] > 0), "significance"] = "Up"
                results_df.loc[(results_df["fdr"] < 0.05) & (results_df["coef"] < 0), "significance"] = "Down"
                results_df.to_csv("Dejagermale_pVals"+x+".csv")
                palette = {{"Up": "red","Down": "blue","Not Significant": "gray"}}
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
                plt.title(f"{{x}} Social Isolation Transcription Factor Activity - Volcano")
                plt.legend(title="Significance", loc="upper right")
                plt.tight_layout()
                plt.savefig(x + "DejagerMaleVolcano.png") #volcano plot!
                plt.show()
"""

base_script_templateNebulaScenic ="""

file_path_target = 'totalAdataAnno012125.h5ad'
# variables = ['In_PV (Basket)', 'In_PV (Chandelier)']
#, 'Ex_L4_5''Ex_L5/6-CC', 'Ex_NRGN','Ast', 'Endo', 'Ex_L2_3', 'Ex_L4', 'Ex_L5', 'Ex_L5_6', 'In_Rosehip', 'In_SST', 'In_VIP', 'Mic', 'Oli', 'OPC'
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
from statsmodels.stats.multitest import multipletests

import statsmodels.formula.api as smf
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
adataOG = ad.read_h5ad(file_path_target)

def fixed_micropool(adata, patient_col="patient_id", pool_size=30):
    pooled_X = []
    pooled_obs = []
    obs_cols = [c for c in adata.obs.columns if c != patient_col]
    print(obs_cols)
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
            pooled_info = {{
                "batch": patient,
                "pool_id": patient+"_pool"+str(i)
            }}
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

x="{variable}"
cellTypeToSubset=x

if True:
    if __name__ == '__main__':
        import os
        import glob
        import pickle
        import numpy as np
        # import ACTIONet as an 
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
        # patients = ["50105301", "10518782", "74284255", "10202345", "15113169", "50101659", "11157783",
        #     "10253148", "3713990", "10490993", "50106730", "50104134", "11444465", "50405330",
        #     "10394182", "50402693", "50405042", "18414513", "44299049", "10101589", "10277308",
        #     "10502798", "11327005"]
        #todo: put files into folder
        rootDirectory = os.path.join(os.environ["SOCISL_OUTPUT_ROOT"], "DeJager", "") 
        os.chdir(rootDirectory)
        folderPath=os.path.join(os.environ["SOCISL_OUTPUT_ROOT"], "") #change path to folder
        db_paths = [folderPath+"hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
        folderPath+"hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"]
        dbs = [RankingDatabase(fname=p, name=os.path.basename(p)) for p in db_paths]
        f_tfs = folderPath+"hg.txt"
        tf_names = load_tf_names(f_tfs)
        CSVPatient = pd.read_csv('dataset_652_basic_03-23-2022.csv')
        adata = adataOG[adataOG.obs['cell_type']==cellTypeToSubset]
        metaDF=adata.obs
        metaDF['batch'] = metaDF['batch'].astype(str)
        CSVPatient['projid'] = CSVPatient['projid'].astype(str)
        metaDF = metaDF.merge(
        CSVPatient[['projid', 'cogdx','msex']],how='left',left_on='batch',right_on='projid')
        metaDF.drop(columns=['projid'], inplace=True)
        adata = adata[metaDF["msex"] == 0].copy()
        adata = fixed_micropool(adata, patient_col="batch", pool_size=50)
        metaDF=adata.obs
        metaDF['batch'] = metaDF['batch'].astype(str)
        CSVPatient['projid'] = CSVPatient['projid'].astype(str)
        metaDF = metaDF.reset_index().merge(
        CSVPatient[['projid', 'cogdx','msex']],how='left',left_on='batch',right_on='projid')
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
                    data = pd.DataFrame({{'AUCell': auc_mtx[regulon],'Condition': metadata['cogdx'],'Batch': metadata['batch']}})
                    print(data.head)
                    plt.figure(figsize=(6,4))
                    sns.regplot(x="Condition",y="AUCell",data=data,scatter_kws={{"s": 10, "alpha": 0.4}},line_kws={{"color": "red"}})
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
                        results.append({{"regulon": regulon,"coef": coef,"pval": pval0}})
                        print(regulon + ":" + str(pval))
                    except Exception as e:
                        results.append({{"error": str(e)}})
                results_df = pd.DataFrame(results).dropna(subset=["coef", "pval"])
                _, pvals_fdr, _, _ = multipletests(results_df["pval"], alpha=0.05, method='fdr_bh')
                results_df["fdr"] = pvals_fdr
                results_df["log_fdr"] = -np.log10(results_df["fdr"])
                results_df["significance"] = "Not Significant"
                results_df.loc[(results_df["fdr"] < 0.05) & (results_df["coef"] > 0), "significance"] = "Down"
                results_df.loc[(results_df["fdr"] < 0.05) & (results_df["coef"] < 0), "significance"] = "Up"
                results_df.to_csv("female_pVals"+x+".csv")
                palette = {{"Up": "red","Down": "blue","Not Significant": "gray"}}
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
                plt.title(f"{{x}} Social Isolation Transcription Factor Activity - Volcano")
                plt.legend(title="Significance", loc="upper right")
                plt.tight_layout()
                plt.savefig(x + "TsaiFemaleVolcano.png") #volcano plot!
                plt.show()

if True:
    if __name__ == '__main__':
        import os
        import glob
        import pickle
        import numpy as np
        # import ACTIONet as an 
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
        # patients = ["50105301", "10518782", "74284255", "10202345", "15113169", "50101659", "11157783",
        #     "10253148", "3713990", "10490993", "50106730", "50104134", "11444465", "50405330",
        #     "10394182", "50402693", "50405042", "18414513", "44299049", "10101589", "10277308",
        #     "10502798", "11327005"]
        #todo: put files into folder
        rootDirectory = os.path.join(os.environ["SOCISL_OUTPUT_ROOT"], "DeJager", "") 
        os.chdir(rootDirectory)
        folderPath=os.path.join(os.environ["SOCISL_OUTPUT_ROOT"], "") #change path to folder
        db_paths = [folderPath+"hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
        folderPath+"hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"]
        dbs = [RankingDatabase(fname=p, name=os.path.basename(p)) for p in db_paths]
        f_tfs = folderPath+"hg.txt"
        tf_names = load_tf_names(f_tfs)
        CSVPatient = pd.read_csv('dataset_652_basic_03-23-2022.csv')
        adata = adataOG[adataOG.obs['cell_type']==cellTypeToSubset]
        metaDF=adata.obs
        metaDF['batch'] = metaDF['batch'].astype(str)
        CSVPatient['projid'] = CSVPatient['projid'].astype(str)
        metaDF = metaDF.merge(
        CSVPatient[['projid', 'cogdx','msex']],how='left',left_on='batch',right_on='projid')
        metaDF.drop(columns=['projid'], inplace=True)
        adata = adata[metaDF["msex"] == 1].copy()
        adata = fixed_micropool(adata, patient_col="batch", pool_size=50)
        metaDF=adata.obs
        metaDF['batch'] = metaDF['batch'].astype(str)
        CSVPatient['projid'] = CSVPatient['projid'].astype(str)
        metaDF = metaDF.reset_index().merge(
        CSVPatient[['projid', 'cogdx','msex']],how='left',left_on='batch',right_on='projid')
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
                    data = pd.DataFrame({{'AUCell': auc_mtx[regulon],'Condition': metadata['cogdx'],'Batch': metadata['batch']}})
                    print(data.head)
                    plt.figure(figsize=(6,4))
                    sns.regplot(x="Condition",y="AUCell",data=data,scatter_kws={{"s": 10, "alpha": 0.4}},line_kws={{"color": "red"}})
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
                        results.append({{"regulon": regulon,"coef": coef,"pval": pval0}})
                        print(regulon + ":" + str(pval))
                    except Exception as e:
                        results.append({{"error": str(e)}})
                results_df = pd.DataFrame(results).dropna(subset=["coef", "pval"])
                _, pvals_fdr, _, _ = multipletests(results_df["pval"], alpha=0.05, method='fdr_bh')
                results_df["fdr"] = pvals_fdr
                results_df["log_fdr"] = -np.log10(results_df["fdr"])
                results_df["significance"] = "Not Significant"
                results_df.loc[(results_df["fdr"] < 0.05) & (results_df["coef"] > 0), "significance"] = "Up"
                results_df.loc[(results_df["fdr"] < 0.05) & (results_df["coef"] < 0), "significance"] = "Down"
                results_df.to_csv("male_pVals"+x+".csv")

                palette = {{"Up": "red","Down": "blue","Not Significant": "gray"}}
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
                plt.title(f"{{x}} Social Isolation Transcription Factor Activity - Volcano")
                plt.legend(title="Significance", loc="upper right")
                plt.tight_layout()
                plt.savefig(x + "TsaiMaleVolcano.png") #volcano plot!
                plt.show()
"""


base_script_templateNebRVS ="""#!/bin/bash

#SBATCH -n 45                    # Number of cores requested
#SBATCH -t 5:00:00                # Runtime in hours
#SBATCH --mem=400G              # GB memory needed (memory PER CORE)
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-user=nkhera@college.harvard.edu
#SBATCH --mail-type=ALL


source "${CONDA_INIT_SCRIPT}"

activate_env "${NEBULA_ENV}"

export HDF5_USE_FILE_LOCKING=FALSE

cd "${SOCISL_OUTPUT_ROOT}/DeJager"

Rscript rvCalcNebulaF{variable}.Rscript

"""
base_script_templateNebRVDS ="""#!/bin/bash

#SBATCH -n 45                    # Number of cores requested
#SBATCH -t 5:00:00                # Runtime in hours
#SBATCH --mem=400G              # GB memory needed (memory PER CORE)
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-user=nkhera@college.harvard.edu
#SBATCH --mail-type=ALL


source "${CONDA_INIT_SCRIPT}"

activate_env "${NEBULA_ENV}"

export HDF5_USE_FILE_LOCKING=FALSE

cd "${SOCISL_OUTPUT_ROOT}/DeJager"

Rscript rvCalcNebulaDF{variable}.Rscript

"""


base_script_templateNebFDej ="""#!/bin/bash

#SBATCH -n 45                    # Number of cores requested
#SBATCH -t 5:00:00                # Runtime in hours
#SBATCH --mem=400G              # GB memory needed (memory PER CORE)
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-user=nkhera@college.harvard.edu
#SBATCH --mail-type=ALL


source "${CONDA_INIT_SCRIPT}"

activate_env "${NEBULA_ENV}"

export HDF5_USE_FILE_LOCKING=FALSE

cd "${SOCISL_OUTPUT_ROOT}/DeJager"

Rscript dej_rvCalcNebulaF{variable}.Rscript


"""



base_script_templateNebMDej ="""#!/bin/bash

#SBATCH -n 45                    # Number of cores requested
#SBATCH -t 5:00:00                # Runtime in hours
#SBATCH --mem=400G              # GB memory needed (memory PER CORE)
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-user=nkhera@college.harvard.edu
#SBATCH --mail-type=ALL


source "${CONDA_INIT_SCRIPT}"

activate_env "${NEBULA_ENV}"

export HDF5_USE_FILE_LOCKING=FALSE

cd "${SOCISL_OUTPUT_ROOT}/DeJager"

Rscript dej_rvCalcNebulaM{variable}.Rscript

sbatch dejAmale_nebulaEnrich{variable}.sh


"""
base_script_templateESM ="""#!/bin/bash

#SBATCH -n 45                    # Number of cores requested
#SBATCH -t 5:00:00                # Runtime in hours
#SBATCH --mem=400G              # GB memory needed (memory PER CORE)
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-user=nkhera@college.harvard.edu
#SBATCH --mail-type=ALL

source "${CONDA_INIT_SCRIPT}"

echo "Hello, 1!"

cd "${SOCISL_OUTPUT_ROOT}/DeJager"

echo "Hello, 2!"

activate_env "${NEBULA_ENV}"

echo "Hello, 3!"

Rscript male_nebulaEnrich{variable}.Rscript

"""

base_script_templateScenicSh ="""#!/bin/bash
#SBATCH -n 16                    # Reduce cores to reduce thread/memory pressure
#SBATCH -t 16:00:00
#SBATCH --mem=600G              # Total memory (not per core!)
#SBATCH -o  %j.out
#SBATCH -e %j.err
#SBATCH --mail-user=nkhera@college.harvard.edu
#SBATCH --mail-type=ALL

source "${CONDA_INIT_SCRIPT}"

activate_env "${SCENIC_ANALYSIS_ENV}"

cd "${SOCISL_OUTPUT_ROOT}/DeJager"

export HDF5_USE_FILE_LOCKING=FALSE

python3 nebulaScenic{variable}.py
"""

base_script_templateScenicDejagerSh ="""#!/bin/bash
#SBATCH -n 16                    # Reduce cores to reduce thread/memory pressure
#SBATCH -t 8:00:00
#SBATCH --mem=600G              # Total memory (not per core!)
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --mail-user=nkhera@college.harvard.edu
#SBATCH --mail-type=ALL

source "${CONDA_INIT_SCRIPT}"

activate_env "${SCENIC_ANALYSIS_ENV}"

cd "${SOCISL_OUTPUT_ROOT}/DeJager"

export HDF5_USE_FILE_LOCKING=FALSE

python3 nebulaScenicDejager{variable}.py
"""


base_script_templateCF ="""#!/bin/bash
#SBATCH -n 40                    # Number of cores requested
#SBATCH -t 8:00:00                # Runtime in hours
#SBATCH --mem=600G              # GB memory needed (memory PER CORE)
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-user=nkhera@college.harvard.edu
#SBATCH --mail-type=ALL

source "${CONDA_INIT_SCRIPT}"

cd "${SOCISL_OUTPUT_ROOT}/DeJager"

#python ver = 3.9

#python3 -m pip install git+https://github.com/yoseflab/Compass.git --upgrade

activate_env "${COMPASS_ANALYSIS_ENV}"

export CPLEX_STUDIO_DIR="${CPLEX_DIR}"
export PATH=$CPLEX_STUDIO_DIR/cplex/python/3.9/x86-64_linux:$PATH
export PYTHONPATH=$CPLEX_STUDIO_DIR/cplex/python/3.9/x86-64_linux:$PYTHONPATH

echo "y" | compass --data dejagerfemale_matrix{variable}.tsv --num-processes 40 --species homo_sapiens --output-dir 2DejagerCompassPF{variable}New

"""

base_script_templateCM ="""#!/bin/bash
#SBATCH -n 40                    # Number of cores requested
#SBATCH -t 8:00:00                # Runtime in hours
#SBATCH --mem=600G              # GB memory needed (memory PER CORE)
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-user=nkhera@college.harvard.edu
#SBATCH --mail-type=ALL

source "${CONDA_INIT_SCRIPT}"

cd "${SOCISL_OUTPUT_ROOT}/DeJager"

#python ver = 3.9

#python3 -m pip install git+https://github.com/yoseflab/Compass.git --upgrade

activate_env "${COMPASS_ANALYSIS_ENV}"

export CPLEX_STUDIO_DIR="${CPLEX_DIR}"
export PATH=$CPLEX_STUDIO_DIR/cplex/python/3.9/x86-64_linux:$PATH
export PYTHONPATH=$CPLEX_STUDIO_DIR/cplex/python/3.9/x86-64_linux:$PYTHONPATH

echo "y" | compass --data dejagermale_matrix{variable}.tsv --num-processes 40 --species homo_sapiens --output-dir 2DejagerCompassPM{variable}New

"""



base_script_templateCFT ="""#!/bin/bash
#SBATCH -n 40                    # Number of cores requested
#SBATCH -t 8:00:00                # Runtime in hours
#SBATCH --mem=600G              # GB memory needed (memory PER CORE)
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-user=nkhera@college.harvard.edu
#SBATCH --mail-type=ALL

source "${CONDA_INIT_SCRIPT}"

cd "${SOCISL_OUTPUT_ROOT}/DeJager"

#python ver = 3.9

#python3 -m pip install git+https://github.com/yoseflab/Compass.git --upgrade

activate_env "${COMPASS_ANALYSIS_ENV}"

export CPLEX_STUDIO_DIR="${CPLEX_DIR}"
export PATH=$CPLEX_STUDIO_DIR/cplex/python/3.9/x86-64_linux:$PATH
export PYTHONPATH=$CPLEX_STUDIO_DIR/cplex/python/3.9/x86-64_linux:$PYTHONPATH

echo "y" | compass --data tsaifemale_matrix{variable}.tsv --num-processes 40 --species homo_sapiens --output-dir 2TsaiCompassPF{variable}New

"""

base_script_templateCMT ="""#!/bin/bash
#SBATCH -n 40                    # Number of cores requested
#SBATCH -t 8:00:00                # Runtime in hours
#SBATCH --mem=600G              # GB memory needed (memory PER CORE)
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-user=nkhera@college.harvard.edu
#SBATCH --mail-type=ALL

source "${CONDA_INIT_SCRIPT}"

cd "${SOCISL_OUTPUT_ROOT}/DeJager"

#python ver = 3.9

#python3 -m pip install git+https://github.com/yoseflab/Compass.git --upgrade

activate_env "${COMPASS_ANALYSIS_ENV}"

export CPLEX_STUDIO_DIR="${CPLEX_DIR}"
export PATH=$CPLEX_STUDIO_DIR/cplex/python/3.9/x86-64_linux:$PATH
export PYTHONPATH=$CPLEX_STUDIO_DIR/cplex/python/3.9/x86-64_linux:$PYTHONPATH

echo "y" | compass --data tsaimale_matrix{variable}.tsv --num-processes 40 --species homo_sapiens --output-dir 2TsaiCompassPM{variable}New

"""

base_script_templateRS ="""#!/bin/bash

#SBATCH -n 45                    # Number of cores requested
#SBATCH -t 1:00:00                # Runtime in hours
#SBATCH --mem=400G              # GB memory needed (memory PER CORE)
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-user=nkhera@college.harvard.edu
#SBATCH --mail-type=ALL

source "${CONDA_INIT_SCRIPT}"

cd "${SOCISL_OUTPUT_ROOT}/DeJager"

activate_env "${NEBULA_ENV}"

Rscript nebulaRRHO2{variable}.Rscript

"""

base_script_templateSS ="""#!/bin/bash
#SBATCH -n 45                    # Number of cores requested
#SBATCH -t 5:00:00                # Runtime in hours
#SBATCH --mem=400G              # GB memory needed (memory PER CORE)
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-user=nkhera@college.harvard.edu
#SBATCH --mail-type=ALL


source "${CONDA_INIT_SCRIPT}"
activate_env "${SCENIC_ANALYSIS_ENV}"

cd "${SOCISL_OUTPUT_ROOT}/Tsai"

python3 nebulaScenic{variable}.py
"""



base_script_templateES ="""#!/bin/bash

#SBATCH -n 45                    # Number of cores requested
#SBATCH -t 5:00:00                # Runtime in hours
#SBATCH --mem=400G              # GB memory needed (memory PER CORE)
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-user=nkhera@college.harvard.edu
#SBATCH --mail-type=ALL

source "${CONDA_INIT_SCRIPT}"

echo "Hello, 1!"

cd "${SOCISL_OUTPUT_ROOT}/DeJager"

echo "Hello, 2!"

activate_env "${NEBULA_ENV}"

echo "Hello, 3!"

Rscript AnebulaEnrich{variable}.Rscript

"""



base_script_templateEMS ="""#!/bin/bash

#SBATCH -n 45                    # Number of cores requested
#SBATCH -t 5:00:00                # Runtime in hours
#SBATCH --mem=400G              # GB memory needed (memory PER CORE)
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-user=nkhera@college.harvard.edu
#SBATCH --mail-type=ALL

source "${CONDA_INIT_SCRIPT}"

echo "Hello, 1!"

cd "${SOCISL_OUTPUT_ROOT}/DeJager"

echo "Hello, 2!"

activate_env "${NEBULA_ENV}"

echo "Hello, 3!"

Rscript Amale_nebulaEnrich{variable}.Rscript

"""



base_script_templateEDS ="""#!/bin/bash

#SBATCH -n 45                    # Number of cores requested
#SBATCH -t 5:00:00                # Runtime in hours
#SBATCH --mem=400G              # GB memory needed (memory PER CORE)
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-user=nkhera@college.harvard.edu
#SBATCH --mail-type=ALL

source "${CONDA_INIT_SCRIPT}"

echo "Hello, 1!"

cd "${SOCISL_OUTPUT_ROOT}/DeJager"

echo "Hello, 2!"

activate_env "${NEBULA_ENV}"

echo "Hello, 3!"

Rscript dejAnebulaEnrich{variable}.Rscript

"""



base_script_templateEDMS ="""#!/bin/bash

#SBATCH -n 45                    # Number of cores requested
#SBATCH -t 5:00:00                # Runtime in hours
#SBATCH --mem=400G              # GB memory needed (memory PER CORE)
#SBATCH -o %j.out               # Standard out goes to this file
#SBATCH -e %j.err               # Standard err goes to this file
#SBATCH --mail-user=nkhera@college.harvard.edu
#SBATCH --mail-type=ALL

source "${CONDA_INIT_SCRIPT}"

echo "Hello, 1!"

cd "${SOCISL_OUTPUT_ROOT}/DeJager"

echo "Hello, 2!"

activate_env "${NEBULA_ENV}"

echo "Hello, 3!"

Rscript dejAmale_nebulaEnrich{variable}.Rscript

"""




# BiocManager::install("rhdf5")

import re

def sanitize_filename(name):
    return re.sub(r'[\\/:"*?<>|]+', '_', name)

def escape_braces_except_variable(template):
    return template.replace('{', '{{').replace('}', '}}').replace('{{variable}}', '{variable}')


def generate_scriptsCompass(variables, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    escaped_template = escape_braces_except_variable(base_script_templateCompassFD)
    for variable in variables:
        safe_variable = sanitize_filename(variable)
        script_content = escaped_template.format(variable=variable)
        script_filename = os.path.join(output_dir, f'nebula_{safe_variable}.Rscript')
        with open(script_filename, 'w') as script_file:
            script_file.write(script_content)

def generate_scriptsNebula(variables, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    escaped_template = escape_braces_except_variable(base_script_templateNebula)
    for variable in variables:
        safe_variable = sanitize_filename(variable)
        script_content = escaped_template.format(variable=variable)
        script_filename = os.path.join(output_dir, f'nebula_{safe_variable}.Rscript')
        with open(script_filename, 'w') as script_file:
            script_file.write(script_content)

def generate_scriptsNebulaRVF(variables, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    escaped_template = escape_braces_except_variable(base_script_templateNebulaRVYay)
    for variable in variables:
        safe_variable = sanitize_filename(variable)
        script_content = escaped_template.format(variable=variable)
        script_filename = os.path.join(output_dir, f'rvCalcNebulaF{safe_variable}.Rscript')
        with open(script_filename, 'w') as script_file:
            script_file.write(script_content)

def generate_scriptsNebulaRVM(variables, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    escaped_template = escape_braces_except_variable(base_script_templateNebulaRVYayM)
    for variable in variables:
        safe_variable = sanitize_filename(variable)
        script_content = escaped_template.format(variable=variable)
        script_filename = os.path.join(output_dir, f'rvCalcNebulaM{safe_variable}.Rscript')
        with open(script_filename, 'w') as script_file:
            script_file.write(script_content)


def generate_scriptsNebulaRVDejF(variables, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    escaped_template = escape_braces_except_variable(base_scriptTemplateFDejN)
    for variable in variables:
        safe_variable = sanitize_filename(variable)
        script_content = escaped_template.format(variable=variable)
        script_filename = os.path.join(output_dir, f'dej_rvCalcNebulaF{safe_variable}.Rscript')
        with open(script_filename, 'w') as script_file:
            script_file.write(script_content)


def generate_scriptsNebulaRVDejM(variables, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    escaped_template = escape_braces_except_variable(base_scriptTemplateMDejN)
    for variable in variables:
        safe_variable = sanitize_filename(variable)
        script_content = escaped_template.format(variable=variable)
        script_filename = os.path.join(output_dir, f'dej_rvCalcNebulaM{safe_variable}.Rscript')
        with open(script_filename, 'w') as script_file:
            script_file.write(script_content)


def generate_scriptsNebulaRVDejFSh(variables, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    escaped_template = escape_braces_except_variable(base_script_templateNebFDej)
    for variable in variables:
        safe_variable = sanitize_filename(variable)
        script_content = escaped_template.format(variable=variable)
        script_filename = os.path.join(output_dir, f'adejFRVscriptNebula_{safe_variable}.sh')
        with open(script_filename, 'w') as script_file:
            script_file.write(script_content)


def generate_scriptsNebulaRVDejMSh(variables, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    escaped_template = escape_braces_except_variable(base_script_templateNebMDej)
    for variable in variables:
        safe_variable = sanitize_filename(variable)
        script_content = escaped_template.format(variable=variable)
        script_filename = os.path.join(output_dir, f'adejMRVscriptNebula_{safe_variable}.sh')
        with open(script_filename, 'w') as script_file:
            script_file.write(script_content)

# def generate_scriptsNebulaRV(variables, output_dir):
#     os.makedirs(output_dir, exist_ok=True)
#     escaped_template = escape_braces_except_variable(base_script_templateNebulaRVYay)
#     for variable in variables:
#         safe_variable = sanitize_filename(variable)
#         script_content = escaped_template.format(variable=variable)
#         script_filename = os.path.join(output_dir, f'rvCalcNebulaM{safe_variable}.Rscript')
#         with open(script_filename, 'w') as script_file:
#             script_file.write(script_content)

def generate_scriptsNebulaSh(variables, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for variable in variables:
        safe_variable = sanitize_filename(variable)
        script_content = base_script_templateNebRVS.format(variable=variable)
        script_filename = os.path.join(output_dir, f'FRVscriptNebula_{safe_variable}.sh')
        with open(script_filename, 'w') as script_file:
            script_file.write(script_content)

def generate_scriptsNebulaShD(variables, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for variable in variables:
        safe_variable = sanitize_filename(variable)
        script_content = base_script_templateNebRVDS.format(variable=variable)
        script_filename = os.path.join(output_dir, f'FRVscriptNebula_{safe_variable}D.sh')
        with open(script_filename, 'w') as script_file:
            script_file.write(script_content)

def generate_scriptsCompassShD(variables, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for variable in variables:
        safe_variable = sanitize_filename(variable)
        script_content = base_script_templateCF.format(variable=variable)
        script_filename = os.path.join(output_dir, f'dejagerFemaleCompass2{safe_variable}D.sh')
        with open(script_filename, 'w') as script_file:
            script_file.write(script_content)

def generate_scriptsCompassShDM(variables, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for variable in variables:
        safe_variable = sanitize_filename(variable)
        script_content = base_script_templateCM.format(variable=variable)
        script_filename = os.path.join(output_dir, f'dejagerMaleCompass2{safe_variable}D.sh')
        with open(script_filename, 'w') as script_file:
            script_file.write(script_content)

def generate_scriptsCompassShT(variables, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for variable in variables:
        safe_variable = sanitize_filename(variable)
        script_content = base_script_templateCFT.format(variable=variable)
        script_filename = os.path.join(output_dir, f'tsaiFemaleCompass2{safe_variable}D.sh')
        with open(script_filename, 'w') as script_file:
            script_file.write(script_content)

def generate_scriptsCompassShTM(variables, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for variable in variables:
        safe_variable = sanitize_filename(variable)
        script_content = base_script_templateCMT.format(variable=variable)
        script_filename = os.path.join(output_dir, f'tsaiMaleCompass2{safe_variable}D.sh')
        with open(script_filename, 'w') as script_file:
            script_file.write(script_content)
            
            
            

def generate_scriptsNebulaD(variables, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    escaped_template = escape_braces_except_variable(base_script_templateNebulaDejager)
    for variable in variables:
        safe_variable = sanitize_filename(variable)
        script_content = escaped_template.format(variable=variable)
        script_filename = os.path.join(output_dir, f'nebula_{safe_variable}D.Rscript')
        with open(script_filename, 'w') as script_file:
            script_file.write(script_content)


def generate_scriptsNebulaM(variables, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    escaped_template = escape_braces_except_variable(base_script_templateNebulaM)
    for variable in variables:
        safe_variable = sanitize_filename(variable)
        script_content = escaped_template.format(variable=variable)
        script_filename = os.path.join(output_dir, f'male_nebula_{safe_variable}.Rscript')
        with open(script_filename, 'w') as script_file:
            script_file.write(script_content)

def generate_scriptsNebulaShM(variables, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for variable in variables:
        safe_variable = sanitize_filename(variable)
        script_content = base_script_templateBAMM.format(variable=variable)
        script_filename = os.path.join(output_dir, f'male_scriptNebula_{safe_variable}.sh')
        with open(script_filename, 'w') as script_file:
            script_file.write(script_content)

def generate_scriptsNebulaShDM(variables, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for variable in variables:
        safe_variable = sanitize_filename(variable)
        script_content = base_script_templateBAMDejagerM.format(variable=variable)
        script_filename = os.path.join(output_dir, f'male_scriptNebula_{safe_variable}D.sh')
        with open(script_filename, 'w') as script_file:
            script_file.write(script_content)

def generate_scriptsNebulaDM(variables, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    escaped_template = escape_braces_except_variable(base_script_templateNebulaDejagerM)
    for variable in variables:
        safe_variable = sanitize_filename(variable)
        script_content = escaped_template.format(variable=variable)
        script_filename = os.path.join(output_dir, f'male_nebula_{safe_variable}D.Rscript')
        with open(script_filename, 'w') as script_file:
            script_file.write(script_content)

def generate_scriptsNebulaE(variables, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for variable in variables:
        safe_variable = sanitize_filename(variable)
        script_content = base_script_templateNebulaEnrichment.format(variable=variable)
        script_filename = os.path.join(output_dir, f'AnebulaEnrich{safe_variable}.Rscript')
        with open(script_filename, 'w') as script_file:
            script_file.write(script_content)

def generate_scriptsNebulaEM(variables, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for variable in variables:
        safe_variable = sanitize_filename(variable)
        script_content = base_script_templateNebulaEnrichmentM.format(variable=variable)
        script_filename = os.path.join(output_dir, f'Amale_nebulaEnrich{safe_variable}.Rscript')
        with open(script_filename, 'w') as script_file:
            script_file.write(script_content)


def generate_scriptsNebulaED(variables, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for variable in variables:
        safe_variable = sanitize_filename(variable)
        script_content = base_script_templateNebulaEnrichmentD.format(variable=variable)
        script_filename = os.path.join(output_dir, f'dejAnebulaEnrich{safe_variable}.Rscript')
        with open(script_filename, 'w') as script_file:
            script_file.write(script_content)

def generate_scriptsNebulaEDM(variables, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for variable in variables:
        safe_variable = sanitize_filename(variable)
        script_content = base_script_templateNebulaEnrichmentMD.format(variable=variable)
        script_filename = os.path.join(output_dir, f'dejAmale_nebulaEnrich{safe_variable}.Rscript')
        with open(script_filename, 'w') as script_file:
            script_file.write(script_content)


def generate_scriptsNebulaES(variables, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for variable in variables:
        safe_variable = sanitize_filename(variable)
        script_content = base_script_templateES.format(variable=variable)
        script_filename = os.path.join(output_dir, f'AnebulaEnrich{safe_variable}.sh')
        with open(script_filename, 'w') as script_file:
            script_file.write(script_content)

def generate_scriptsNebulaEMS(variables, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for variable in variables:
        safe_variable = sanitize_filename(variable)
        script_content = base_script_templateEMS.format(variable=variable)
        script_filename = os.path.join(output_dir, f'Amale_nebulaEnrich{safe_variable}.sh')
        with open(script_filename, 'w') as script_file:
            script_file.write(script_content)


def generate_scriptsNebulaEDS(variables, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for variable in variables:
        safe_variable = sanitize_filename(variable)
        script_content = base_script_templateEDS.format(variable=variable)
        script_filename = os.path.join(output_dir, f'dejAnebulaEnrich{safe_variable}.sh')
        with open(script_filename, 'w') as script_file:
            script_file.write(script_content)

def generate_scriptsNebulaEDMS(variables, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for variable in variables:
        safe_variable = sanitize_filename(variable)
        script_content = base_script_templateEDMS.format(variable=variable)
        script_filename = os.path.join(output_dir, f'dejAmale_nebulaEnrich{safe_variable}.sh')
        with open(script_filename, 'w') as script_file:
            script_file.write(script_content)

def generate_scriptsNebulaS(variables, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for variable in variables:
        safe_variable = sanitize_filename(variable)
        script_content = base_script_templateNebulaScenic.format(variable=variable)
        script_filename = os.path.join(output_dir, f'nebulaScenic{safe_variable}.py')
        with open(script_filename, 'w') as script_file:
            script_file.write(script_content)

def generate_scriptsNebulaSD(variables, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for variable in variables:
        safe_variable = sanitize_filename(variable)
        script_content = base_script_templateNebulaScenicDejager.format(variable=variable)
        script_filename = os.path.join(output_dir, f'nebulaScenicDejager{safe_variable}.py')
        with open(script_filename, 'w') as script_file:
            script_file.write(script_content)

def generate_scriptsNebulaSSh(variables, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for variable in variables:
        safe_variable = sanitize_filename(variable)
        script_content = base_script_templateScenicSh.format(variable=variable)
        script_filename = os.path.join(output_dir, f'nebulaScenic{safe_variable}.sh')
        with open(script_filename, 'w') as script_file:
            script_file.write(script_content)

def generate_scriptsNebulaSDSh(variables, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for variable in variables:
        safe_variable = sanitize_filename(variable)
        script_content = base_script_templateScenicDejagerSh.format(variable=variable)
        script_filename = os.path.join(output_dir, f'DnebulaScenicDejager{safe_variable}.sh')
        with open(script_filename, 'w') as script_file:
            script_file.write(script_content)


def generate_scriptsNebulaSS(variables, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for variable in variables:
        safe_variable = sanitize_filename(variable)
        script_content = base_script_templateSS.format(variable=safe_variable)
        script_filename = os.path.join(output_dir, f'nebulaScenic{safe_variable}.sh')
        with open(script_filename, 'w') as script_file:
            script_file.write(script_content)

def generate_scriptsNebulaR(variables, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for variable in variables:
        safe_variable = sanitize_filename(variable)
        script_content = base_script_templateNebulaRRHO2.format(variable=variable)
        script_filename = os.path.join(output_dir, f'nebulaRRHO2{safe_variable}.Rscript')
        with open(script_filename, 'w') as script_file:
            script_file.write(script_content)

def generate_scriptsNebulaRS(variables, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    for variable in variables:
        safe_variable = sanitize_filename(variable)
        script_content = base_script_templateRS.format(variable=safe_variable)
        script_filename = os.path.join(output_dir, f'nebulaRRHO2{safe_variable}.sh')
        with open(script_filename, 'w') as script_file:
            script_file.write(script_content)
# def generate_scriptsNebulaSS(variables, output_dir):
#     os.makedirs(output_dir, exist_ok=True)
#     for variable in variables:
#         safe_variable = sanitize_filename(variable)
#         script_content = base_script_templateNebulaScenic.format(variable=variable)
#         script_filename = os.path.join(output_dir, f'nebulaScenic{safe_variable}.sh')
#         with open(script_filename, 'w') as script_file:
#             script_file.write(script_content)


# def generate_scripts1Accompaniment(variables, output_dir):
#     os.makedirs(output_dir, exist_ok=True)

#     for variable in variables:
#         script_content = base_script_template1Accompaniment.format(variable=variable)
#         script_lines = script_content.split('\n')
#         script_content_without_first_line = '\n'.join(script_lines[1:])
#         script_filename = os.path.join(output_dir, f'script1Run_{variable}.sh')

#         with open(script_filename, 'w') as script_file:
#             script_file.write(script_content_without_first_line)



# def generate_scripts1Accompaniment(variables, output_dir):
#     os.makedirs(output_dir, exist_ok=True)

#     for variable in variables:
#         script_content = base_script_template1Accompaniment.format(variable=variable)
#         script_filename = os.path.join(output_dir, f'script1Run_{variable}.sh')

#         # Split script_content into lines and find the "!benbash" line
#         lines = script_content.split('\n')
#         start_writing = False
#         with open(script_filename, 'w') as script_file:
#             for line in lines:
#                 if start_writing:
#                     script_file.write(line + '\n')
#                 elif line.startswith('#!/bin/bash'):
#                     start_writing = True
#                     script_file.write(line + '\n')  # Include the "!benbash" line as well


# def generate_scripts2Accompaniment(variables, output_dir):
#     os.makedirs(output_dir, exist_ok=True)

#     for variable in variables:
#         script_content = base_script_template2Accompaniment.format(variable=variable)
#         script_filename = os.path.join(output_dir, f'script2Run_{variable}.sh')

#         with open(script_filename, 'w') as script_file:
#             script_file.write(script_content)

# def generate_scripts3Accompaniment(variables, output_dir):
#     os.makedirs(output_dir, exist_ok=True)

#     for variable in variables:
#         script_content = base_script_template3Accompaniment.format(variable=variable)
#         script_filename = os.path.join(output_dir, f'script3Run_{variable}.sh')

#         with open(script_filename, 'w') as script_file:
#             script_file.write(script_content)

# def generate_scripts3(variables, output_dir):
#     os.makedirs(output_dir, exist_ok=True)

#     for variable in variables:
#         script_content = base_script_template3.format(variable=variable)
#         script_filename = os.path.join(output_dir, f'script3_{variable}.py')

#         with open(script_filename, 'w') as script_file:
#             script_file.write(script_content)

# List of variables to iterate over
# source "${CONDA_INIT_SCRIPT}"


# cd "${DEJAGER_FASTQS}"

# activate_env "${BCFTOOLS_ENV}"


import os

# python3 notebookBaseCreate.py
#
output_dir = os.path.join(os.environ['SOCISL_OUTPUT_ROOT'], 'DeJager', '')

variables = ['Ast', 'Endo', 'Ex_L2_3', 'Ex_L4', 'Ex_L4_5', 'Ex_L5', 'Ex_L5_6', 'Ex_L5/6-CC', 'Ex_NRGN', 'In_PV(Basket)', 'In_PV(Chandelier)', 'In_Rosehip', 'In_SST', 'In_VIP', 'Mic', 'Oli', 'OPC']
# generate_scriptsNebulaR(variables, output_dir)
# generate_scriptsNebulaRS(variables, output_dir)

# generate_scriptsNebulaRVDejMSh(variables, output_dir)
# generate_scriptsNebulaRVDejFSh(variables, output_dir)
# generate_scriptsNebulaSh(variables, output_dir)
# generate_scriptsNebulaShD(variables, output_dir)
generate_scriptsCompassShT(variables, output_dir)
generate_scriptsCompassShTM(variables, output_dir)
# generate_scriptsNebulaRVF(variables, output_dir)
# generate_scriptsNebulaRVM(variables, output_dir)
# generate_scriptsNebulaRVDejF(variables, output_dir)
# generate_scriptsNebulaRVDejM(variables, output_dir)

# print("Hi")
# generate_scriptsNebulaES(variables, output_dir)
# print("Hi2")
# generate_scriptsNebulaE(variables, output_dir)
# print("Hi3")
# generate_scriptsNebulaESM(variables, output_dir)
# print("Hi2")
# generate_scriptsNebulaEM(variables, output_dir)
# print("Hi3")

# Directory to save the generated scripts
# os.makedirs(output_dir)

# Generate the scripts
# generate_scriptsNebulaD(variables, output_dir)
# generate_scriptsNebulaSh(variables, output_dir)
# generate_scriptsNebulaRV(variables, output_dir)
# generate_scriptsNebulaShD(variables, output_dir)
# generate_scriptsNebulaRVD(variables, output_dir)
# generate_scriptsNebulaDM(variables, output_dir)
# generate_scriptsNebulaM(variables, output_dir)
# generate_scriptsNebulaR(variables, output_dir)
# generate_scriptsNebulaRS(variables, output_dir)
# generate_scriptsNebulaE(variables, output_dir)
# generate_scriptsNebulaEM(variables, output_dir)
# generate_scriptsNebulaES(variables, output_dir)
# generate_scriptsNebulaEDS(variables, output_dir)
# # generate_scriptsNebulaEMS(variables, output_dir)
# generate_scriptsNebulaEDMS(variables, output_dir)
# # generate_scriptsNebulaE(variables, output_dir)
# generate_scriptsNebulaED(variables, output_dir)
# # generate_scriptsNebulaEM(variables, output_dir)
# generate_scriptsNebulaEDM(variables, output_dir)
# generate_scriptsNebulaE(variables, output_dir)
# generate_scriptsNebulaS(variables, output_dir)
# generate_scriptsNebulaSD(variables, output_dir)
# generate_scriptsNebulaSSh(variables, output_dir)
# generate_scriptsNebulaSDSh(variables, output_dir)
# generate_scriptsNebulaES(variables, output_dir)
# generate_scriptsNebulaESM(variables, output_dir)
# generate_scriptsNebulaRVF(variables, output_dir)
# generate_scriptsNebulaRVM(variables, output_dir)
# generate_scriptsNebulaS(variables, output_dir)
# generate_scriptsNebulaSD(variables, output_dir)
# generate_scriptsNebulaSSh(variables, output_dir)
# generate_scriptsNebulaSDSh(variables, output_dir)
# # generate_scriptsNebulaRVDejF(variables, output_dir)
# generate_scriptsNebulaRVDejM(variables, output_dir)
# generate_scriptsNebulaRVDejFSh(variables, output_dir)
# generate_scriptsNebulaRVDejMSh(variables, output_dir)
# generate_scriptsCompassShD(variables, output_dir)
# generate_scriptsCompassShDM(variables, output_dir)
# generate_scriptsCompassShT(variables, output_dir)
# generate_scriptsCompassShTM(variables, output_dir)
# generate_scriptsFreemuxNP(libraries, output_dir)
# generate_scriptsBAM(libraries, output_dir)
# generate_scriptsSample(libraries, output_dir)

# generate_scriptsDemux(libraries, output_dir)
# generate_scripts2(libraries, output_dir)
# generate_scripts1Accompaniment(libraries, output_dir)
# generate_scripts2(libraries, output_dir)

#now why dont we run all the scripts
# generate_scripts2Accompaniment(libraries, output_dir)
# generate_scripts3Accompaniment(libraries, output_dir)
