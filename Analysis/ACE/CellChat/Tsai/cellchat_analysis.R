#!/usr/bin/env Rscript
# =============================================================================
# ACE Cell-Cell Communication Analysis via CellChat
#
# Compares ligand-receptor interactions between ACE-high and ACE-low groups
# using CellChat v2. Individuals are split by median of tot_adverse_exp
# within a given sex.
#
# Usage:
#   Rscript cellchat_analysis.R \
#     --sex Male \
#     --input-h5ad /path/to/tsai_annotated.h5ad \
#     --pheno-csv /path/to/ACE_scores.csv \
#     --output-dir /path/to/output \
#     [--phenotype tot_adverse_exp] \
#     [--smoke]
# =============================================================================

# --- Argument parsing --------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

sex         <- NULL
input_h5ad  <- NULL
pheno_csv   <- NULL
output_dir  <- NULL
phenotype   <- "tot_adverse_exp"
smoke_mode  <- FALSE

i <- 1
while (i <= length(args)) {
    if (args[i] == "--sex") {
        i <- i + 1; sex <- args[i]
    } else if (args[i] == "--input-h5ad") {
        i <- i + 1; input_h5ad <- args[i]
    } else if (args[i] == "--pheno-csv") {
        i <- i + 1; pheno_csv <- args[i]
    } else if (args[i] == "--output-dir") {
        i <- i + 1; output_dir <- args[i]
    } else if (args[i] == "--phenotype") {
        i <- i + 1; phenotype <- args[i]
    } else if (args[i] == "--smoke") {
        smoke_mode <- TRUE
    } else {
        stop(paste("Unknown argument:", args[i]))
    }
    i <- i + 1
}

# Validate required args
stopifnot(!is.null(sex), !is.null(input_h5ad), !is.null(pheno_csv),
          !is.null(output_dir))
stopifnot(sex %in% c("Male", "Female"))

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== CellChat Analysis ===\n")
cat("Sex:       ", sex, "\n")
cat("Phenotype: ", phenotype, "\n")
cat("Input:     ", input_h5ad, "\n")
cat("Output:    ", output_dir, "\n")
cat("Smoke:     ", smoke_mode, "\n\n")

# --- Load libraries ----------------------------------------------------------
suppressPackageStartupMessages({
    library(CellChat)
    library(zellkonverter)
    library(SingleCellExperiment)
    library(Seurat)
    library(dplyr)
    library(ggplot2)
    library(patchwork)
})

# --- Load data ---------------------------------------------------------------
cat("Loading h5ad...\n")
sce <- readH5AD(input_h5ad)

# Map sex code
sex_code <- ifelse(sex == "Male", 1, 0)

# --- Merge phenotype and filter ----------------------------------------------
meta <- as.data.frame(colData(sce))

# Identify patient ID column
pid_col <- if ("projid" %in% colnames(meta)) "projid" else if ("patient_id" %in% colnames(meta)) "patient_id" else "sample_id"
meta[[pid_col]] <- as.character(meta[[pid_col]])

# Shared ACE phenotype loader (dedups pooled Tsai+DeJager CSV before merge)
script_path <- sub("^--file=", "",
                   grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
shared_loader <- file.path(dirname(script_path), "..", "..", "_shared", "load_ace_phenotype.R")
source(shared_loader)
cat("Loading phenotype data via shared loader...\n")
pheno <- load_ace_for_projids(pheno_csv, unique(meta[[pid_col]]))
pheno$projid <- rownames(pheno)

# Merge phenotype
meta <- merge(meta, pheno[, c("projid", phenotype, "msex", "age_death", "pmi", "niareagansc")],
              by.x = pid_col, by.y = "projid", all.x = TRUE, suffixes = c("", ".pheno"))

# Filter to specified sex and non-NA phenotype
meta <- meta[!is.na(meta$msex) & meta$msex == sex_code, ]
meta <- meta[!is.na(meta[[phenotype]]), ]

if (nrow(meta) == 0) {
    cat("WARNING: No cells remaining after filtering. Exiting.\n")
    quit(save = "no", status = 0)
}

# Split by median of phenotype
median_val <- median(meta[[phenotype]], na.rm = TRUE)
meta$ace_group <- ifelse(meta[[phenotype]] > median_val, "ACE_high", "ACE_low")

cat("Median", phenotype, "=", median_val, "\n")
cat("ACE_high patients:", length(unique(meta[meta$ace_group == "ACE_high", pid_col])), "\n")
cat("ACE_low patients: ", length(unique(meta[meta$ace_group == "ACE_low", pid_col])), "\n")
cat("ACE_high cells:   ", sum(meta$ace_group == "ACE_high"), "\n")
cat("ACE_low cells:    ", sum(meta$ace_group == "ACE_low"), "\n\n")

# Subset SCE to filtered cells
sce <- sce[, rownames(meta)]

# In smoke mode, subsample to 5000 cells per group
if (smoke_mode) {
    cat("SMOKE MODE: subsampling to 5000 cells per group\n")
    set.seed(42)
    high_idx <- which(meta$ace_group == "ACE_high")
    low_idx <- which(meta$ace_group == "ACE_low")
    high_sample <- sample(high_idx, min(5000, length(high_idx)))
    low_sample <- sample(low_idx, min(5000, length(low_idx)))
    keep <- sort(c(high_sample, low_sample))
    sce <- sce[, keep]
    meta <- meta[keep, ]
}

# --- Build CellChat objects --------------------------------------------------
run_cellchat <- function(sce_sub, meta_sub, group_label) {
    cat("Building CellChat object for", group_label, "...\n")

    # Get normalized expression matrix
    expr <- as.matrix(counts(sce_sub))
    # CPM normalize
    lib_sizes <- colSums(expr)
    lib_sizes[lib_sizes == 0] <- 1
    expr <- t(t(expr) / lib_sizes) * 1e6
    expr <- log1p(expr)

    # Cell type labels
    cell_labels <- meta_sub$cell_type
    names(cell_labels) <- rownames(meta_sub)

    # Create CellChat object
    cellchat <- createCellChat(object = expr, meta = meta_sub,
                                group.by = "cell_type")

    # Set L-R database
    CellChatDB <- CellChatDB.human
    cellchat@DB <- CellChatDB

    # Preprocessing
    cat("  Identifying overexpressed genes...\n")
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)

    # Communication inference
    cat("  Computing communication probability...\n")
    cellchat <- computeCommunProb(cellchat, type = "triMean")
    cellchat <- filterCommunication(cellchat, min.cells = 10)

    # Pathway-level
    cat("  Computing pathway-level communication...\n")
    cellchat <- computeCommunProbPathway(cellchat)

    # Network analysis
    cat("  Aggregating network...\n")
    cellchat <- aggregateNet(cellchat)

    cat("  Done with", group_label, "\n\n")
    return(cellchat)
}

# Split cells by ACE group
high_cells <- rownames(meta[meta$ace_group == "ACE_high", ])
low_cells <- rownames(meta[meta$ace_group == "ACE_low", ])

cellchat_high <- run_cellchat(sce[, high_cells], meta[high_cells, ], "ACE_high")
cellchat_low <- run_cellchat(sce[, low_cells], meta[low_cells, ], "ACE_low")

# --- Save individual CellChat objects ----------------------------------------
cat("Saving CellChat objects...\n")
saveRDS(cellchat_high, file.path(output_dir, "cellchat_high.rds"))
saveRDS(cellchat_low, file.path(output_dir, "cellchat_low.rds"))

# --- Merge and compare -------------------------------------------------------
cat("Merging and comparing CellChat objects...\n")
object.list <- list(ACE_high = cellchat_high, ACE_low = cellchat_low)
cellchat_merged <- mergeCellChat(object.list, add.names = names(object.list))

# --- Extract differential interactions --------------------------------------
cat("Extracting differential interactions...\n")

# Net interaction count comparison
net_high <- cellchat_high@net$count
net_low <- cellchat_low@net$count

# Ensure same dimensions
common_ct <- intersect(rownames(net_high), rownames(net_low))
net_diff <- net_high[common_ct, common_ct] - net_low[common_ct, common_ct]

# Flatten to data frame
diff_df <- expand.grid(sender = common_ct, receiver = common_ct,
                       stringsAsFactors = FALSE)
diff_df$count_high <- as.vector(net_high[common_ct, common_ct])
diff_df$count_low <- as.vector(net_low[common_ct, common_ct])
diff_df$count_diff <- diff_df$count_high - diff_df$count_low
diff_df <- diff_df[diff_df$count_high > 0 | diff_df$count_low > 0, ]
diff_df <- diff_df[order(-abs(diff_df$count_diff)), ]

write.csv(diff_df, file.path(output_dir, "differential_interactions.csv"),
          row.names = FALSE)

# --- Pathway-level changes ---------------------------------------------------
cat("Extracting pathway-level changes...\n")

# Get pathway probabilities
path_high <- cellchat_high@netP$pathways
path_low <- cellchat_low@netP$pathways
common_paths <- intersect(path_high, path_low)

pathway_results <- data.frame(
    pathway = character(),
    prob_high = numeric(),
    prob_low = numeric(),
    prob_diff = numeric(),
    stringsAsFactors = FALSE
)

for (path in common_paths) {
    tryCatch({
        p_high <- sum(cellchat_high@netP$prob[, , path], na.rm = TRUE)
        p_low <- sum(cellchat_low@netP$prob[, , path], na.rm = TRUE)
        pathway_results <- rbind(pathway_results, data.frame(
            pathway = path,
            prob_high = p_high,
            prob_low = p_low,
            prob_diff = p_high - p_low,
            stringsAsFactors = FALSE
        ))
    }, error = function(e) {
        cat("  Skipping pathway", path, ":", conditionMessage(e), "\n")
    })
}

pathway_results <- pathway_results[order(-abs(pathway_results$prob_diff)), ]
write.csv(pathway_results, file.path(output_dir, "pathway_changes.csv"),
          row.names = FALSE)

# --- Focus axes analysis -----------------------------------------------------
cat("Analyzing focused interaction axes...\n")

focus_axes <- list(
    list(sender = "Mic", receiver = "In-PV_Basket", label = "Mic_to_PVBasket"),
    list(sender = "Mic", receiver = "Exc", label = "Mic_to_Exc"),
    list(sender = "Mic", receiver = "Inh", label = "Mic_to_Inh"),
    list(sender = "Ast", receiver = "In-PV_Basket", label = "Ast_to_PVBasket"),
    list(sender = "Ast", receiver = "Exc", label = "Ast_to_Exc"),
    list(sender = "Oli", receiver = "Exc", label = "Oli_to_Exc"),
    list(sender = "Exc", receiver = "Mic", label = "Exc_to_Mic"),
    list(sender = "In-PV_Basket", receiver = "Mic", label = "PVBasket_to_Mic")
)

focus_results <- data.frame(
    axis = character(),
    sender = character(),
    receiver = character(),
    n_interactions_high = numeric(),
    n_interactions_low = numeric(),
    n_diff = numeric(),
    stringsAsFactors = FALSE
)

for (ax in focus_axes) {
    if (ax$sender %in% common_ct && ax$receiver %in% common_ct) {
        n_high <- net_high[ax$sender, ax$receiver]
        n_low <- net_low[ax$sender, ax$receiver]
        focus_results <- rbind(focus_results, data.frame(
            axis = ax$label,
            sender = ax$sender,
            receiver = ax$receiver,
            n_interactions_high = n_high,
            n_interactions_low = n_low,
            n_diff = n_high - n_low,
            stringsAsFactors = FALSE
        ))
    }
}

write.csv(focus_results, file.path(output_dir, "focus_axes_results.csv"),
          row.names = FALSE)

# --- Save merged object for visualization ------------------------------------
saveRDS(cellchat_merged, file.path(output_dir, "cellchat_merged.rds"))

cat("\n=== CellChat analysis complete ===\n")
cat("Differential interactions:", nrow(diff_df), "cell-type pairs\n")
cat("Pathway changes:         ", nrow(pathway_results), "pathways\n")
cat("Focus axes:              ", nrow(focus_results), "axes\n")
cat("Output directory:        ", output_dir, "\n")
