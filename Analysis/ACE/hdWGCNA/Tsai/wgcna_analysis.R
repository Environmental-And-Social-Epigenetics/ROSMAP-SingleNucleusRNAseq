#!/usr/bin/env Rscript
# =============================================================================
# ACE hdWGCNA Co-Expression Network Analysis
#
# Identifies co-expressed gene modules at single-cell resolution using
# hdWGCNA, then tests module-trait associations with ACE phenotypes.
#
# Usage:
#   Rscript wgcna_analysis.R \
#     --cell-type Mic \
#     --sex Male \
#     --phenotype tot_adverse_exp \
#     --input-h5ad /path/to/Mic.h5ad \
#     --pheno-csv /path/to/ACE_scores.csv \
#     --output-dir /path/to/output \
#     --deg-results-dir /path/to/DEG/results \
#     [--metacell-k 25] \
#     [--smoke]
# =============================================================================

# --- Argument parsing --------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

cell_type       <- NULL
sex             <- NULL
phenotype       <- "tot_adverse_exp"
input_h5ad      <- NULL
pheno_csv       <- NULL
output_dir      <- NULL
deg_results_dir <- NULL
metacell_k      <- 25
smoke_mode      <- FALSE

i <- 1
while (i <= length(args)) {
    if (args[i] == "--cell-type") {
        i <- i + 1; cell_type <- args[i]
    } else if (args[i] == "--sex") {
        i <- i + 1; sex <- args[i]
    } else if (args[i] == "--phenotype") {
        i <- i + 1; phenotype <- args[i]
    } else if (args[i] == "--input-h5ad") {
        i <- i + 1; input_h5ad <- args[i]
    } else if (args[i] == "--pheno-csv") {
        i <- i + 1; pheno_csv <- args[i]
    } else if (args[i] == "--output-dir") {
        i <- i + 1; output_dir <- args[i]
    } else if (args[i] == "--deg-results-dir") {
        i <- i + 1; deg_results_dir <- args[i]
    } else if (args[i] == "--metacell-k") {
        i <- i + 1; metacell_k <- as.integer(args[i])
    } else if (args[i] == "--smoke") {
        smoke_mode <- TRUE
    } else {
        stop(paste("Unknown argument:", args[i]))
    }
    i <- i + 1
}

stopifnot(!is.null(cell_type), !is.null(sex), !is.null(input_h5ad),
          !is.null(pheno_csv), !is.null(output_dir))
stopifnot(sex %in% c("Male", "Female"))

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== hdWGCNA Analysis ===\n")
cat("Cell type:   ", cell_type, "\n")
cat("Sex:         ", sex, "\n")
cat("Phenotype:   ", phenotype, "\n")
cat("Input:       ", input_h5ad, "\n")
cat("Metacell k:  ", metacell_k, "\n")
cat("Smoke mode:  ", smoke_mode, "\n\n")

# --- Load libraries ----------------------------------------------------------
suppressPackageStartupMessages({
    library(Seurat)
    library(hdWGCNA)
    library(WGCNA)
    library(zellkonverter)
    library(SingleCellExperiment)
    library(dplyr)
    library(ggplot2)
})

# Allow multi-threading for WGCNA
allowWGCNAThreads()

# --- Load and prepare data ---------------------------------------------------
cat("Loading h5ad...\n")
# Use the native R reader (reader="R"), matching the DEG pipeline. The default
# (python/basilisk) reader mis-handles these large CSR splits and produces a
# dimnames/Dim mismatch ("'X' matrix could not be converted to R"). The R reader
# reads the same files correctly (as the DEG step does).
sce <- readH5AD(input_h5ad, reader = "R")

# Map sex code
sex_code <- ifelse(sex == "Male", 1, 0)

# Get metadata. Preserve the cell barcode as an explicit column: merge() below
# reorders rows and resets rownames, which would break `sce[, rownames(meta)]`.
meta <- as.data.frame(colData(sce))
meta$.barcode <- colnames(sce)
pid_col <- if ("projid" %in% colnames(meta)) {
    "projid"
} else if ("patient_id" %in% colnames(meta)) {
    "patient_id"
} else {
    "sample_id"
}
meta[[pid_col]] <- as.character(meta[[pid_col]])

# Shared ACE phenotype loader (dedups pooled Tsai+DeJager CSV before merge)
script_path <- sub("^--file=", "",
                   grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
shared_loader <- file.path(dirname(script_path), "..", "..", "_shared", "load_ace_phenotype.R")
source(shared_loader)
cat("Loading phenotype data via shared loader...\n")
pheno <- load_ace_for_projids(pheno_csv, unique(meta[[pid_col]]))
pheno$projid <- rownames(pheno)

meta <- merge(meta, pheno[, c("projid", phenotype, "msex", "age_death", "pmi",
                                "niareagansc")],
              by.x = pid_col, by.y = "projid", all.x = TRUE,
              suffixes = c("", ".pheno"))

# Filter by sex and non-NA phenotype
meta <- meta[!is.na(meta$msex) & meta$msex == sex_code, ]
meta <- meta[!is.na(meta[[phenotype]]), ]

n_patients <- length(unique(meta[[pid_col]]))
cat("Patients after filtering:", n_patients, "\n")
cat("Cells after filtering:   ", nrow(meta), "\n")

if (n_patients < 10) {
    cat("WARNING: Fewer than 10 patients. Skipping hdWGCNA.\n")
    quit(save = "no", status = 0)
}

# Subset SCE to the retained cells by barcode (merge() reset meta's rownames, so
# use the preserved .barcode column, not rownames(meta)). This also re-aligns the
# SCE column order to match `meta` row order so downstream metadata assignment by
# position (seurat_obj@meta.data[[...]] <- meta$...) stays correct.
sce <- sce[, meta$.barcode]

# In smoke mode, subsample
if (smoke_mode) {
    set.seed(42)
    n_keep <- min(5000, ncol(sce))
    keep_idx <- sample(ncol(sce), n_keep)
    sce <- sce[, keep_idx]
    meta <- meta[keep_idx, ]
    cat("SMOKE MODE: subsampled to", n_keep, "cells\n")
}

# --- Convert to Seurat -------------------------------------------------------
cat("Converting to Seurat object...\n")
seurat_obj <- as.Seurat(sce, counts = "X", data = NULL)
# as.Seurat() names the assay "originalexp" (the SCE default), but hdWGCNA's
# SetDatExpr() below expects assay = "RNA". Rename so the pipeline finds it.
if (!"RNA" %in% Assays(seurat_obj)) {
    seurat_obj <- RenameAssays(seurat_obj, assay.name = DefaultAssay(seurat_obj),
                               new.assay.name = "RNA")
    DefaultAssay(seurat_obj) <- "RNA"
}
seurat_obj@meta.data[[pid_col]] <- meta[[pid_col]]
seurat_obj@meta.data[[phenotype]] <- meta[[phenotype]]
seurat_obj@meta.data$age_death <- meta$age_death
seurat_obj@meta.data$pmi <- meta$pmi
seurat_obj@meta.data$niareagansc <- meta$niareagansc

# Normalize
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 3000)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj, npcs = 30)

# Set cell type identity (single cell type, but needed for hdWGCNA)
seurat_obj$cell_type <- cell_type
Idents(seurat_obj) <- "cell_type"

# --- hdWGCNA setup -----------------------------------------------------------
cat("Setting up hdWGCNA...\n")
seurat_obj <- SetupForWGCNA(
    seurat_obj,
    gene_select = "fraction",
    fraction = 0.05,  # genes expressed in at least 5% of cells
    wgcna_name = "ACE_WGCNA"
)

# --- Construct metacells -----------------------------------------------------
cat("Constructing metacells (k=", metacell_k, ")...\n")
seurat_obj <- MetacellsByGroups(
    seurat_obj,
    group.by = c("cell_type", pid_col),
    k = metacell_k,
    max_shared = 10,
    ident.group = "cell_type"
)

seurat_obj <- NormalizeMetacells(seurat_obj)

# --- WGCNA network construction ---------------------------------------------
cat("Selecting soft threshold power...\n")
seurat_obj <- SetDatExpr(seurat_obj, group_name = cell_type,
                          group.by = "cell_type",
                          assay = "RNA", slot = "data")

# Test soft thresholds
pdf(file.path(output_dir, "soft_threshold_plot.pdf"), width = 8, height = 5)
seurat_obj <- TestSoftPowers(seurat_obj, networkType = "signed")
plot_list <- PlotSoftPowers(seurat_obj)
print(wrap_plots(plot_list, ncol = 2))
dev.off()

cat("Constructing co-expression network...\n")
seurat_obj <- ConstructNetwork(
    seurat_obj,
    tom_name = "ACE_TOM",
    overwrite_tom = TRUE
)

# Module dendrogram
pdf(file.path(output_dir, "dendrogram.pdf"), width = 12, height = 6)
PlotDendrogram(seurat_obj, main = paste(cell_type, sex, "- Gene Dendrogram"))
dev.off()

# --- Module eigengenes -------------------------------------------------------
cat("Computing module eigengenes...\n")
seurat_obj <- ModuleEigengenes(seurat_obj, group.by.vars = pid_col)

# Get module assignments
modules <- GetModules(seurat_obj)
write.csv(modules, file.path(output_dir, "module_assignments.csv"),
          row.names = FALSE)

cat("Modules identified:", length(unique(modules$module)) - 1, "\n")  # -1 for grey

# --- Module eigengenes per patient -------------------------------------------
cat("Aggregating module eigengenes to patient level...\n")
me_df <- GetMEs(seurat_obj)
me_df[[pid_col]] <- seurat_obj@meta.data[[pid_col]]

# Average MEs across metacells per patient
me_patient <- me_df %>%
    group_by(across(all_of(pid_col))) %>%
    summarise(across(starts_with("ME"), mean, na.rm = TRUE), .groups = "drop")

# Merge with phenotype
me_patient <- merge(me_patient, pheno[, c("projid", phenotype, "age_death",
                                            "pmi", "niareagansc")],
                     by.x = pid_col, by.y = "projid")

write.csv(me_patient, file.path(output_dir, "module_eigengenes.csv"),
          row.names = FALSE)

# --- Module-trait correlations -----------------------------------------------
cat("Testing module-trait correlations...\n")
me_cols <- grep("^ME", colnames(me_patient), value = TRUE)
me_cols <- me_cols[me_cols != "MEgrey"]  # exclude grey (unassigned)

trait_results <- data.frame(
    module = character(),
    n_genes = integer(),
    correlation = numeric(),
    pvalue = numeric(),
    stringsAsFactors = FALSE
)

for (me_col in me_cols) {
    module_name <- sub("^ME", "", me_col)
    n_genes <- sum(modules$module == module_name, na.rm = TRUE)

    ct <- cor.test(me_patient[[me_col]], me_patient[[phenotype]],
                   method = "pearson")

    trait_results <- rbind(trait_results, data.frame(
        module = module_name,
        n_genes = n_genes,
        correlation = ct$estimate,
        pvalue = ct$p.value,
        stringsAsFactors = FALSE
    ))
}

# FDR correction
trait_results$padj <- p.adjust(trait_results$pvalue, method = "BH")
trait_results <- trait_results[order(trait_results$pvalue), ]

write.csv(trait_results, file.path(output_dir, "module_trait_correlations.csv"),
          row.names = FALSE)

cat("Significant modules (padj < 0.05):", sum(trait_results$padj < 0.05), "\n")

# --- Module-DEG overlap (if DEG results available) ---------------------------
if (!is.null(deg_results_dir) && dir.exists(deg_results_dir)) {
    cat("Computing module-DEG overlap...\n")

    # Find DEG results file for this cell type and sex
    sex_label <- ifelse(sex == "Male", "Male", "Fem")
    deg_cell_type <- sub("^broad_", "", cell_type)
    rda_pattern <- paste0("deseqAnalysisACE_", phenotype, "_", deg_cell_type, "_",
                          sex_label, "\\.rda$")
    rda_files <- list.files(deg_results_dir, pattern = rda_pattern,
                            full.names = TRUE)

    if (length(rda_files) > 0) {
        env <- new.env()
        load(rda_files[1], envir = env)

        # Get DESeq2 results object. Our .rda files store the result as "res"
        # (or "res_ace" for the two-contrast male arms); older fixtures used
        # resF/resM. Prefer res / res_ace, fall back to sex-specific names.
        res_obj <- if (exists("res", envir = env)) env$res
                   else if (exists("res_ace", envir = env)) env$res_ace
                   else if (sex == "Female" && exists("resF", envir = env)) env$resF
                   else if (exists("resM", envir = env)) env$resM
                   else NULL
        if (!is.null(res_obj)) {
            deg_genes <- rownames(res_obj)[!is.na(res_obj$padj) &
                                            res_obj$padj < 0.05]

            overlap_results <- data.frame(
                module = character(),
                n_module_genes = integer(),
                n_deg_in_module = integer(),
                fraction_deg = numeric(),
                enrichment_p = numeric(),
                stringsAsFactors = FALSE
            )

            all_genes <- unique(modules$gene_name)
            n_total <- length(all_genes)
            n_deg_total <- length(intersect(deg_genes, all_genes))

            for (mod in unique(modules$module)) {
                if (mod == "grey") next
                mod_genes <- modules$gene_name[modules$module == mod]
                n_mod <- length(mod_genes)
                n_overlap <- length(intersect(mod_genes, deg_genes))

                # Hypergeometric test
                hp <- phyper(n_overlap - 1, n_deg_total,
                             n_total - n_deg_total, n_mod,
                             lower.tail = FALSE)

                overlap_results <- rbind(overlap_results, data.frame(
                    module = mod,
                    n_module_genes = n_mod,
                    n_deg_in_module = n_overlap,
                    fraction_deg = n_overlap / max(n_mod, 1),
                    enrichment_p = hp,
                    stringsAsFactors = FALSE
                ))
            }

            overlap_results$enrichment_padj <- p.adjust(
                overlap_results$enrichment_p, method = "BH")
            write.csv(overlap_results,
                      file.path(output_dir, "module_deg_overlap.csv"),
                      row.names = FALSE)
        }
    }
}

# --- Module GO enrichment ----------------------------------------------------
cat("Running GO enrichment for modules...\n")
tryCatch({
    library(clusterProfiler)
    library(org.Hs.eg.db)

    go_results_all <- data.frame()

    for (mod in unique(modules$module)) {
        if (mod == "grey") next
        mod_genes <- modules$gene_name[modules$module == mod]

        tryCatch({
            ego <- enrichGO(
                gene = mod_genes,
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.1
            )

            if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
                ego_df <- as.data.frame(ego)
                ego_df$module <- mod
                go_results_all <- rbind(go_results_all, ego_df)
            }
        }, error = function(e) {
            cat("  GO enrichment failed for module", mod, ":",
                conditionMessage(e), "\n")
        })
    }

    if (nrow(go_results_all) > 0) {
        write.csv(go_results_all,
                  file.path(output_dir, "module_go_enrichment.csv"),
                  row.names = FALSE)
    }
}, error = function(e) {
    cat("GO enrichment skipped:", conditionMessage(e), "\n")
})

# --- Save WGCNA object -------------------------------------------------------
cat("Saving WGCNA object...\n")
saveRDS(seurat_obj, file.path(output_dir, "wgcna_object.rds"))

# --- Summary -----------------------------------------------------------------
cat("\n=== hdWGCNA Analysis Complete ===\n")
cat("Cell type:           ", cell_type, "\n")
cat("Sex:                 ", sex, "\n")
cat("Modules identified:  ", length(unique(modules$module)) - 1, "\n")
cat("Significant modules: ", sum(trait_results$padj < 0.05), "\n")
cat("Output:              ", output_dir, "\n")
