#!/usr/bin/env Rscript
# DESeq2 analysis for broad_Exc using pre-computed pseudobulk from Python.
# Reads CSV count matrices + metadata produced by pseudobulk_broad_Exc.py.
#
# Usage:
#   Rscript run_deseq_broad_Exc.R --integration projid
#   Rscript run_deseq_broad_Exc.R --integration batch

library(DESeq2)
library(dplyr)

# ── Parse arguments ──────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
integration <- "batch"

i <- 1
while (i <= length(args)) {
  if (args[i] == "--integration" && i < length(args)) {
    integration <- args[i + 1]
    i <- i + 2
  } else {
    i <- i + 1
  }
}

cat("Integration:", integration, "\n")

script_dir <- getwd()
results_base <- file.path(script_dir, paste0("results_", integration))
phenotypes <- c("tot_adverse_exp", "early_hh_ses", "ace_aggregate")

for (sex_label in c("Fem", "Male")) {
  cat("\n=== Sex:", sex_label, "===\n")

  counts_path <- file.path(results_base, paste0("pseudobulk_counts_Exc_", sex_label, ".csv"))
  meta_path <- file.path(results_base, paste0("pseudobulk_meta_Exc_", sex_label, ".csv"))

  if (!file.exists(counts_path) || !file.exists(meta_path)) {
    cat("  SKIP: pseudobulk files not found\n")
    next
  }

  # Load pseudobulk count matrix (genes x patients)
  counts <- read.csv(counts_path, row.names = 1, check.names = FALSE)
  meta <- read.csv(meta_path, stringsAsFactors = FALSE)
  meta$patient_id <- as.character(meta$patient_id)

  cat("  Loaded:", nrow(counts), "genes x", ncol(counts), "patients\n")

  # Align metadata to count matrix columns
  meta <- meta[match(colnames(counts), meta$patient_id), ]
  rownames(meta) <- meta$patient_id

  # Scale covariates (matching aceDegT.Rscript)
  meta$age_death <- as.numeric(scale(as.numeric(meta$age_death)))
  meta$pmi <- as.numeric(scale(as.numeric(meta$pmi)))

  # Save pseudobulk RDS (as a simple list, for consistency)
  for (pvar in phenotypes) {
    out_dir <- file.path(results_base, pvar)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    cat("  Phenotype:", pvar, "\n")

    # Filter patients with non-NA phenotype
    pvals <- as.numeric(meta[[pvar]])
    keep <- !is.na(pvals)
    if (sum(keep) < 5) {
      cat("    SKIP: too few non-NA values\n")
      next
    }

    counts_sub <- counts[, keep, drop = FALSE]
    meta_sub <- meta[keep, , drop = FALSE]

    # Ensure numeric
    meta_sub$niareagansc <- as.numeric(meta_sub$niareagansc)
    meta_sub[[pvar]] <- as.numeric(meta_sub[[pvar]])

    design_formula <- as.formula(paste0("~ age_death + pmi + ", pvar, " + niareagansc"))

    tryCatch({
      dds <- DESeqDataSetFromMatrix(
        countData = as.matrix(counts_sub),
        colData = meta_sub,
        design = design_formula
      )
      dds <- DESeq(dds)
      res <- results(dds, name = pvar)
      save(res, file = file.path(out_dir,
        paste0("deseqAnalysisACE_", pvar, "_Exc_", sex_label, ".rda")))
      cat("    DESeq2 done:", sum(res$padj < 0.05, na.rm = TRUE), "DEGs (padj < 0.05)\n")
    }, error = function(e) {
      cat("    ERROR:", conditionMessage(e), "\n")
    })
  }
}

cat("\nbroad_Exc DEG analysis complete.\n")
