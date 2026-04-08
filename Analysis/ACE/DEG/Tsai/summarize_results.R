#!/usr/bin/env Rscript
# summarize_results.R -- Aggregate DESeq2 results into CSV and JSON summaries
#
# Usage:
#   Rscript summarize_results.R --results-dir results_batch --output-dir output/
#
# Outputs:
#   deg_summary.csv  -- flat table for figures and human inspection
#   deg_summary.json -- structured data for Claude agent report generation

library(dplyr)
library(jsonlite)

# ── Parse arguments ──────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
results_dir <- NULL
output_dir <- NULL

i <- 1
while (i <= length(args)) {
  if (args[i] == "--results-dir" && i < length(args)) {
    results_dir <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--output-dir" && i < length(args)) {
    output_dir <- args[i + 1]; i <- i + 2
  } else {
    i <- i + 1
  }
}

if (is.null(results_dir)) stop("ERROR: --results-dir required")
if (is.null(output_dir)) output_dir <- results_dir

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Results directory:", results_dir, "\n")
cat("Output directory:", output_dir, "\n")

# ── Walk .rda files and extract statistics ────────────────────────────────────
phenotype_dirs <- list.dirs(results_dir, recursive = FALSE, full.names = TRUE)
integration <- basename(results_dir)
# Strip "results_" prefix if present
integration <- sub("^results_", "", integration)

rows <- list()
for (pheno_dir in phenotype_dirs) {
  phenotype <- basename(pheno_dir)
  rda_files <- list.files(pheno_dir, pattern = "\\.rda$", full.names = TRUE)

  for (f in rda_files) {
    tryCatch({
      load(f)
      bn <- sub("\\.rda$", "", basename(f))
      # Pattern: deseqAnalysisACE_{phenotype}_{celltype}_{sex}
      suffix <- sub(paste0("deseqAnalysisACE_", phenotype, "_"), "", bn)
      sex <- ifelse(grepl("_Fem$", suffix), "Female", "Male")
      cell_type <- sub("_(Fem|Male)$", "", suffix)

      n_tested <- sum(!is.na(res$padj))
      n_sig <- sum(res$padj < 0.05, na.rm = TRUE)
      n_up <- sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE)
      n_down <- sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)

      # Top gene by padj
      top_gene <- ""
      top_padj <- NA
      if (n_sig > 0) {
        idx <- which.min(res$padj)
        top_gene <- rownames(res)[idx]
        top_padj <- res$padj[idx]
      }

      rows[[length(rows) + 1]] <- data.frame(
        integration = integration, phenotype = phenotype,
        cell_type = cell_type, sex = sex,
        n_tested = n_tested, n_sig = n_sig, n_up = n_up, n_down = n_down,
        top_gene = top_gene, top_padj = top_padj,
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      cat("  WARN: failed to load", f, ":", conditionMessage(e), "\n")
    })
  }
}

df <- do.call(rbind, rows)
cat("Summarized", nrow(df), "result files\n")

# ── Write CSV ─────────────────────────────────────────────────────────────────
csv_path <- file.path(output_dir, "deg_summary.csv")
write.csv(df, csv_path, row.names = FALSE)
cat("Wrote:", csv_path, "\n")

# ── Build JSON summary ────────────────────────────────────────────────────────
# Totals by phenotype x sex
totals <- df %>%
  group_by(phenotype, sex) %>%
  summarize(
    total_degs = sum(n_sig), total_up = sum(n_up), total_down = sum(n_down),
    .groups = "drop"
  ) %>%
  arrange(phenotype, sex)

# Top cell types per phenotype x sex
top_celltypes <- df %>%
  filter(n_sig > 0) %>%
  group_by(phenotype, sex) %>%
  arrange(desc(n_sig)) %>%
  slice_head(n = 10) %>%
  ungroup()

# Grand totals
grand <- list(
  integration = integration,
  total_result_files = nrow(df),
  total_degs = sum(df$n_sig),
  total_up = sum(df$n_up),
  total_down = sum(df$n_down),
  pct_down = round(sum(df$n_down) / max(sum(df$n_sig), 1) * 100, 1)
)

summary_json <- list(
  metadata = grand,
  totals_by_phenotype_sex = totals,
  top_celltypes = top_celltypes,
  full_table = df
)

json_path <- file.path(output_dir, "deg_summary.json")
write_json(summary_json, json_path, pretty = TRUE, auto_unbox = TRUE)
cat("Wrote:", json_path, "\n")
cat("Done.\n")
