#!/usr/bin/env Rscript
# run_go_enrichment.R -- Automated GO enrichment for top cell types
#
# Reads deg_summary.json to auto-select which cell types deserve GO analysis,
# then runs clusterProfiler enrichGO for each.
#
# Usage:
#   Rscript run_go_enrichment.R \
#     --results-dir results_batch \
#     --summary-json results_batch/deg_summary.json \
#     --output-dir results_batch/go_enrichment \
#     --min-degs 50 --top-n 5

library(dplyr)
library(jsonlite)
library(clusterProfiler)
library(org.Hs.eg.db)

# ── Parse arguments ──────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
results_dir <- NULL
summary_json_path <- NULL
output_dir <- NULL
min_degs <- 50
top_n <- 5

i <- 1
while (i <= length(args)) {
  if (args[i] == "--results-dir" && i < length(args)) {
    results_dir <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--summary-json" && i < length(args)) {
    summary_json_path <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--output-dir" && i < length(args)) {
    output_dir <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--min-degs" && i < length(args)) {
    min_degs <- as.integer(args[i + 1]); i <- i + 2
  } else if (args[i] == "--top-n" && i < length(args)) {
    top_n <- as.integer(args[i + 1]); i <- i + 2
  } else {
    i <- i + 1
  }
}

if (is.null(results_dir)) stop("ERROR: --results-dir required")
if (is.null(summary_json_path)) {
  summary_json_path <- file.path(results_dir, "deg_summary.json")
}
if (is.null(output_dir)) output_dir <- file.path(results_dir, "go_enrichment")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
fig_dir <- file.path(output_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE)

cat("Results dir:", results_dir, "\n")
cat("Summary JSON:", summary_json_path, "\n")
cat("Output dir:", output_dir, "\n")
cat("Min DEGs:", min_degs, "| Top N:", top_n, "\n\n")

# ── Read summary and select cell types ────────────────────────────────────────
summary <- fromJSON(summary_json_path)
full <- as.data.frame(summary$full_table)

# Select top cell types per phenotype x sex x direction
targets <- full %>%
  filter(n_sig >= min_degs) %>%
  group_by(phenotype, sex) %>%
  arrange(desc(n_sig)) %>%
  slice_head(n = top_n) %>%
  ungroup()

cat("Selected", nrow(targets), "cell type x phenotype x sex combinations for GO:\n")
print(targets[, c("phenotype", "cell_type", "sex", "n_sig", "n_up", "n_down")])
cat("\n")

# ── Run GO enrichment ─────────────────────────────────────────────────────────
go_results_all <- list()

for (r in seq_len(nrow(targets))) {
  row <- targets[r, ]
  ct <- row$cell_type
  sex <- row$sex
  pheno <- row$phenotype
  sex_label <- ifelse(sex == "Female", "Fem", "Male")

  rda_path <- file.path(results_dir, pheno,
    paste0("deseqAnalysisACE_", pheno, "_", ct, "_", sex_label, ".rda"))

  if (!file.exists(rda_path)) {
    cat("  SKIP (file not found):", rda_path, "\n")
    next
  }

  cat(sprintf("  %s %s %s ...\n", ct, sex, pheno))
  load(rda_path)
  df_res <- as.data.frame(res) %>% mutate(gene = rownames(res)) %>% filter(!is.na(padj))

  for (direction in c("down", "up")) {
    if (direction == "down") {
      sig_genes <- df_res %>% filter(padj < 0.05, log2FoldChange < 0) %>% pull(gene)
      dir_label <- "Downregulated"
    } else {
      sig_genes <- df_res %>% filter(padj < 0.05, log2FoldChange > 0) %>% pull(gene)
      dir_label <- "Upregulated"
    }

    if (length(sig_genes) < 5) next

    tryCatch({
      ego <- enrichGO(
        gene = sig_genes, universe = df_res$gene,
        OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
        ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05
      )

      n_go <- nrow(ego@result %>% filter(p.adjust < 0.05))
      cat(sprintf("    %s: %d sig GO terms\n", dir_label, n_go))

      if (n_go > 0) {
        # Save CSV
        tag <- paste0(ct, "_", sex_label, "_", pheno, "_", direction)
        csv_path <- file.path(output_dir, paste0("go_", tag, ".csv"))
        write.csv(ego@result, csv_path, row.names = FALSE)

        # Save dotplot
        pdf_path <- file.path(fig_dir, paste0("go_dotplot_", tag, ".pdf"), width = 10, height = 8)
        pdf(file.path(fig_dir, paste0("go_dotplot_", tag, ".pdf")), width = 10, height = 8)
        print(dotplot(ego, showCategory = 20,
          title = paste0("GO BP: ", dir_label, " in ", ct, " (", sex, ")\n", pheno)))
        dev.off()

        # Collect top terms for JSON
        top_terms <- ego@result %>%
          filter(p.adjust < 0.05) %>%
          head(10) %>%
          select(ID, Description, GeneRatio, pvalue, p.adjust, Count)

        go_results_all[[length(go_results_all) + 1]] <- list(
          cell_type = ct, sex = sex, phenotype = pheno, direction = direction,
          n_input_genes = length(sig_genes), n_go_terms = n_go,
          top_terms = top_terms
        )
      }
    }, error = function(e) {
      cat(sprintf("    WARN: GO failed for %s %s: %s\n", ct, direction, conditionMessage(e)))
    })
  }
}

# ── Write GO summary JSON ─────────────────────────────────────────────────────
json_path <- file.path(output_dir, "go_summary.json")
write_json(go_results_all, json_path, pretty = TRUE, auto_unbox = TRUE)
cat("\nWrote:", json_path, "\n")
cat("GO enrichment complete.\n")
