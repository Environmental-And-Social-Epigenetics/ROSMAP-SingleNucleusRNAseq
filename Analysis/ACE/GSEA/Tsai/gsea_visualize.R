#!/usr/bin/env Rscript
# gsea_visualize.R -- Visualization of GSEA enrichment results
#
# Usage:
#   Rscript gsea_visualize.R \
#     --results-dir /path/to/gsea/results \
#     --output-dir  /path/to/figures \
#     --phenotype   tot_adverse_exp \
#     --sex         Female
#
# Reads gsea_summary.csv produced by gsea_analysis.R and generates:
#   1. Bubble heatmap of GO-BP terms across cell types
#   2. Per-cell-type barplots of top enriched terms
#   3. Pathway convergence table (terms shared across cell types)

library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# ── Parse arguments ──────────────────────────────────────────────────────────

args <- commandArgs(trailingOnly = TRUE)
results_dir <- NULL
output_dir  <- NULL
phenotype   <- NULL
sex         <- NULL

i <- 1
while (i <= length(args)) {
  if (args[i] == "--results-dir" && i < length(args)) {
    results_dir <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--output-dir" && i < length(args)) {
    output_dir <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--phenotype" && i < length(args)) {
    phenotype <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--sex" && i < length(args)) {
    sex <- args[i + 1]; i <- i + 2
  } else {
    i <- i + 1
  }
}

if (is.null(results_dir)) stop("ERROR: --results-dir is required.")
if (is.null(output_dir))  stop("ERROR: --output-dir is required.")
if (is.null(phenotype))   stop("ERROR: --phenotype is required.")
if (is.null(sex))         stop("ERROR: --sex is required.")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== ACE GSEA Visualization ===\n")
cat("Results dir:", results_dir, "\n")
cat("Output dir: ", output_dir, "\n")
cat("Phenotype:  ", phenotype, "\n")
cat("Sex:        ", sex, "\n\n")

# ── Load summary ─────────────────────────────────────────────────────────────

summary_path <- file.path(results_dir, "gsea_summary.csv")
if (!file.exists(summary_path)) {
  stop("ERROR: gsea_summary.csv not found in ", results_dir)
}

df <- read.csv(summary_path, stringsAsFactors = FALSE)

if (nrow(df) == 0) {
  cat("Summary file is empty -- no enrichment results to visualize.\n")
  quit(status = 0)
}

cat("Loaded", nrow(df), "enriched gene sets across",
    length(unique(df$cell_type)), "cell types and",
    length(unique(df$database)), "databases.\n\n")

# ── Helper: truncate long labels ─────────────────────────────────────────────

truncate_label <- function(x, max_len = 60) {
  ifelse(nchar(x) > max_len, paste0(substr(x, 1, max_len - 3), "..."), x)
}

# ── Figure 1: Bubble heatmap (GO-BP terms shared across cell types) ──────────

cat("── Figure 1: Bubble heatmap ──\n")

bp_df <- df %>%
  filter(database == "Biological_Process_noRedundant")

if (nrow(bp_df) > 0) {
  # Filter to terms significant (FDR < 0.05) in at least 2 cell types
  shared_terms <- bp_df %>%
    filter(FDR < 0.05) %>%
    group_by(description) %>%
    summarise(n_celltypes = n_distinct(cell_type), .groups = "drop") %>%
    filter(n_celltypes >= 2) %>%
    pull(description)

  if (length(shared_terms) == 0) {
    # Relax: FDR < 0.2 in at least 2 cell types
    shared_terms <- bp_df %>%
      filter(FDR < 0.2) %>%
      group_by(description) %>%
      summarise(n_celltypes = n_distinct(cell_type), .groups = "drop") %>%
      filter(n_celltypes >= 2) %>%
      arrange(desc(n_celltypes)) %>%
      slice_head(n = 30) %>%
      pull(description)
    cat("  Relaxed threshold to FDR < 0.2 for shared terms.\n")
  }

  if (length(shared_terms) > 0) {
    # Limit to top 40 terms for readability
    if (length(shared_terms) > 40) {
      # Prioritize by mean absolute NES
      term_importance <- bp_df %>%
        filter(description %in% shared_terms) %>%
        group_by(description) %>%
        summarise(mean_abs_nes = mean(abs(NES), na.rm = TRUE), .groups = "drop") %>%
        arrange(desc(mean_abs_nes)) %>%
        slice_head(n = 40)
      shared_terms <- term_importance$description
    }

    plot_df <- bp_df %>%
      filter(description %in% shared_terms) %>%
      mutate(description = truncate_label(description))

    # Order terms by hierarchical clustering on NES matrix
    nes_mat <- plot_df %>%
      select(description, cell_type, NES) %>%
      pivot_wider(names_from = cell_type, values_from = NES, values_fill = 0) %>%
      tibble::column_to_rownames("description")

    if (nrow(nes_mat) > 2 && ncol(nes_mat) > 1) {
      term_order <- rownames(nes_mat)[hclust(dist(nes_mat))$order]
      plot_df$description <- factor(plot_df$description, levels = rev(term_order))
    }

    p1 <- ggplot(plot_df, aes(x = cell_type, y = description)) +
      geom_point(aes(size = leadingEdgeNum, color = NES)) +
      scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                            name = "NES") +
      scale_size_continuous(name = "Leading\nEdge Genes", range = c(2, 8)) +
      theme_bw(base_size = 11) +
      theme(
        axis.text.y = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        plot.title  = element_text(size = 12, face = "bold"),
        legend.position = "right"
      ) +
      labs(
        x     = "Cell Type",
        y     = "GO Biological Process",
        title = paste0("GSEA: GO-BP Terms Shared Across Cell Types\n",
                       phenotype, " / ", sex)
      )

    fig1_path <- file.path(output_dir,
                           paste0("bubble_heatmap_", phenotype, "_", sex, ".pdf"))
    n_terms <- length(unique(plot_df$description))
    fig_height <- max(8, 0.3 * n_terms + 2)
    ggsave(fig1_path, plot = p1, width = 10, height = fig_height, limitsize = FALSE)
    cat("  Saved:", basename(fig1_path), "\n")
  } else {
    cat("  No GO-BP terms shared across 2+ cell types. Skipping bubble heatmap.\n")
  }
} else {
  cat("  No GO-BP results found. Skipping bubble heatmap.\n")
}

# ── Figure 2: Per-cell-type barplots (top 15 terms by FDR) ──────────────────

cat("── Figure 2: Per-cell-type barplots ──\n")

cell_types <- unique(df$cell_type)

for (ct in cell_types) {
  ct_df <- df %>%
    filter(cell_type == ct, FDR < 0.2) %>%
    arrange(FDR) %>%
    slice_head(n = 15) %>%
    mutate(
      description = truncate_label(description, 50),
      direction   = ifelse(NES > 0, "Upregulated", "Downregulated"),
      description = factor(description, levels = rev(description))
    )

  if (nrow(ct_df) == 0) {
    cat("  ", ct, ": no significant terms (FDR < 0.2). Skipping.\n")
    next
  }

  p2 <- ggplot(ct_df, aes(x = NES, y = description, fill = direction)) +
    geom_col() +
    scale_fill_manual(values = c("Upregulated" = "firebrick3",
                                 "Downregulated" = "steelblue3"),
                      name = "NES Direction") +
    theme_bw(base_size = 11) +
    theme(
      axis.text.y = element_text(size = 9),
      plot.title  = element_text(size = 12, face = "bold")
    ) +
    labs(
      x     = "Normalized Enrichment Score (NES)",
      y     = NULL,
      title = paste0(ct, " -- Top Enriched Terms\n", phenotype, " / ", sex)
    )

  fig2_path <- file.path(output_dir,
                         paste0("barplot_", ct, "_", phenotype, "_", sex, ".pdf"))
  fig_height <- max(5, 0.35 * nrow(ct_df) + 2)
  ggsave(fig2_path, plot = p2, width = 9, height = fig_height)
  cat("  Saved:", basename(fig2_path), "\n")
}

# ── Figure 3: Pathway convergence table ─────────────────────────────────────

cat("── Figure 3: Pathway convergence table ──\n")

convergence <- bp_df %>%
  filter(FDR < 0.05) %>%
  group_by(geneSet, description) %>%
  summarise(
    n_cell_types     = n_distinct(cell_type),
    cell_types       = paste(sort(unique(cell_type)), collapse = ", "),
    mean_NES         = round(mean(NES, na.rm = TRUE), 3),
    min_FDR          = signif(min(FDR, na.rm = TRUE), 3),
    .groups          = "drop"
  ) %>%
  arrange(desc(n_cell_types), min_FDR)

convergence_path <- file.path(output_dir,
                              paste0("convergence_table_", phenotype, "_", sex, ".csv"))
write.csv(convergence, convergence_path, row.names = FALSE)
cat("  Saved:", basename(convergence_path), "(", nrow(convergence), "terms )\n")

# Also create a PDF table for the top terms
if (nrow(convergence) > 0) {
  top_conv <- convergence %>%
    filter(n_cell_types >= 2) %>%
    slice_head(n = 30) %>%
    mutate(description = truncate_label(description, 55))

  if (nrow(top_conv) > 0) {
    p3 <- ggplot(top_conv, aes(x = n_cell_types, y = reorder(description, n_cell_types))) +
      geom_col(aes(fill = mean_NES)) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                           name = "Mean NES") +
      theme_bw(base_size = 11) +
      theme(
        axis.text.y = element_text(size = 9),
        plot.title  = element_text(size = 12, face = "bold")
      ) +
      labs(
        x     = "Number of Cell Types (FDR < 0.05)",
        y     = NULL,
        title = paste0("GO-BP Pathway Convergence\n", phenotype, " / ", sex)
      )

    fig3_path <- file.path(output_dir,
                           paste0("convergence_plot_", phenotype, "_", sex, ".pdf"))
    fig_height <- max(5, 0.35 * nrow(top_conv) + 2)
    ggsave(fig3_path, plot = p3, width = 9, height = fig_height)
    cat("  Saved:", basename(fig3_path), "\n")
  }
}

cat("\n=== Visualization complete ===\n")
