#!/usr/bin/env Rscript
# =============================================================================
# ACE CellChat Visualization
#
# Generates publication-quality figures from CellChat differential analysis.
#
# Usage:
#   Rscript cellchat_visualize.R \
#     --results-dir /path/to/results \
#     --output-dir /path/to/figures
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

results_dir <- NULL
output_dir  <- NULL

i <- 1
while (i <= length(args)) {
    if (args[i] == "--results-dir") {
        i <- i + 1; results_dir <- args[i]
    } else if (args[i] == "--output-dir") {
        i <- i + 1; output_dir <- args[i]
    }
    i <- i + 1
}

stopifnot(!is.null(results_dir), !is.null(output_dir))
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
    library(CellChat)
    library(ggplot2)
    library(patchwork)
    library(dplyr)
})

# --- Load data ---------------------------------------------------------------
cat("Loading results...\n")
diff_df <- read.csv(file.path(results_dir, "differential_interactions.csv"))
pathway_df <- read.csv(file.path(results_dir, "pathway_changes.csv"))
focus_df <- read.csv(file.path(results_dir, "focus_axes_results.csv"))

# --- 1. Differential interaction heatmap -------------------------------------
cat("Generating differential interaction heatmap...\n")
if (nrow(diff_df) > 0) {
    # Pivot to matrix
    cell_types <- sort(unique(c(diff_df$sender, diff_df$receiver)))
    mat <- matrix(0, nrow = length(cell_types), ncol = length(cell_types),
                  dimnames = list(cell_types, cell_types))
    for (r in seq_len(nrow(diff_df))) {
        mat[diff_df$sender[r], diff_df$receiver[r]] <- diff_df$count_diff[r]
    }

    pdf(file.path(output_dir, "heatmap_differential_interactions.pdf"),
        width = 10, height = 8)
    heatmap(mat, scale = "none", main = "Differential Interactions (High - Low)",
            col = colorRampPalette(c("blue", "white", "red"))(100))
    dev.off()
}

# --- 2. Pathway changes barplot ----------------------------------------------
cat("Generating pathway changes barplot...\n")
if (nrow(pathway_df) > 0) {
    top_paths <- pathway_df %>%
        arrange(desc(abs(prob_diff))) %>%
        head(30)

    p <- ggplot(top_paths, aes(x = reorder(pathway, prob_diff), y = prob_diff,
                                fill = prob_diff > 0)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        scale_fill_manual(values = c("TRUE" = "#D73027", "FALSE" = "#4575B4"),
                          labels = c("TRUE" = "Higher in ACE-high",
                                     "FALSE" = "Higher in ACE-low"),
                          name = "Direction") +
        labs(x = "Signaling Pathway", y = "Communication Probability Difference",
             title = "Top Altered Signaling Pathways (ACE-high vs ACE-low)") +
        theme_bw(base_size = 12) +
        theme(legend.position = "bottom")

    ggsave(file.path(output_dir, "barplot_pathway_changes.pdf"), p,
           width = 10, height = 8)
}

# --- 3. Focus axes barplot ---------------------------------------------------
cat("Generating focus axes barplot...\n")
if (nrow(focus_df) > 0) {
    p <- ggplot(focus_df, aes(x = reorder(axis, n_diff), y = n_diff,
                               fill = n_diff > 0)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        scale_fill_manual(values = c("TRUE" = "#D73027", "FALSE" = "#4575B4"),
                          labels = c("TRUE" = "More in ACE-high",
                                     "FALSE" = "More in ACE-low"),
                          name = "Direction") +
        labs(x = "Interaction Axis", y = "Interaction Count Difference",
             title = "Focused Cell-Cell Communication Axes") +
        theme_bw(base_size = 14)

    ggsave(file.path(output_dir, "barplot_focus_axes.pdf"), p,
           width = 8, height = 5)
}

# --- 4. Load merged object for chord diagrams if available -------------------
merged_file <- file.path(results_dir, "cellchat_merged.rds")
if (file.exists(merged_file)) {
    cat("Loading merged CellChat object for chord diagrams...\n")
    cellchat_merged <- readRDS(merged_file)

    tryCatch({
        pdf(file.path(output_dir, "chord_interaction_comparison.pdf"),
            width = 12, height = 6)
        par(mfrow = c(1, 2))
        netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE,
                                   measure = "count")
        netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE,
                                   measure = "weight")
        dev.off()
    }, error = function(e) {
        cat("  Chord diagram failed:", conditionMessage(e), "\n")
    })
}

cat("\n=== Visualization complete ===\n")
cat("Figures saved to:", output_dir, "\n")
