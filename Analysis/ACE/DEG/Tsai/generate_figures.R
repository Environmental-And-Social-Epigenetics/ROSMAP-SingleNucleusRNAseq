#!/usr/bin/env Rscript
# generate_figures.R -- Generate publication-quality figures from DEG results
#
# Usage:
#   Rscript generate_figures.R \
#     --results-dir results_batch \
#     --summary-csv results_batch/deg_summary.csv \
#     --output-dir results_batch/figures \
#     --integration-label "Derived-batch"

library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)
library(tibble)

# ── Parse arguments ──────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
results_dir <- NULL
summary_csv <- NULL
output_dir <- NULL
integration_label <- NULL

i <- 1
while (i <= length(args)) {
  if (args[i] == "--results-dir" && i < length(args)) {
    results_dir <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--summary-csv" && i < length(args)) {
    summary_csv <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--output-dir" && i < length(args)) {
    output_dir <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--integration-label" && i < length(args)) {
    integration_label <- args[i + 1]; i <- i + 2
  } else {
    i <- i + 1
  }
}

# Defaults for backward compatibility (run from the Tsai directory)
if (is.null(results_dir)) results_dir <- "results_projid"
if (is.null(summary_csv)) summary_csv <- file.path(results_dir, "deg_summary.csv")
if (is.null(output_dir)) output_dir <- file.path(results_dir, "figures")
if (is.null(integration_label)) {
  integration_label <- sub("^results_", "", basename(results_dir))
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Results dir:", results_dir, "\n")
cat("Summary CSV:", summary_csv, "\n")
cat("Output dir:", output_dir, "\n")
cat("Integration:", integration_label, "\n\n")

# ── Load summary data ────────────────────────────────────────────────────────
df <- read.csv(summary_csv, stringsAsFactors = FALSE)

# Handle both old format (ct_sex column) and new format (cell_type + sex columns)
if (!"cell_type" %in% names(df) && "ct_sex" %in% names(df)) {
  df$sex <- ifelse(grepl("_Fem$", df$ct_sex), "Female", "Male")
  df$cell_type <- gsub("^(adverse_exp_|hh_ses_|aggregate_)", "", df$ct_sex)
  df$cell_type <- gsub("_(Fem|Male)$", "", df$cell_type)
}

ct_order <- c("Exc", "Inh", "Ast", "Mic", "Oli", "OPC", "Endo",
              "Ex-L2_3", "Ex-L4", "Ex-L4_5", "Ex-L5", "Ex-L5_6", "Ex-L5_6-CC", "Ex-NRGN",
              "In-PV_Basket", "In-PV_Chandelier", "In-Rosehip", "In-SST", "In-VIP")

pheno_labels <- c(
  "tot_adverse_exp" = "Total Adverse\nExperiences",
  "early_hh_ses"    = "Early Household\nSES",
  "ace_aggregate"   = "ACE Aggregate\nScore"
)

df$cell_type <- factor(df$cell_type, levels = ct_order)
df$phenotype_label <- pheno_labels[df$phenotype]

# ── FIGURE 1: Heatmaps per phenotype ─────────────────────────────────────────
for (pheno in unique(df$phenotype)) {
  sub <- df %>% filter(phenotype == pheno)
  mat <- sub %>%
    select(cell_type, sex, n_sig) %>%
    pivot_wider(names_from = sex, values_from = n_sig, values_fill = list(n_sig = 0L)) %>%
    filter(!is.na(cell_type)) %>%
    arrange(cell_type) %>%
    column_to_rownames("cell_type") %>%
    as.matrix()

  pdf(file.path(output_dir, paste0("heatmap_", pheno, ".pdf")), width = 5, height = 8)
  pheatmap(log10(mat + 1),
    cluster_rows = FALSE, cluster_cols = FALSE,
    color = colorRampPalette(c("white", "#FEE08B", "#FC8D59", "#D73027"))(100),
    display_numbers = mat, number_format = "%.0f",
    main = paste0("DEGs (padj<0.05) - ", gsub("_", " ", pheno), "\n", integration_label),
    fontsize = 10, fontsize_number = 8, angle_col = 0, border_color = "grey80")
  dev.off()
  cat("Heatmap done:", pheno, "\n")
}

# ── FIGURE 2: Bar plot by phenotype x sex ─────────────────────────────────────
totals <- df %>%
  group_by(phenotype_label, sex) %>%
  summarize(total_degs = sum(n_sig), .groups = "drop")

pdf(file.path(output_dir, "barplot_phenotype_sex.pdf"), width = 8, height = 5)
print(ggplot(totals, aes(x = phenotype_label, y = total_degs, fill = sex)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("Female" = "#E78AC3", "Male" = "#66C2A5")) +
  labs(x = "", y = "Total Significant DEGs (padj < 0.05)",
       title = "Differential Expression by ACE Phenotype and Sex",
       subtitle = paste(integration_label, "integration")) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top", plot.title = element_text(face = "bold")))
dev.off()
cat("Barplot phenotype done\n")

# ── FIGURE 3: Up/Down bar for tot_adverse_exp ─────────────────────────────────
updown <- df %>%
  filter(phenotype == "tot_adverse_exp") %>%
  select(cell_type, sex, n_up, n_down) %>%
  pivot_longer(cols = c(n_up, n_down), names_to = "dir", values_to = "count") %>%
  mutate(
    dir = ifelse(dir == "n_up", "Upregulated", "Downregulated"),
    count_signed = ifelse(dir == "Downregulated", -count, count)
  )

pdf(file.path(output_dir, "barplot_updown_tot_adverse.pdf"), width = 12, height = 6)
print(ggplot(updown, aes(x = cell_type, y = count_signed, fill = dir)) +
  geom_bar(stat = "identity") +
  facet_wrap(~sex, ncol = 1) +
  scale_fill_manual(values = c("Upregulated" = "#D73027", "Downregulated" = "#4575B4")) +
  coord_flip() +
  labs(x = "", y = "Number of DEGs", fill = "",
       title = "Up/Down-regulated DEGs for Total Adverse Experiences",
       subtitle = paste(integration_label, "| padj < 0.05")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top", plot.title = element_text(face = "bold")))
dev.off()
cat("Up/Down barplot done\n")

# ── FIGURE 4: Volcano plots for top 4 male cell types ────────────────────────
top4 <- df %>%
  filter(phenotype == "tot_adverse_exp", sex == "Male") %>%
  arrange(desc(n_sig)) %>%
  head(4) %>%
  pull(cell_type) %>%
  as.character()

volcano_list <- list()
for (ct in top4) {
  fname <- file.path(results_dir, "tot_adverse_exp",
    paste0("deseqAnalysisACE_tot_adverse_exp_", ct, "_Male.rda"))
  if (file.exists(fname)) {
    load(fname)
    vdf <- as.data.frame(res) %>%
      mutate(
        gene = rownames(res),
        sig = ifelse(padj < 0.05 & log2FoldChange > 0, "Up",
              ifelse(padj < 0.05 & log2FoldChange < 0, "Down", "NS")),
        cell_type = ct
      )
    volcano_list[[length(volcano_list) + 1]] <- vdf
  }
}

if (length(volcano_list) > 0) {
  volc_df <- bind_rows(volcano_list)
  volc_df$sig <- factor(volc_df$sig, levels = c("Down", "NS", "Up"))
  volc_df$cell_type <- factor(volc_df$cell_type, levels = top4)

  pdf(file.path(output_dir, "volcano_top4_male.pdf"), width = 12, height = 10)
  print(ggplot(volc_df, aes(x = log2FoldChange, y = -log10(pvalue), color = sig)) +
    geom_point(alpha = 0.4, size = 0.8) +
    facet_wrap(~cell_type, ncol = 2, scales = "free_y") +
    scale_color_manual(values = c("Down" = "#4575B4", "NS" = "grey70", "Up" = "#D73027")) +
    labs(x = "log2 Fold Change", y = "-log10(p-value)", color = "",
         title = "Volcano Plots: Males, Total Adverse Experiences",
         subtitle = integration_label) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top", plot.title = element_text(face = "bold"),
          strip.text = element_text(face = "bold", size = 12)))
  dev.off()
  cat("Volcano top4 done\n")
}

# ── FIGURE 5: Cell type vulnerability dot plot ────────────────────────────────
ct_stats <- list()
for (ct in ct_order) {
  fname <- file.path(results_dir, "tot_adverse_exp",
    paste0("deseqAnalysisACE_tot_adverse_exp_", ct, "_Male.rda"))
  if (file.exists(fname)) {
    load(fname)
    min_p <- min(res$padj, na.rm = TRUE)
    ns <- sum(res$padj < 0.05, na.rm = TRUE)
    med_lfc <- median(abs(res$log2FoldChange[res$padj < 0.05]), na.rm = TRUE)
    ct_stats[[length(ct_stats) + 1]] <- data.frame(
      cell_type = ct, n_sig = ns, min_padj = min_p,
      med_abs_lfc = ifelse(is.na(med_lfc), 0, med_lfc))
  }
}

if (length(ct_stats) > 0) {
  ct_df <- bind_rows(ct_stats)
  ct_df$cell_type <- factor(ct_df$cell_type, levels = rev(ct_order))

  pdf(file.path(output_dir, "dotplot_vulnerability.pdf"), width = 8, height = 7)
  print(ggplot(ct_df, aes(x = med_abs_lfc, y = cell_type, size = n_sig, color = -log10(min_padj))) +
    geom_point() +
    scale_size_continuous(range = c(1, 12), name = "N DEGs") +
    scale_color_gradient(low = "#FEE08B", high = "#D73027", name = "-log10(min padj)") +
    labs(x = "Median |log2FC| of Significant DEGs", y = "",
         title = "Cell Type Vulnerability to Adverse Childhood Experiences",
         subtitle = paste("Males | tot_adverse_exp |", integration_label)) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold")))
  dev.off()
  cat("Dotplot vulnerability done\n")
}

cat("\nAll figures saved to:", output_dir, "\n")
cat(paste(list.files(output_dir, pattern = "\\.pdf$"), collapse = "\n"), "\n")
