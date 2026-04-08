library(dplyr)
library(ggplot2)
library(patchwork)
library(readr)
library(sccomp)
library(tidyr)

args <- commandArgs(trailingOnly = TRUE)
integration <- "derived_batch"
output_root <- Sys.getenv("ACE_PROP_OUTPUT_ROOT", unset = "")

i <- 1
while (i <= length(args)) {
  if (args[i] == "--integration" && i < length(args)) {
    integration <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--output-root" && i < length(args)) {
    output_root <- args[i + 1]
    i <- i + 2
  } else {
    i <- i + 1
  }
}

if (output_root == "") {
  stop("Missing --output-root or ACE_PROP_OUTPUT_ROOT.")
}

cat("=== sccomp Visualization ===\n")
cat("Integration:", integration, "\n")
cat("Output root:", output_root, "\n\n")

data_dir <- file.path(output_root, "data")
results_base <- file.path(output_root, paste0("results_", integration))
fig_dir <- file.path(output_root, "figures", integration)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

find_rds_files <- function(base_dir) {
  rds_files <- list.files(base_dir, pattern = "\\.rds$", recursive = TRUE, full.names = TRUE)
  rds_files[grepl("^sccomp_", basename(rds_files))]
}

safe_ggsave <- function(plot, filename, width = 10, height = 8) {
  tryCatch({
    ggsave(filename, plot, width = width, height = height, dpi = 150)
    cat("  Saved:", filename, "\n")
  }, error = function(e) {
    cat("  Warning: could not save", filename, ":", conditionMessage(e), "\n")
  })
}

all_results <- tibble()

cat("--- Generating per-model plots ---\n")
rds_files <- find_rds_files(results_base)
cat("Found", length(rds_files), "model result files\n")

for (rds_path in rds_files) {
  model_name <- tools::file_path_sans_ext(basename(rds_path))
  cat("\nPlotting:", model_name, "\n")

  result <- tryCatch(readRDS(rds_path), error = function(e) {
    cat("  Could not load:", conditionMessage(e), "\n")
    NULL
  })
  if (is.null(result)) {
    next
  }

  model_fig_dir <- file.path(fig_dir, "per_model")
  dir.create(model_fig_dir, recursive = TRUE, showWarnings = FALSE)

  tryCatch({
    p <- result |> sccomp_boxplot()
    safe_ggsave(p, file.path(model_fig_dir, paste0(model_name, "_boxplot.png")), width = 14, height = 10)
  }, error = function(e) cat("  boxplot error:", conditionMessage(e), "\n"))

  tryCatch({
    p <- result |> plot_1D_intervals()
    safe_ggsave(p, file.path(model_fig_dir, paste0(model_name, "_intervals_1D.png")), width = 10, height = 8)
  }, error = function(e) cat("  1D intervals error:", conditionMessage(e), "\n"))

  tryCatch({
    p <- result |> plot_2D_intervals()
    safe_ggsave(p, file.path(model_fig_dir, paste0(model_name, "_intervals_2D.png")), width = 10, height = 8)
  }, error = function(e) cat("  2D intervals error:", conditionMessage(e), "\n"))
}

cat("\n--- Generating summary heatmap ---\n")
csv_files <- list.files(results_base, pattern = "_results\\.csv$", recursive = TRUE, full.names = TRUE)
if (length(csv_files) > 0) {
  all_results <- bind_rows(lapply(csv_files, function(path) {
    df <- read_csv(path, show_col_types = FALSE)
    parts <- strsplit(tools::file_path_sans_ext(basename(path)), "_")[[1]]
    parts <- parts[parts != "sccomp" & parts != "results"]
    if (length(parts) >= 3) {
      df$resolution <- parts[length(parts)]
      df$sex_stratum <- parts[length(parts) - 1]
      df$phenotype_var <- paste(parts[1:(length(parts) - 2)], collapse = "_")
    }
    df
  }))

  heatmap_data <- all_results |>
    filter(resolution == "fine") |>
    select(cell_group, phenotype_var, sex_stratum, c_effect, c_FDR) |>
    mutate(sig_label = case_when(
      c_FDR < 0.01 ~ "**",
      c_FDR < 0.05 ~ "*",
      c_FDR < 0.1 ~ ".",
      TRUE ~ ""
    ))

  if (nrow(heatmap_data) > 0) {
    for (sex_val in unique(heatmap_data$sex_stratum)) {
      plot_data <- heatmap_data |> filter(sex_stratum == sex_val)
      if (nrow(plot_data) == 0) {
        next
      }

      p <- ggplot(plot_data, aes(x = phenotype_var, y = cell_group, fill = c_effect)) +
        geom_tile(color = "white", linewidth = 0.5) +
        geom_text(aes(label = sig_label), size = 4, color = "black") +
        scale_fill_gradient2(
          low = "#2166AC",
          mid = "white",
          high = "#B2182B",
          midpoint = 0,
          name = "Effect\n(logit scale)"
        ) +
        labs(
          title = paste0("sccomp Composition Effects - ", integration, " integration"),
          subtitle = paste0("Sex: ", sex_val, " | * FDR<0.05  ** FDR<0.01  . FDR<0.1"),
          x = "Phenotype",
          y = "Cell Type"
        ) +
        theme_minimal(base_size = 12) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid = element_blank()
        )

      safe_ggsave(p, file.path(fig_dir, paste0("heatmap_", sex_val, ".png")), width = 12, height = 8)
    }
  }
} else {
  cat("  No results CSVs found, skipping heatmap\n")
}

cat("\n--- Generating proportion barplots ---\n")
counts_file <- file.path(data_dir, paste0("cell_counts_fine_", integration, ".csv"))
meta_file <- file.path(data_dir, paste0("metadata_", integration, ".csv"))

if (file.exists(counts_file) && file.exists(meta_file)) {
  counts <- read_csv(counts_file, show_col_types = FALSE) |>
    mutate(sample = as.character(sample))
  meta <- read_csv(meta_file, show_col_types = FALSE) |>
    mutate(sample = as.character(projid))

  sample_totals <- counts |>
    group_by(sample) |>
    summarise(total = sum(count), .groups = "drop")

  prop_data <- counts |>
    left_join(sample_totals, by = "sample") |>
    mutate(proportion = count / total) |>
    left_join(meta |> select(sample, tot_adverse_exp), by = "sample") |>
    filter(!is.na(tot_adverse_exp))

  sample_order <- prop_data |>
    distinct(sample, tot_adverse_exp) |>
    arrange(tot_adverse_exp) |>
    pull(sample)

  prop_data$sample <- factor(prop_data$sample, levels = sample_order)

  p <- ggplot(prop_data, aes(x = sample, y = proportion, fill = cell_group)) +
    geom_bar(stat = "identity", width = 1) +
    labs(
      title = paste0("Cell Type Proportions by Sample - ", integration),
      x = "Samples (ordered by increasing ACE score)",
      y = "Proportion",
      fill = "Cell Type"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank()
    )

  safe_ggsave(p, file.path(fig_dir, "proportion_barplot_by_ace.png"), width = 16, height = 8)
} else {
  cat("  Count/metadata files not found, skipping barplots\n")
}

cat("\n--- Generating sensitivity comparison ---\n")
if (nrow(all_results) > 0) {
  sensitivity_data <- all_results |>
    filter(phenotype_var == "tot_adverse_exp", sex_stratum == "all") |>
    select(cell_group, resolution, c_effect, c_lower, c_upper, c_FDR)

  if (nrow(sensitivity_data) > 0 && n_distinct(sensitivity_data$resolution) > 1) {
    p <- ggplot(sensitivity_data, aes(
      x = c_effect,
      y = reorder(cell_group, c_effect),
      xmin = c_lower,
      xmax = c_upper,
      color = resolution
    )) +
      geom_pointrange(position = position_dodge(width = 0.5), size = 0.6) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      scale_color_manual(values = c("fine" = "#E41A1C", "broad" = "#377EB8")) +
      labs(
        title = paste0("Fine vs Broad Resolution - tot_adverse_exp (", integration, ")"),
        x = "Composition Effect (logit scale)",
        y = "Cell Type",
        color = "Resolution"
      ) +
      theme_minimal(base_size = 12)

    safe_ggsave(p, file.path(fig_dir, "sensitivity_fine_vs_broad.png"), width = 10, height = 8)
  }
}

cat("\n--- Generating forest plots ---\n")
if (nrow(all_results) > 0) {
  primary_phenotypes <- c("tot_adverse_exp", "early_hh_ses", "ace_aggregate")
  forest_data <- all_results |>
    filter(phenotype_var %in% primary_phenotypes, resolution == "fine", sex_stratum == "all")

  if (nrow(forest_data) > 0) {
    p <- ggplot(forest_data, aes(
      x = c_effect,
      y = reorder(cell_group, c_effect),
      xmin = c_lower,
      xmax = c_upper,
      color = c_FDR < 0.05
    )) +
      geom_pointrange(size = 0.5) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      scale_color_manual(
        values = c("FALSE" = "grey60", "TRUE" = "#E41A1C"),
        labels = c("NS", "FDR < 0.05"),
        name = ""
      ) +
      facet_wrap(~ phenotype_var, scales = "free_x") +
      labs(
        title = paste0("Cell Type Composition Effects - ", integration),
        x = "Effect Size (logit scale)",
        y = "Cell Type"
      ) +
      theme_minimal(base_size = 11) +
      theme(strip.text = element_text(face = "bold"))

    safe_ggsave(p, file.path(fig_dir, "forest_primary_phenotypes.png"), width = 14, height = 8)
  }
}

cat("\n=== Visualization complete ===\n")
cat("Figures saved to:", fig_dir, "\n")
