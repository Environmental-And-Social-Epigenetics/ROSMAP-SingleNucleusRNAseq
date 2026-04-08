library(sccomp)
library(dplyr)
library(readr)
library(tidyr)

ace_components <- c(
  "emotional_neglect",
  "family_pro_sep",
  "financial_need",
  "parental_intimidation",
  "parental_violence"
)

args <- commandArgs(trailingOnly = TRUE)
integration <- "library_id"
sex_stratum <- "all"
resolution <- "fine"
output_root <- Sys.getenv("ACE_PROP_OUTPUT_ROOT", unset = "")
smoke_mode <- FALSE

i <- 1
while (i <= length(args)) {
  if (args[i] == "--integration" && i < length(args)) {
    integration <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--sex" && i < length(args)) {
    sex_stratum <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--resolution" && i < length(args)) {
    resolution <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--output-root" && i < length(args)) {
    output_root <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--smoke") {
    smoke_mode <- TRUE
    i <- i + 1
  } else {
    i <- i + 1
  }
}

if (output_root == "") {
  stop("Missing --output-root or ACE_PROP_OUTPUT_ROOT.")
}

cat("=== sccomp Cell Type Proportion Analysis ===\n")
cat("Integration:", integration, "\n")
cat("Sex stratum:", sex_stratum, "\n")
cat("Resolution: ", resolution, "\n")
cat("Output root:", output_root, "\n\n")

data_dir <- file.path(output_root, "data")
results_base <- file.path(output_root, paste0("results_", integration))
dir.create(results_base, recursive = TRUE, showWarnings = FALSE)

counts_file <- file.path(data_dir, paste0("cell_counts_", resolution, "_", integration, ".csv"))
meta_file <- file.path(data_dir, paste0("metadata_", integration, ".csv"))

if (!file.exists(counts_file)) stop("Count file not found: ", counts_file)
if (!file.exists(meta_file)) stop("Metadata file not found: ", meta_file)

counts <- read_csv(counts_file, show_col_types = FALSE)
meta <- read_csv(meta_file, show_col_types = FALSE) |>
  mutate(sample = as.character(projid))

required_meta <- c(
  "sample",
  "tot_adverse_exp",
  "early_hh_ses",
  "ace_aggregate",
  "msex",
  "age_death",
  "pmi",
  "niareagansc",
  ace_components
)
missing_meta <- setdiff(required_meta, names(meta))
if (length(missing_meta) > 0) {
  stop("Metadata file is missing required columns: ", paste(missing_meta, collapse = ", "))
}

dat <- counts |>
  mutate(sample = as.character(sample)) |>
  left_join(meta, by = "sample")

dat_ace <- dat |> filter(!is.na(tot_adverse_exp))
n_samples <- n_distinct(dat_ace$sample)
cat("Samples with ACE data:", n_samples, "\n")

if (sex_stratum == "male") {
  dat_ace <- dat_ace |> filter(msex == 1)
  cat("Filtered to males:", n_distinct(dat_ace$sample), "samples\n")
} else if (sex_stratum == "female") {
  dat_ace <- dat_ace |> filter(msex == 0)
  cat("Filtered to females:", n_distinct(dat_ace$sample), "samples\n")
}

if (smoke_mode) {
  if (n_distinct(dat_ace$sample) == 0) {
    stop("Smoke mode found no samples with ACE data.")
  }

  smoke_dir <- file.path(results_base, "smoke")
  dir.create(smoke_dir, recursive = TRUE, showWarnings = FALSE)
  smoke_summary <- tibble(
    integration = integration,
    sex = sex_stratum,
    resolution = resolution,
    samples = n_distinct(dat_ace$sample),
    cell_groups = n_distinct(dat_ace$cell_group),
    total_count = sum(dat_ace$count, na.rm = TRUE)
  )
  write_csv(smoke_summary, file.path(smoke_dir, "sccomp_smoke_summary.csv"))
  cat("Smoke summary written\n")
  quit(save = "no", status = 0)
}

if (n_distinct(dat_ace$sample) < 20) {
  cat("SKIP: fewer than 20 samples (", n_distinct(dat_ace$sample), ") -- insufficient for sccomp\n")
  quit(save = "no", status = 0)
}

scale_col <- function(df, col) {
  values <- as.numeric(df[[col]])
  if (sum(!is.na(values)) > 1) {
    df[[col]] <- as.numeric(scale(values))
  }
  df
}

dat_ace <- dat_ace |>
  scale_col("age_death") |>
  scale_col("pmi")

covars <- if (sex_stratum == "all") {
  "age_death + msex + pmi + niareagansc"
} else {
  "age_death + pmi + niareagansc"
}

models <- list(
  list(name = "primary", phenotype = "tot_adverse_exp", subdir = "primary", diff_var = FALSE),
  list(name = "ses", phenotype = "early_hh_ses", subdir = "ses", diff_var = FALSE),
  list(name = "aggregate", phenotype = "ace_aggregate", subdir = "aggregate", diff_var = FALSE),
  list(name = "variability", phenotype = "tot_adverse_exp", subdir = "variability", diff_var = TRUE)
)

for (component in ace_components) {
  models[[length(models) + 1]] <- list(
    name = component,
    phenotype = component,
    subdir = "components",
    diff_var = FALSE
  )
}

run_sccomp_model <- function(dat, phenotype, covars, diff_var, out_dir, label) {
  cat("\n--- Model:", label, "---\n")

  dat_sub <- dat |> filter(!is.na(.data[[phenotype]]))
  n <- n_distinct(dat_sub$sample)
  cat("  Samples with non-NA", phenotype, ":", n, "\n")
  if (n < 20) {
    cat("  SKIP: too few samples\n")
    return(NULL)
  }

  dat_sub <- dat_sub |> scale_col(phenotype)

  for (column in c("age_death", "pmi", "niareagansc", "msex", phenotype)) {
    if (column %in% names(dat_sub)) {
      dat_sub[[column]] <- as.numeric(dat_sub[[column]])
    }
  }

  formula_comp <- as.formula(paste0("~ ", phenotype, " + ", covars))
  formula_var <- if (diff_var) as.formula(paste0("~ ", phenotype)) else ~ 1

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  file_prefix <- paste0("sccomp_", phenotype, "_", sex_stratum, "_", resolution)

  result <- tryCatch({
    dat_sub |>
      sccomp_estimate(
        formula_composition = formula_comp,
        formula_variability = formula_var,
        .sample = sample,
        .cell_group = cell_group,
        .count = count,
        bimodal_mean_variability_association = TRUE,
        cores = parallel::detectCores(),
        mcmc_seed = 42,
        verbose = FALSE
      ) |>
      sccomp_remove_outliers(
        cores = parallel::detectCores(),
        verbose = FALSE
      ) |>
      sccomp_test()
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    writeLines(conditionMessage(e), file.path(out_dir, paste0(file_prefix, "_ERROR.txt")))
    return(NULL)
  })

  if (is.null(result)) {
    return(NULL)
  }

  saveRDS(result, file.path(out_dir, paste0(file_prefix, ".rds")))

  res_tbl <- tryCatch({
    result |>
      as_tibble() |>
      filter(parameter == phenotype)
  }, error = function(e) {
    cat("  Warning: could not extract results table:", conditionMessage(e), "\n")
    NULL
  })

  if (!is.null(res_tbl) && nrow(res_tbl) > 0) {
    write_csv(res_tbl, file.path(out_dir, paste0(file_prefix, "_results.csv")))
  }

  fc_tbl <- tryCatch({
    result |>
      sccomp_proportional_fold_change(
        formula_composition = formula_comp,
        from = -1,
        to = 1
      )
  }, error = function(e) {
    cat("  Warning: could not compute fold change:", conditionMessage(e), "\n")
    NULL
  })

  if (!is.null(fc_tbl)) {
    write_csv(fc_tbl, file.path(out_dir, paste0(file_prefix, "_fc.csv")))
  }

  outliers <- tryCatch(attr(result, "outliers"), error = function(e) NULL)
  if (!is.null(outliers) && nrow(outliers) > 0) {
    write_csv(outliers, file.path(out_dir, paste0(file_prefix, "_outliers.csv")))
  }

  invisible(result)
}

for (model in models) {
  model_dir <- file.path(results_base, model$subdir)
  actual_out <- if (resolution == "broad") {
    file.path(results_base, "sensitivity")
  } else {
    model_dir
  }

  run_sccomp_model(
    dat = dat_ace,
    phenotype = model$phenotype,
    covars = covars,
    diff_var = model$diff_var,
    out_dir = actual_out,
    label = paste(model$name, sex_stratum, resolution, sep = " / ")
  )
}

cat("\n=== sccomp analysis complete ===\n")
cat("Results saved under:", results_base, "\n")
