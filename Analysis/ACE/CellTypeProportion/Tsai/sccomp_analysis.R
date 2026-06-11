# Canonical, cohort-aware ACE cell-type proportion (sccomp) analysis.
#
# Single source of truth for BOTH cohorts. The DeJager sccomp_analysis.R is a
# thin wrapper that source()s this file (mirroring the DEG parity pattern where
# aceDegDJ.Rscript source()s aceDegT.Rscript). Cohort-specific bits stay
# parameterized:
#   * --cohort {tsai,dejager} (or ACE_PROP_COHORT) selects the default
#     integration label (tsai -> derived_batch, dejager -> library_id);
#   * --output-root (or ACE_PROP_OUTPUT_ROOT) selects the cohort output tree.
# The cell-grouping id column is already normalized to `sample`/`projid` by each
# cohort's prep_counts.py (Tsai keys on projid/sample_id, DeJager on patient_id,
# both written out as the metadata `projid` column and the counts `sample`
# column), so the statistical core below reads identical columns for both.
#
# Statistical core (identical across cohorts): the modern, maintained sccomp
# pipeline
#   sccomp_estimate() |> sccomp_remove_outliers() |> sccomp_test()
# plus sccomp_proportional_fold_change(). The deprecated monolithic sccomp_glm()
# is NOT used.

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
cohort <- tolower(Sys.getenv("ACE_PROP_COHORT", unset = "tsai"))
integration <- NULL
sex_stratum <- "all"
resolution <- "fine"
output_root <- Sys.getenv("ACE_PROP_OUTPUT_ROOT", unset = "")
smoke_mode <- FALSE
arm <- NULL   # male AD-model arm (e.g. MaleContAD): males-only, single primary
              # tot_adverse_exp composition model with that arm's AD covariates.

i <- 1
while (i <= length(args)) {
  if (args[i] == "--cohort" && i < length(args)) {
    cohort <- tolower(args[i + 1])
    i <- i + 2
  } else if (args[i] == "--integration" && i < length(args)) {
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
  } else if (args[i] == "--arm" && i < length(args)) {
    arm <- args[i + 1]
    i <- i + 2
  } else if (args[i] == "--smoke") {
    smoke_mode <- TRUE
    i <- i + 1
  } else {
    i <- i + 1
  }
}

# Cohort-specific default integration label (overridable via --integration).
if (is.null(integration)) {
  integration <- switch(
    cohort,
    tsai = "derived_batch",
    dejager = "library_id",
    stop("Unknown --cohort '", cohort, "'. Known: tsai, dejager.")
  )
}

# A male AD-model arm implies males-only and pulls its covariate set from the
# shared single-source-of-truth spec.
arm_spec_obj <- NULL
if (!is.null(arm)) {
  this_dir <- tryCatch(dirname(sub("^--file=", "",
    grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), error = function(e) ".")
  shared_R <- file.path(this_dir, "..", "..", "_shared", "arm_covariates.R")
  if (!file.exists(shared_R)) {
    # Portable fallback: resolve from REPO_ROOT (exported by config/paths.sh).
    repo_root <- Sys.getenv("REPO_ROOT", "")
    if (nzchar(repo_root)) {
      shared_R <- file.path(repo_root, "Analysis", "ACE", "_shared", "arm_covariates.R")
    }
  }
  source(shared_R)
  arm_spec_obj <- arm_spec(arm)
  sex_stratum <- "male"   # male AD-model arms are males-only
}

if (output_root == "") {
  stop("Missing --output-root or ACE_PROP_OUTPUT_ROOT.")
}

cat("=== sccomp Cell Type Proportion Analysis ===\n")
cat("Cohort:     ", cohort, "\n")
cat("Integration:", integration, "\n")
cat("Sex stratum:", sex_stratum, "\n")
cat("Resolution: ", resolution, "\n")
cat("Output root:", output_root, "\n\n")

data_dir <- file.path(output_root, "data")
# Arm runs go to a dedicated results dir so they don't collide with the baseline.
results_base <- if (!is.null(arm)) {
  file.path(output_root, paste0("results_", integration, "_", arm))
} else {
  file.path(output_root, paste0("results_", integration))
}
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
  mutate(sample = as.character(sample), count = as.integer(count)) |>
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

n_samples <- n_distinct(dat_ace$sample)
if (smoke_mode) {
  if (n_samples == 0) {
    stop("Smoke mode found no samples with ACE data.")
  }

  smoke_dir <- file.path(results_base, "smoke")
  dir.create(smoke_dir, recursive = TRUE, showWarnings = FALSE)
  smoke_summary <- tibble(
    cohort = cohort,
    integration = integration,
    sex = sex_stratum,
    resolution = resolution,
    samples = n_samples,
    cell_groups = n_distinct(dat_ace$cell_group),
    total_count = sum(dat_ace$count, na.rm = TRUE)
  )
  write_csv(smoke_summary, file.path(smoke_dir, "sccomp_smoke_summary.csv"))
  cat("Smoke summary written\n")
  quit(save = "no", status = 0)
}

if (n_samples < 20) {
  cat("SKIP: fewer than 20 samples (", n_samples, ") -- insufficient for sccomp\n")
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

if (!is.null(arm)) {
  # Per-arm AD-confounding mode: males-only, single primary tot_adverse_exp
  # composition model whose covariates match the DEG arm. AD_binary derived from
  # niareagansc when the arm needs it; interaction term added for MaleAceByAD.
  if (arm_spec_obj$needs_ad_binary) {
    dat_ace$AD_binary <- as.integer(dat_ace$niareagansc %in% c(1, 2))
  }
  missing_ad <- setdiff(arm_spec_obj$ad_covars, c(names(dat_ace), "AD_binary"))
  if (length(missing_ad) > 0) {
    stop("Metadata missing AD covariate(s) for arm ", arm, ": ",
         paste(missing_ad, collapse = ", "))
  }
  # scale continuous AD covariates (amylsqrt/tangsqrt) like age_death/pmi
  for (cv in arm_spec_obj$ad_covars) {
    if (cv != "AD_binary" && cv %in% names(dat_ace)) dat_ace <- scale_col(dat_ace, cv)
  }
  ad_terms <- arm_spec_obj$ad_covars
  covars <- paste(c("age_death", "pmi", ad_terms), collapse = " + ")
  # MaleAceByAD: add the AD_binary x phenotype interaction in the composition formula.
  arm_interaction <- isTRUE(arm_spec_obj$interaction)

  # sccomp needs complete covariates: drop samples missing any AD covariate
  # (e.g. ~2 males lack amylsqrt). Drop whole samples (all their cell_group rows).
  ad_cols_present <- intersect(c(ad_terms, "niareagansc"), names(dat_ace))
  if (length(ad_cols_present) > 0) {
    bad_samples <- dat_ace |>
      filter(if_any(all_of(ad_cols_present), is.na)) |>
      pull(sample) |>
      unique()
    if (length(bad_samples) > 0) {
      dat_ace <- dat_ace |> filter(!sample %in% bad_samples)
      cat("Dropped", length(bad_samples), "samples missing AD covariate(s) for arm", arm, "\n")
    }
  }

  models <- list(
    list(name = "primary", phenotype = "tot_adverse_exp", subdir = "primary", diff_var = FALSE)
  )
} else {
  covars <- if (sex_stratum == "all") {
    "age_death + msex + pmi + niareagansc"
  } else {
    "age_death + pmi + niareagansc"
  }
  arm_interaction <- FALSE

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

  for (column in c("age_death", "pmi", "niareagansc", "msex", phenotype,
                   "amylsqrt", "tangsqrt", "AD_binary")) {
    if (column %in% names(dat_sub)) {
      dat_sub[[column]] <- as.numeric(dat_sub[[column]])
    }
  }

  rhs <- paste0(phenotype, " + ", covars)
  if (isTRUE(arm_interaction)) {
    rhs <- paste0(rhs, " + AD_binary:", phenotype)
  }
  formula_comp <- as.formula(paste0("~ ", rhs))
  formula_var <- if (diff_var) as.formula(paste0("~ ", phenotype)) else ~ 1

  cat("  formula_composition:", deparse(formula_comp), "\n")
  cat("  formula_variability:", deparse(formula_var), "\n")

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  file_prefix <- paste0("sccomp_", phenotype, "_", sex_stratum, "_", resolution)

  result <- tryCatch({
    fit <- dat_sub |>
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

    cat("  Model fitted successfully\n")
    fit
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
    cat("  Results saved:", nrow(res_tbl), "cell types\n")
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

  outliers <- tryCatch({
    attr(result, "outliers")
  }, error = function(e) NULL)

  if (!is.null(outliers) && nrow(outliers) > 0) {
    write_csv(outliers, file.path(out_dir, paste0(file_prefix, "_outliers.csv")))
  }

  invisible(result)
}

cat("\n========================================\n")
cat("Running", length(models), "models\n")
cat("========================================\n")

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
