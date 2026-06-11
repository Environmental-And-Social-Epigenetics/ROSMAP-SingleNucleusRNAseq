#!/usr/bin/env Rscript
# =============================================================================
# Per-arm WGCNA module association (no module re-detection).
#
# hdWGCNA module detection is phenotype-independent and already wrote, per cell
# type (males):
#   module_eigengenes.csv   (per-patient module eigengenes + projid)
#   module_assignments.csv  (gene -> module)
# This script loads those, then for a given male AD-model arm:
#   1. module-trait association: ME ~ phenotype + age_death + pmi [+ arm AD covars]
#      [+ AD_binary:phenotype], reporting the phenotype term (partial assoc.).
#   2. module-DEG overlap: hypergeometric overlap of each module's genes with the
#      arm's ACE DEGs (padj < 0.05) from that arm's DEG .rda.
#
# Usage:
#   Rscript wgcna_associate.R \
#     --modules-dir   .../results_derived_batch/tot_adverse_exp/Male_<CT> \
#     --pheno-csv     $ACE_SCORES_CSV \
#     --arm           MaleContAD \
#     --phenotype     tot_adverse_exp \
#     --cell-type     <CT> \
#     --deg-results-dir .../results_derived_batch_<ARM>_AllCellTypes/tot_adverse_exp \
#     --output-dir    .../results_derived_batch_<ARM>/tot_adverse_exp/Male_<CT>
# =============================================================================

suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
modules_dir <- NULL; pheno_csv <- NULL; arm <- NULL
phenotype <- "tot_adverse_exp"; cell_type <- NULL
deg_results_dir <- NULL; output_dir <- NULL
i <- 1
while (i <= length(args)) {
  k <- args[i]
  if (k == "--modules-dir")           { i <- i + 1; modules_dir <- args[i] }
  else if (k == "--pheno-csv")        { i <- i + 1; pheno_csv <- args[i] }
  else if (k == "--arm")              { i <- i + 1; arm <- args[i] }
  else if (k == "--phenotype")        { i <- i + 1; phenotype <- args[i] }
  else if (k == "--cell-type")        { i <- i + 1; cell_type <- args[i] }
  else if (k == "--deg-results-dir")  { i <- i + 1; deg_results_dir <- args[i] }
  else if (k == "--output-dir")       { i <- i + 1; output_dir <- args[i] }
  else stop(paste("Unknown argument:", k))
  i <- i + 1
}
stopifnot(!is.null(modules_dir), !is.null(pheno_csv), !is.null(arm),
          !is.null(cell_type), !is.null(output_dir))

# Shared per-arm covariate spec
this_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(this_dir) || !nzchar(this_dir)) this_dir <- "."
shared_R <- file.path(this_dir, "..", "..", "_shared", "arm_covariates.R")
if (!file.exists(shared_R)) {
  # Portable fallback: resolve from REPO_ROOT (exported by config/paths.sh).
  repo_root <- Sys.getenv("REPO_ROOT", "")
  if (nzchar(repo_root)) {
    shared_R <- file.path(repo_root, "Analysis", "ACE", "_shared", "arm_covariates.R")
  }
}
source(shared_R)
spec <- arm_spec(arm)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

me_path  <- file.path(modules_dir, "module_eigengenes.csv")
mod_path <- file.path(modules_dir, "module_assignments.csv")
if (!file.exists(me_path))  stop("ERROR: module_eigengenes.csv not found in ", modules_dir)
if (!file.exists(mod_path)) stop("ERROR: module_assignments.csv not found in ", modules_dir)

me_patient <- read.csv(me_path, check.names = FALSE)
modules    <- read.csv(mod_path, check.names = FALSE)

# pid column in the MEs file (projid or patient_id)
pid_col <- if ("projid" %in% colnames(me_patient)) "projid" else
           if ("patient_id" %in% colnames(me_patient)) "patient_id" else colnames(me_patient)[1]
me_patient[[pid_col]] <- as.character(me_patient[[pid_col]])

# Phenotype + covariates (males already implied by the Male_<CT> module build)
pheno <- read.csv(pheno_csv)
pheno$projid <- as.character(pheno$projid)
pheno <- pheno[!duplicated(pheno$projid), ]
if (spec$needs_ad_binary) pheno$AD_binary <- as.integer(pheno$niareagansc %in% c(1, 2))

cov_cols <- unique(c(phenotype, "age_death", "pmi", spec$ad_covars,
                     if (spec$interaction) "AD_binary" else NULL))
miss <- setdiff(cov_cols, colnames(pheno))
if (length(miss) > 0) stop("ERROR: phenotype CSV missing columns for arm ", arm, ": ",
                           paste(miss, collapse = ", "))

# Re-merge arm covariates onto the per-patient MEs (drop any stale covariate cols first)
me_patient <- me_patient[, !colnames(me_patient) %in% setdiff(cov_cols, pid_col), drop = FALSE]
me_patient <- merge(me_patient, pheno[, c("projid", cov_cols)],
                    by.x = pid_col, by.y = "projid")

# --- Module-trait association (partial, arm-covariate-adjusted) --------------
me_cols <- grep("^ME", colnames(me_patient), value = TRUE)
me_cols <- me_cols[me_cols != "MEgrey"]
rhs <- arm_rhs(arm, phenotype)   # phenotype + age_death + pmi [+ ad covars] [+ interaction]

trait_rows <- list()
for (me_col in me_cols) {
  module_name <- sub("^ME", "", me_col)
  n_genes <- sum(modules$module == module_name, na.rm = TRUE)
  df <- me_patient[, c(me_col, cov_cols)]
  df[[phenotype]] <- as.numeric(df[[phenotype]])
  df <- df[stats::complete.cases(df), , drop = FALSE]
  if (nrow(df) < 5) next
  fml <- as.formula(paste0("`", me_col, "` ~ ", rhs))
  fit <- tryCatch(stats::lm(fml, data = df), error = function(e) NULL)
  if (is.null(fit)) next
  co <- summary(fit)$coefficients
  if (!phenotype %in% rownames(co)) next
  trait_rows[[length(trait_rows) + 1]] <- data.frame(
    module = module_name, n_genes = n_genes,
    coef = co[phenotype, "Estimate"],
    stderr = co[phenotype, "Std. Error"],
    pvalue = co[phenotype, "Pr(>|t|)"],
    n_obs = nrow(df), stringsAsFactors = FALSE
  )
}
trait_results <- if (length(trait_rows)) do.call(rbind, trait_rows) else
  data.frame(module = character(), n_genes = integer(), coef = numeric(),
             stderr = numeric(), pvalue = numeric(), n_obs = integer())
if (nrow(trait_results)) {
  trait_results$padj <- p.adjust(trait_results$pvalue, method = "BH")
  trait_results <- trait_results[order(trait_results$pvalue), ]
}
write.csv(trait_results, file.path(output_dir, "module_trait_correlations.csv"), row.names = FALSE)
cat("Module-trait (arm=", arm, "): ", sum(trait_results$padj < 0.05, na.rm = TRUE),
    " modules padj<0.05\n", sep = "")

# --- Module-DEG overlap (arm's ACE DEGs) ------------------------------------
if (!is.null(deg_results_dir) && dir.exists(deg_results_dir)) {
  deg_ct <- sub("^broad_", "", cell_type)
  rda <- file.path(deg_results_dir,
                   paste0("deseqAnalysisACE_", phenotype, "_", deg_ct, "_", spec$deg_suffix, ".rda"))
  if (file.exists(rda)) {
    env <- new.env(); load(rda, envir = env)
    res_obj <- if (exists(spec$deg_obj, envir = env)) get(spec$deg_obj, envir = env) else NULL
    if (!is.null(res_obj)) {
      deg_genes <- rownames(res_obj)[!is.na(res_obj$padj) & res_obj$padj < 0.05]
      all_genes <- unique(modules$gene_name)
      n_total <- length(all_genes)
      n_deg_total <- length(intersect(deg_genes, all_genes))
      ov <- list()
      for (mod in unique(modules$module)) {
        if (mod == "grey") next
        mg <- modules$gene_name[modules$module == mod]
        n_mod <- length(mg); n_ov <- length(intersect(mg, deg_genes))
        hp <- phyper(n_ov - 1, n_deg_total, n_total - n_deg_total, n_mod, lower.tail = FALSE)
        ov[[length(ov) + 1]] <- data.frame(
          module = mod, n_module_genes = n_mod, n_deg_in_module = n_ov,
          fraction_deg = n_ov / max(n_mod, 1), enrichment_p = hp, stringsAsFactors = FALSE)
      }
      if (length(ov)) {
        ov <- do.call(rbind, ov)
        ov$enrichment_padj <- p.adjust(ov$enrichment_p, method = "BH")
        write.csv(ov, file.path(output_dir, "module_deg_overlap.csv"), row.names = FALSE)
        cat("Module-DEG overlap: ", sum(ov$enrichment_padj < 0.05), " modules enriched (",
            n_deg_total, " DEGs in universe)\n", sep = "")
      }
    } else {
      cat("WARN: object '", spec$deg_obj, "' not in ", basename(rda), "\n", sep = "")
    }
  } else {
    cat("WARN: DEG rda not found: ", rda, "\n", sep = "")
  }
}

cat("Done: wgcna_associate ", cell_type, " / ", arm, "\n", sep = "")
