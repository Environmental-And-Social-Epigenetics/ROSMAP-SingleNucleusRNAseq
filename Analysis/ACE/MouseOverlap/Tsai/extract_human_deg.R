#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
})

parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  out <- list(
    results_root = Sys.getenv("ACE_HUMAN_DEG_ROOT", unset = ""),
    output = "",
    integrations = "derived_batch,projid",
    phenotypes = "tot_adverse_exp,ace_aggregate,early_hh_ses",
    cell_types = "Exc,Inh,In-PV_Basket,In-PV_Chandelier"
  )
  i <- 1
  while (i <= length(args)) {
    if (args[i] == "--results-root" && i < length(args)) {
      out$results_root <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--output" && i < length(args)) {
      out$output <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--integrations" && i < length(args)) {
      out$integrations <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--phenotypes" && i < length(args)) {
      out$phenotypes <- args[i + 1]
      i <- i + 2
    } else if (args[i] == "--cell-types" && i < length(args)) {
      out$cell_types <- args[i + 1]
      i <- i + 2
    } else {
      stop("Unknown or incomplete argument: ", args[i])
    }
  }
  out
}

split_csv_arg <- function(x) {
  x <- trimws(unlist(strsplit(x, ",")))
  x[nzchar(x)]
}

integration_dirs <- function(results_root, integration) {
  if (integration %in% c("derived_batch", "batch")) {
    labels <- c("derived_batch", "batch")
  } else {
    labels <- integration
  }
  file.path(results_root, paste0("results_", labels))
}

extract_one <- function(path, integration, phenotype) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  if (!exists("res", envir = env)) {
    stop("No object named 'res' in ", path)
  }

  res <- get("res", envir = env)
  df <- as.data.frame(res)
  df$gene_symbol <- rownames(df)

  bn <- sub("\\.rda$", "", basename(path))
  prefix <- paste0("deseqAnalysisACE_", phenotype, "_")
  suffix <- sub(prefix, "", bn, fixed = TRUE)
  sex <- sub("^.*_(Fem|Male)$", "\\1", suffix)
  cell_type <- sub("_(Fem|Male)$", "", suffix)
  if (!sex %in% c("Fem", "Male")) {
    stop("Could not parse sex from filename: ", basename(path))
  }

  keep <- c("gene_symbol", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
  missing <- setdiff(keep, colnames(df))
  if (length(missing) > 0) {
    stop("Missing DESeq2 columns in ", path, ": ", paste(missing, collapse = ", "))
  }

  out <- df[keep]
  out$integration <- integration
  out$phenotype <- phenotype
  out$cell_type <- cell_type
  out$sex <- sex
  out$human_sig <- !is.na(out$padj) & out$padj < 0.05
  out$human_direction <- ifelse(
    out$log2FoldChange > 0, "ACE_up",
    ifelse(out$log2FoldChange < 0, "ACE_down", "zero")
  )

  out[, c(
    "integration", "phenotype", "cell_type", "sex", "gene_symbol",
    "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj",
    "human_sig", "human_direction"
  )]
}

args <- parse_args()
if (!nzchar(args$results_root)) {
  stop("ERROR: --results-root is required or ACE_HUMAN_DEG_ROOT must be set.")
}
if (!nzchar(args$output)) {
  stop("ERROR: --output is required.")
}

integrations <- split_csv_arg(args$integrations)
phenotypes <- split_csv_arg(args$phenotypes)
cell_types <- split_csv_arg(args$cell_types)
rows <- list()

cat("Human DEG results root:", args$results_root, "\n")
for (integration in integrations) {
  canonical_integration <- ifelse(integration == "batch", "derived_batch", integration)
  dirs <- integration_dirs(args$results_root, integration)
  result_dir <- dirs[file.exists(dirs)][1]
  if (is.na(result_dir)) {
    warning("No result directory found for integration ", integration)
    next
  }

  for (phenotype in phenotypes) {
    pheno_dir <- file.path(result_dir, phenotype)
    if (!dir.exists(pheno_dir)) {
      warning("No phenotype directory found: ", pheno_dir)
      next
    }

    files <- list.files(
      pheno_dir,
      pattern = paste0("^deseqAnalysisACE_", phenotype, "_.*\\.rda$"),
      full.names = TRUE
    )
    if (length(cell_types) > 0) {
      keep_patterns <- paste0("deseqAnalysisACE_", phenotype, "_", cell_types, "_(Fem|Male)\\.rda$")
      keep <- Reduce(`|`, lapply(keep_patterns, function(pattern) grepl(pattern, basename(files))))
      files <- files[keep]
    }
    cat("  ", canonical_integration, " / ", phenotype, ": ", length(files), " files\n", sep = "")
    for (f in files) {
      rows[[length(rows) + 1]] <- extract_one(f, canonical_integration, phenotype)
    }
  }
}

if (length(rows) == 0) {
  stop("No human DEG result rows extracted.")
}

flat <- do.call(rbind, rows)
dir.create(dirname(args$output), recursive = TRUE, showWarnings = FALSE)
write.csv(flat, args$output, row.names = FALSE)

cat("Wrote:", args$output, "\n")
cat("Rows:", nrow(flat), "\n")
cat("Significant rows:", sum(flat$human_sig, na.rm = TRUE), "\n")
