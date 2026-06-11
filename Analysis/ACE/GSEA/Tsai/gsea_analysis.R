#!/usr/bin/env Rscript
# gsea_analysis.R -- Ranked GSEA via WebGestaltR for ACE DEG results
#
# Usage:
#   Rscript gsea_analysis.R \
#     --deg-results-dir /path/to/deg/results \
#     --phenotype tot_adverse_exp \
#     --sex Female \
#     --output-dir /path/to/output \
#     [--smoke]
#
# Reads DESeq2 .rda files, computes gene ranks, and runs GSEA against 8
# pathway/network databases.  Produces per-cell-type RDS results, ranked gene
# CSVs, and a combined summary CSV.

library(WebGestaltR)
library(dplyr)

# ── GSEA databases ───────────────────────────────────────────────────────────

GSEA_DATABASES <- c(
  "geneontology_Biological_Process_noRedundant",
  "geneontology_Cellular_Component_noRedundant",
  "geneontology_Molecular_Function_noRedundant",
  "pathway_KEGG",
  "pathway_Panther",
  "pathway_Reactome",
  "pathway_Wikipathway",
  "network_Transcription_Factor_target"
)

# Short names for output files (strip common prefixes)
db_short_name <- function(db) {
  gsub("^geneontology_|^pathway_|^network_", "", db)
}

# ── Parse arguments ──────────────────────────────────────────────────────────

args <- commandArgs(trailingOnly = TRUE)
deg_results_dir <- NULL
phenotype       <- NULL
sex             <- NULL
output_dir      <- NULL
smoke_mode      <- FALSE
arm             <- NULL   # male AD-model arm suffix, e.g. "MaleContAD" (overrides --sex naming)
obj_name        <- NULL   # DESeq2 object name in the .rda (default "res")
celltypes_arg   <- NULL   # optional comma-list to restrict to specific cell types

i <- 1
while (i <= length(args)) {
  if (args[i] == "--deg-results-dir" && i < length(args)) {
    deg_results_dir <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--phenotype" && i < length(args)) {
    phenotype <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--sex" && i < length(args)) {
    sex <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--output-dir" && i < length(args)) {
    output_dir <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--arm" && i < length(args)) {
    arm <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--obj-name" && i < length(args)) {
    obj_name <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--celltypes" && i < length(args)) {
    celltypes_arg <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--smoke") {
    smoke_mode <- TRUE; i <- i + 1
  } else {
    i <- i + 1
  }
}

# ── Validate inputs ─────────────────────────────────────────────────────────

if (is.null(deg_results_dir)) stop("ERROR: --deg-results-dir is required.")
if (is.null(phenotype))       stop("ERROR: --phenotype is required.")
# --sex is required only for the baseline (non-arm) naming; male arms are males-only.
if (is.null(sex) && is.null(arm)) stop("ERROR: --sex (or --arm) is required.")
if (is.null(sex)) sex <- "Male"
if (is.null(output_dir))      stop("ERROR: --output-dir is required.")

if (!dir.exists(deg_results_dir)) {
  stop("ERROR: DEG results directory does not exist: ", deg_results_dir)
}

# The DESeq2 results object is named "res" by default; male AD arms with two
# contrasts (BinaryAD/AceByAD) store the ACE result as "res_ace".
if (is.null(obj_name)) {
  obj_name <- if (!is.null(arm) && arm %in% c("MaleBinaryAD", "MaleAceByAD")) "res_ace" else "res"
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== ACE GSEA Analysis ===\n")
cat("DEG results dir:", deg_results_dir, "\n")
cat("Phenotype:      ", phenotype, "\n")
cat("Sex:            ", sex, "\n")
cat("Output dir:     ", output_dir, "\n")
cat("Smoke mode:     ", smoke_mode, "\n")
cat("\n")

# ── Discover .rda files ──────────────────────────────────────────────────────

# Determine the filename label. Baseline sex-stratified files use "Fem"/"Male";
# male AD-model arms name the .rda by the ARM only (e.g. "_MaleContAD.rda",
# "_MaleBinaryAD.rda"). The two-contrast arms (BinaryAD/AceByAD) store BOTH
# results in that single .rda -- the ACE result is the object "res_ace"
# (handled by obj_name above). The "_ACEmain" suffix exists only on the CSVs,
# NOT the .rda, so the rda label is always just the arm name.
if (!is.null(arm)) {
  sex_label <- arm
} else {
  sex_label <- ifelse(sex == "Female", "Fem", "Male")
}
rda_pattern <- paste0("^deseqAnalysisACE_", phenotype, "_.*_", sex_label, "\\.rda$")
rda_files <- list.files(deg_results_dir, pattern = rda_pattern, full.names = TRUE)
# For male arms the suffix may contain "_" (e.g. "_ACEmain"); ensure we only match
# files ending in exactly that suffix (the pattern above already anchors with $).
cat("Using filename label (suffix):", sex_label, "  object:", obj_name, "\n")

if (length(rda_files) == 0) {
  stop("ERROR: No .rda files found matching pattern '", rda_pattern,
       "' in ", deg_results_dir)
}

# Optional restriction to specific cell types (comma list). The cell type is the
# filename segment between the phenotype prefix and the "_<sex_label>.rda" suffix.
if (!is.null(celltypes_arg)) {
  keep_cts <- trimws(strsplit(celltypes_arg, ",")[[1]])
  file_ct <- sub(paste0("_", sex_label, "\\.rda$"), "",
                 sub(paste0("^deseqAnalysisACE_", phenotype, "_"), "",
                     basename(rda_files)))
  rda_files <- rda_files[file_ct %in% keep_cts]
  if (length(rda_files) == 0) {
    stop("ERROR: --celltypes filter removed all files. Requested: ",
         paste(keep_cts, collapse = ", "))
  }
  cat("Restricted to cell types:", paste(keep_cts, collapse = ", "), "\n")
}

cat("Found", length(rda_files), ".rda file(s):\n")
for (f in rda_files) cat("  ", basename(f), "\n")
cat("\n")

# In smoke mode, only process the first file
if (smoke_mode) {
  rda_files <- rda_files[1]
  cat("[SMOKE] Limiting to first .rda file:", basename(rda_files[1]), "\n\n")
}

# ── Extract cell type from filename ──────────────────────────────────────────

extract_celltype <- function(filename, phenotype, sex_label) {
  # Pattern: deseqAnalysisACE_{phenotype}_{celltype}_{sex_label}.rda
  base <- tools::file_path_sans_ext(basename(filename))
  prefix <- paste0("deseqAnalysisACE_", phenotype, "_")
  suffix <- paste0("_", sex_label)
  ct <- sub(paste0("^", prefix), "", base)
  ct <- sub(paste0(suffix, "$"), "", ct)
  return(ct)
}

# ── Compute gene rank ────────────────────────────────────────────────────────

compute_gene_rank <- function(res_obj) {
  # res_obj is a DESeqResults S4 object
  df <- data.frame(
    gene           = rownames(res_obj),
    log2FoldChange = res_obj$log2FoldChange,
    pvalue         = res_obj$pvalue,
    stringsAsFactors = FALSE
  )

  # Compute signed rank statistic

  df$rank <- sign(df$log2FoldChange) * -log10(df$pvalue)

  # Remove genes with NA or infinite rank
  df <- df[is.finite(df$rank), ]

  # Sort descending by rank
  df <- df[order(-df$rank), ]

  return(data.frame(gene = df$gene, rank = df$rank, stringsAsFactors = FALSE))
}

# ── Run GSEA for a single cell type ─────────────────────────────────────────

run_gsea_celltype <- function(rda_path, celltype, databases, output_dir) {
  cat("── Processing cell type:", celltype, "──\n")

  # Load the DESeq2 results object
  env <- new.env()
  load(rda_path, envir = env)

  if (!exists(obj_name, envir = env)) {
    warning("Object '", obj_name, "' not found in ", basename(rda_path),
            ". Available: ", paste(ls(env), collapse = ", "), ". Skipping.")
    return(NULL)
  }

  res_obj <- get(obj_name, envir = env)
  cat("  Loaded", nrow(res_obj), "genes from", basename(rda_path), "\n")

  # Compute gene rank
  gsea_rank <- compute_gene_rank(res_obj)
  cat("  Ranked genes (after filtering):", nrow(gsea_rank), "\n")

  if (nrow(gsea_rank) < 50) {
    warning("Fewer than 50 ranked genes for ", celltype, ". Skipping.")
    return(NULL)
  }

  # Save ranked gene list
  rank_csv <- file.path(output_dir, paste0(celltype, "_ranked_genes.csv"))
  write.csv(gsea_rank, rank_csv, row.names = FALSE)
  cat("  Saved ranked genes:", basename(rank_csv), "\n")

  # Run GSEA for each database
  all_db_results <- list()

  skip_existing <- identical(Sys.getenv("ACE_GSEA_SKIP_EXISTING", unset = "0"), "1")

  for (db in databases) {
    db_short <- db_short_name(db)
    rds_path <- file.path(output_dir, paste0(celltype, "_", db_short, ".rds"))

    if (skip_existing && file.exists(rds_path) && file.info(rds_path)$size > 0) {
      cat("  Skip (exists): ", db_short, "\n")
      next
    }

    cat("  Running GSEA: ", db_short, " ... ")

    tryCatch({
      result <- WebGestaltR(
        enrichMethod      = "GSEA",
        organism           = "hsapiens",
        minNum             = 10,
        fdrThr             = 0.2,
        enrichDatabase     = db,
        enrichDatabaseType = "genesymbol",
        interestGene       = gsea_rank,
        interestGeneType   = "genesymbol",
        sigMethod          = "fdr",
        isOutput           = FALSE
      )

      if (is.null(result) || nrow(result) == 0) {
        cat("no significant results\n")
      } else {
        cat(nrow(result), "enriched sets\n")
        saveRDS(result, rds_path)

        # Annotate with metadata for summary
        result$cell_type <- celltype
        result$sex       <- sex
        result$database  <- db_short
        all_db_results[[db_short]] <- result
      }
    }, error = function(e) {
      cat("ERROR:", conditionMessage(e), "\n")
    })
  }

  return(all_db_results)
}

# ── Main loop ────────────────────────────────────────────────────────────────

databases_to_run <- GSEA_DATABASES
if (smoke_mode) {
  databases_to_run <- GSEA_DATABASES[1:min(2, length(GSEA_DATABASES))]
  cat("[SMOKE] Limiting to", length(databases_to_run), "databases\n\n")
}

all_results <- list()

for (rda_file in rda_files) {
  celltype <- extract_celltype(rda_file, phenotype, sex_label)
  ct_results <- run_gsea_celltype(rda_file, celltype, databases_to_run, output_dir)

  if (!is.null(ct_results)) {
    all_results <- c(all_results, ct_results)
  }
}

# ── Combine into summary CSV ────────────────────────────────────────────────

cat("\n── Building summary ──\n")

summary_rows <- lapply(all_results, function(df) {
  if (is.null(df) || nrow(df) == 0) return(NULL)

  # Select and rename columns to a standard schema
  df %>%
    select(
      cell_type,
      sex,
      database,
      geneSet,
      description,
      NES             = normalizedEnrichmentScore,
      pValue,
      FDR             = FDR,
      size,
      leadingEdgeNum
    )
})

summary_rows <- summary_rows[!sapply(summary_rows, is.null)]

# When restricted to a single cell type (per-CT fan-out), name the summary by
# cell type so concurrent per-CT jobs don't clobber a shared gsea_summary.csv.
# A separate merge step (merge_gsea_summaries.R) combines them.
summary_name <- if (!is.null(celltypes_arg) &&
                    length(trimws(strsplit(celltypes_arg, ",")[[1]])) == 1) {
  paste0("gsea_summary_", trimws(strsplit(celltypes_arg, ",")[[1]])[1], ".csv")
} else {
  "gsea_summary.csv"
}

if (length(summary_rows) > 0) {
  summary_df <- bind_rows(summary_rows)
  summary_path <- file.path(output_dir, summary_name)
  write.csv(summary_df, summary_path, row.names = FALSE)
  cat("Summary written:", summary_path, "\n")
  cat("  Total enriched sets:", nrow(summary_df), "\n")
  cat("  Cell types:", paste(unique(summary_df$cell_type), collapse = ", "), "\n")
  cat("  Databases:", paste(unique(summary_df$database), collapse = ", "), "\n")
} else {
  cat("No significant enrichment results across any cell type / database.\n")
  # Write empty summary so downstream scripts do not error
  summary_df <- data.frame(
    cell_type = character(), sex = character(), database = character(),
    geneSet = character(), description = character(), NES = numeric(),
    pValue = numeric(), FDR = numeric(), size = integer(),
    leadingEdgeNum = integer()
  )
  write.csv(summary_df, file.path(output_dir, summary_name), row.names = FALSE)
}

cat("\n=== GSEA analysis complete ===\n")
