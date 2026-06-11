#!/usr/bin/env Rscript
# Dump DESeq2 result tables from .rda files to CSV. Called by integrate.py.
#
# Usage:
#   Rscript _dump_deg_tables.R --deg-dir /path/results_derived_batch/<phenotype>/ --out-dir /path/cache

args <- commandArgs(trailingOnly = TRUE)

deg_dir <- NULL
out_dir <- NULL

i <- 1
while (i <= length(args)) {
    if (args[i] == "--deg-dir") { i <- i + 1; deg_dir <- args[i] }
    else if (args[i] == "--out-dir") { i <- i + 1; out_dir <- args[i] }
    else stop(paste("Unknown arg:", args[i]))
    i <- i + 1
}
stopifnot(!is.null(deg_dir), !is.null(out_dir))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rda_files <- list.files(deg_dir, pattern = "^deseqAnalysisACE_.*\\.rda$", full.names = TRUE)
cat("Found", length(rda_files), ".rda files in", deg_dir, "\n")

for (rda in rda_files) {
    env <- new.env()
    load(rda, envir = env)
    obj <- ls(env)
    if (length(obj) != 1) {
        cat("  skip", basename(rda), "(expected one object, got", length(obj), ")\n")
        next
    }
    res <- get(obj[1], envir = env)
    df <- tryCatch(as.data.frame(res), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) {
        cat("  skip", basename(rda), "(could not coerce to data.frame)\n")
        next
    }
    if (!"gene" %in% colnames(df)) {
        df$gene <- rownames(df)
    }
    out_csv <- file.path(out_dir, sub("\\.rda$", ".csv", basename(rda)))
    write.csv(df, out_csv, row.names = FALSE)
    cat("  wrote", basename(out_csv), "(", nrow(df), "genes )\n")
}
