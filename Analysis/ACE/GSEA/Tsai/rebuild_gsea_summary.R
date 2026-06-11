#!/usr/bin/env Rscript
# Rebuild gsea_summary.csv for an arm from the per-(cell_type x database) .rds
# files on disk. Needed because the in-job summary only includes databases that
# ran fresh (ACE_GSEA_SKIP_EXISTING=1 skips already-computed DBs, so they were
# not re-aggregated). The .rds files hold the full results; this just collects them.
#
# Usage: Rscript rebuild_gsea_summary.R --results-dir <arm>/tot_adverse_exp [--sex Male]

suppressPackageStartupMessages(library(dplyr))
args <- commandArgs(trailingOnly = TRUE)
results_dir <- NULL; sex <- "Male"
i <- 1
while (i <= length(args)) {
  if (args[i] == "--results-dir") { i <- i + 1; results_dir <- args[i] }
  else if (args[i] == "--sex")    { i <- i + 1; sex <- args[i] }
  i <- i + 1
}
stopifnot(!is.null(results_dir))
if (!dir.exists(results_dir)) stop("ERROR: dir not found: ", results_dir)

rds_files <- list.files(results_dir, pattern = "\\.rds$", full.names = TRUE)
# A per-DB result file is named <celltype>_<db_short>.rds (db_short has no spaces).
# Ranked-gene CSVs are .csv, so .rds here are all enrichment results.
if (length(rds_files) == 0) { cat("No .rds files in ", results_dir, "\n"); quit(save="no", status=0) }

db_shorts <- c("Biological_Process_noRedundant","Cellular_Component_noRedundant",
               "Molecular_Function_noRedundant","KEGG","Panther","Reactome",
               "Wikipathway","Transcription_Factor_target")

rows <- list()
for (f in rds_files) {
  base <- sub("\\.rds$", "", basename(f))
  # strip the db_short suffix to recover cell type
  db <- NA_character_; ct <- base
  for (d in db_shorts) {
    if (endsWith(base, paste0("_", d))) { db <- d; ct <- substr(base, 1, nchar(base) - nchar(d) - 1); break }
  }
  if (is.na(db)) next
  r <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(r) || !is.data.frame(r) || nrow(r) == 0) next
  need <- c("geneSet","description","normalizedEnrichmentScore","pValue","FDR","size","leadingEdgeNum")
  if (!all(need %in% colnames(r))) next
  rows[[length(rows)+1]] <- data.frame(
    cell_type = ct, sex = sex, database = db,
    geneSet = r$geneSet, description = r$description,
    NES = r$normalizedEnrichmentScore, pValue = r$pValue, FDR = r$FDR,
    size = r$size, leadingEdgeNum = r$leadingEdgeNum,
    stringsAsFactors = FALSE)
}

out <- file.path(results_dir, "gsea_summary.csv")
if (length(rows) == 0) {
  empty <- data.frame(cell_type=character(),sex=character(),database=character(),
    geneSet=character(),description=character(),NES=numeric(),pValue=numeric(),
    FDR=numeric(),size=integer(),leadingEdgeNum=integer())
  write.csv(empty, out, row.names = FALSE)
  cat("Rebuilt (empty): ", out, "\n"); quit(save="no", status=0)
}
df <- bind_rows(rows)
df <- df[order(df$cell_type, df$database, df$FDR), ]
write.csv(df, out, row.names = FALSE)
cat("Rebuilt: ", out, " (", nrow(df), " enriched sets across ",
    length(unique(df$cell_type)), " cell types, ",
    length(unique(df$database)), " databases)\n", sep = "")
