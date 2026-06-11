#' Shared ACE phenotype loader for Tsai / DeJager R workflows.
#'
#' The canonical CSV pools Tsai + DeJager and has duplicate projids.
#' Use load_ace_for_projids() before any projid-indexed join/merge.

ACE_REQUIRED_COLUMNS <- c(
  "projid", "tot_adverse_exp", "early_hh_ses",
  "msex", "age_death", "pmi", "niareagansc"
)

.add_ace_aggregate <- function(pheno) {
  if ("ace_aggregate" %in% colnames(pheno)) return(pheno)
  if (!all(c("tot_adverse_exp", "early_hh_ses") %in% colnames(pheno))) return(pheno)
  z_adv <- scale(pheno$tot_adverse_exp)[, 1]
  z_ses <- scale(pheno$early_hh_ses)[, 1]
  pheno$ace_aggregate <- z_adv - z_ses
  pheno
}

load_ace_raw <- function(pheno_csv) {
  pheno <- read.csv(pheno_csv, stringsAsFactors = FALSE, check.names = FALSE)
  # Drop the R-written rowname column ("" header becomes "X" under base R's
  # read.csv). Carrying it breaks dedup because every row has a unique index.
  name_col_candidates <- c("X", "")
  drop_cols <- intersect(name_col_candidates, colnames(pheno))
  if (length(drop_cols) > 0) pheno <- pheno[, setdiff(colnames(pheno), drop_cols), drop = FALSE]
  missing_cols <- setdiff(ACE_REQUIRED_COLUMNS, colnames(pheno))
  if (length(missing_cols) > 0) {
    stop(
      "Phenotype CSV missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  pheno$projid <- as.character(pheno$projid)
  .add_ace_aggregate(pheno)
}

#' Load phenotypes restricted to the supplied projids (e.g. cohort obs_names).
#'
#' Filters to rows whose projid is in `projids`, drops exact duplicate rows,
#' and stops if any projid remains duplicated after dedup.
#'
#' @param pheno_csv path to ACE phenotype CSV
#' @param projids character vector of projids to keep
#' @param require_all logical; stop if any projid in `projids` is missing
#' @return data.frame with rownames == projids (ordered to match `projids`)
load_ace_for_projids <- function(pheno_csv, projids, require_all = FALSE) {
  projids <- as.character(projids)
  pheno <- load_ace_raw(pheno_csv)

  sub <- pheno[pheno$projid %in% projids, , drop = FALSE]
  sub <- unique(sub)

  dup_ids <- sub$projid[duplicated(sub$projid) | duplicated(sub$projid, fromLast = TRUE)]
  if (length(dup_ids) > 0) {
    offenders <- unique(dup_ids)
    stop(
      "Phenotype CSV has conflicting rows for projid(s): ",
      paste(utils::head(offenders, 10), collapse = ", "),
      if (length(offenders) > 10) " ..." else "",
      " (exact duplicates were dropped; these remain after dedup)."
    )
  }

  rownames(sub) <- sub$projid
  sub <- sub[projids, , drop = FALSE]
  rownames(sub) <- projids

  if (require_all) {
    missing_pid <- projids[rowSums(!is.na(sub)) == 0]
    if (length(missing_pid) > 0) {
      stop(
        length(missing_pid), " projid(s) from cohort not found in phenotype CSV: ",
        paste(utils::head(missing_pid, 10), collapse = ", "),
        if (length(missing_pid) > 10) " ..." else ""
      )
    }
  }
  sub
}
