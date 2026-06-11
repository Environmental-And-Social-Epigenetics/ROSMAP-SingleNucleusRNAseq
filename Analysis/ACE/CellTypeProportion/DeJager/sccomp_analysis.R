# Thin DeJager wrapper around the canonical, cohort-aware sccomp analysis.
# Mirrors the DEG parity pattern (aceDegDJ.Rscript source()s aceDegT.Rscript).
# Default the cohort to dejager (-> integration library_id) unless the caller
# already passed --cohort / set ACE_PROP_COHORT; then source the single source
# of truth in ../Tsai/sccomp_analysis.R, which honours the same --cohort,
# --integration, --sex, --resolution, --output-root, --arm and --smoke args.
if (!nzchar(Sys.getenv("ACE_PROP_COHORT")) &&
    !any(commandArgs(trailingOnly = TRUE) == "--cohort")) {
  Sys.setenv(ACE_PROP_COHORT = "dejager")
}

# Resolve this wrapper's directory robustly. Under `Rscript file.R` (how the
# launchers invoke it) the path comes from the --file= arg; fall back to the
# source() ofile, then to the working directory.
this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
script_dir <- if (!is.na(this_file) && nzchar(this_file)) {
  dirname(normalizePath(this_file, mustWork = FALSE))
} else {
  of <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(of) && nzchar(of)) dirname(normalizePath(of, mustWork = FALSE)) else getwd()
}

source(file.path(script_dir, "..", "Tsai", "sccomp_analysis.R"))
