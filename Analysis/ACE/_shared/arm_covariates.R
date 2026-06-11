# Single source of truth for the non-ANCOVA male ACE AD-model arms (R side).
# Mirrors arm_covariates.py / arm_covariates.sh. Used by GSEA (DEG suffix/object) and
# WGCNA (module-trait covariate set). All arms males-only, phenotype tot_adverse_exp,
# integration derived_batch. AD_binary = 1 if niareagansc in {1,2}.

ARM_SPECS <- list(
  MaleNoADadj = list(
    ad_covars = character(0),
    interaction = FALSE, needs_ad_binary = FALSE,
    deg_obj = "res", deg_suffix = "MaleNoADadj", nia_filter = c(1, 2, 3, 4)
  ),
  MaleNiaReagan = list(
    ad_covars = c("niareagansc"),
    interaction = FALSE, needs_ad_binary = FALSE,
    deg_obj = "res", deg_suffix = "MaleNiaReagan", nia_filter = c(1, 2, 3, 4)
  ),
  MaleBinaryAD = list(
    ad_covars = c("AD_binary"),
    interaction = FALSE, needs_ad_binary = TRUE,
    deg_obj = "res_ace", deg_suffix = "MaleBinaryAD", nia_filter = c(1, 2, 3, 4)
  ),
  MaleContAD = list(
    ad_covars = c("amylsqrt", "tangsqrt"),
    interaction = FALSE, needs_ad_binary = FALSE,
    deg_obj = "res", deg_suffix = "MaleContAD", nia_filter = c(1, 2, 3, 4)
  ),
  MaleAceByAD = list(
    ad_covars = c("AD_binary"),
    interaction = TRUE, needs_ad_binary = TRUE,
    deg_obj = "res_ace", deg_suffix = "MaleAceByAD", nia_filter = c(1, 2, 3, 4)
  )
)

NON_ANCOVA_ARMS <- names(ARM_SPECS)
SIGNAL_CELLTYPES <- c("In-PV_Basket", "Ast", "Mic", "Oli", "OPC", "Inh")

arm_spec <- function(arm) {
  if (!arm %in% names(ARM_SPECS)) {
    stop("Unknown arm '", arm, "'. Known: ", paste(NON_ANCOVA_ARMS, collapse = ", "))
  }
  ARM_SPECS[[arm]]
}

# Build a model-formula RHS string for an arm:
#   <pheno> + age_death + pmi [+ ad_covars...] [+ AD_binary:<pheno>]
arm_rhs <- function(arm, phenotype) {
  spec <- arm_spec(arm)
  terms <- c(phenotype, "age_death", "pmi", spec$ad_covars)
  if (spec$interaction) terms <- c(terms, paste0("AD_binary:", phenotype))
  paste(terms, collapse = " + ")
}
