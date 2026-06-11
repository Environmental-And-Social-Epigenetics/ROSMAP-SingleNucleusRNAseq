# Downstream analysis suite — non-ANCOVA male ACE AD-model arms

Runs GSEA, module enrichment (hdWGCNA), SCENIC, and TF activity for each of the
5 non-ANCOVA male AD-model arms, over the signal-bearing cell types, so each
model's ACE signal can be compared at pathway / module / regulon / TF level.

- **Arms:** MaleNoADadj, MaleNiaReagan, MaleBinaryAD, MaleContAD, MaleAceByAD
- **Cell types:** In-PV_Basket, Ast, Mic, Oli, OPC, Inh (males only)
- **Phenotype:** tot_adverse_exp; **integration:** derived_batch
- **DEG (A):** already complete — used as input, not re-run.

## Per-arm covariate model (single source of truth)

`_shared/arm_covariates.{sh,py,R}` defines, per arm, the AD covariates added to the
baseline `tot_adverse_exp + age_death + pmi`, the DEG `.rda` object/suffix, and the
AD_binary derivation. All four analyses import it so association formulas match each
DEG arm exactly:

| Arm | extra covariates |
|---|---|
| MaleNoADadj | — |
| MaleNiaReagan | niareagansc |
| MaleBinaryAD | AD_binary |
| MaleContAD | amylsqrt + tangsqrt |
| MaleAceByAD | AD_binary + AD_binary:phenotype |

## Run

```bash
# everything (B GSEA, C module, D SCENIC, E TF):
bash run_downstream_male_arms.sh
# subset:
ANALYSES="B E" bash run_downstream_male_arms.sh
# dry run / reuse heavy builds:
SMOKE_FLAG=--smoke bash run_downstream_male_arms.sh
SKIP_BUILD=1 bash run_downstream_male_arms.sh
```

Or per analysis: `GSEA/Tsai/aceGseaT_male_arms.sh`,
`hdWGCNA/Tsai/aceWgcnaT_male_arms.sh`, `SCENIC/Tsai/aceScenicT_male_arms.sh`,
`TFActivity/Tsai/aceTfActT_male_arms.sh`.

## Build-once / associate-per-arm

The phenotype-independent heavy step runs ONCE per cell type (males); only the cheap
per-arm association re-runs. Existing heavy outputs are reused.

- **SCENIC:** GRN→cisTarget→AUCell once (`results_derived_batch/.../Male_<CT>/auc_matrix.csv`,
  mostly already present; only OPC missing). `scenic_associate.py` runs the per-arm OLS over
  the cached AUCell (micropool id `<projid>_pool<N>`). Stage 2 gated `afterok` on Stage 1.
- **hdWGCNA:** module detection once (`module_eigengenes.csv` + `module_assignments.csv`).
  `wgcna_associate.R` runs the per-arm module-trait regression (ME ~ pheno + arm covars) and
  module-DEG overlap (vs that arm's ACE DEGs). Stage 2 gated `afterok`.
- **GSEA / TF:** no heavy stage. GSEA reads the arm's DEG `.rda` (object `res`/`res_ace`).
  TF recomputes decoupler activities (cheap) then applies the arm OLS.

## Outputs

`${ACE_OUTPUT_ROOT}/<ANALYSIS>/Tsai/results_derived_batch_<ARM>/tot_adverse_exp/...`
- GSEA: `gsea_summary.csv` + per-cell-type ranked genes / per-DB RDS.
- TF: `tf_summary.csv`, `pathway_summary.csv`, per-CT `tf_{dorothea,collectri}_*` / `pathway_progeny_*`.
- SCENIC: `<CT>/regression_results.csv` (regulon ~ ACE, arm-adjusted).
- hdWGCNA: `Male_<CT>/module_trait_correlations.csv` + `module_deg_overlap.csv`.

## Notes / caveats
- Partition `pi_lhtsai,pi_manoli` only. Each job activates its env and calls the env's
  Rscript/python by full path (avoids login-shell PATH shadowing).
- `Inh` (DEG/GSEA/TF label) == `broad_Inh` (SCENIC/WGCNA h5ad/dir). Mapping handled in launchers.
- `wgcna_analysis.R` DEG-overlap object lookup was fixed to use `res`/`res_ace`
  (older fixtures used resF/resM).
- For cross-arm reading, the headline comparison is whether top pathways/regulons/modules
  persist from MaleNoADadj → MaleContAD (mirrors the DEG attenuation story).
