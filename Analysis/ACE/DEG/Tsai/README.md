# ACE DEG - Tsai

Canonical entry point:

```bash
bash aceDegT.sh
```

This launcher:

1. builds per-cell-type split `.h5ad` inputs from the canonical Tsai integrated object
2. runs pseudobulk DEG for the official ACE phenotype trio
3. writes outputs to `${ANALYSIS_OUTPUT_ROOT}/ACE/DEG/Tsai`

Useful helpers:

- `prep_celltype_splits.py`
- `run_deg.sh`
- `smoke_test.sh`

Canonical integration: `derived_batch`

Optional sensitivity integration: `projid`

Legacy exploratory helpers live under `legacy/` and are not part of the
supported workflow or smoke-test surface.

## Male AD-confounding model suite

Males-only (`msex == 1`) DEG arms for `tot_adverse_exp`, testing whether the ACE
signal is independent of Alzheimer's pathology. All arms share the same N=115 male
cohort (`niareagansc %in% {1,2,3,4}`) so the comparison contrasts the *adjustment*,
not the *sample*. Baseline covariates: `age_death + pmi`.

| Arm | Design (`+ tot_adverse_exp`) | Output dir suffix |
| --- | --- | --- |
| MaleNoADadj | `~ age_death + pmi` | `_MaleNoADadj_AllCellTypes` |
| **MaleNiaReagan** (Model 0) | `+ niareagansc` (raw ordinal) | `_MaleNiaReagan_AllCellTypes` |
| MaleBinaryAD (Model 2) | `+ AD_binary` (niareagansc∈{1,2}) | `_MaleBinaryAD_AllCellTypes` |
| MaleContAD (Model 1) | `+ amylsqrt + tangsqrt` | `_MaleContAD_AllCellTypes` |
| MaleAceByAD | `+ AD_binary + AD_binary:ACE` | `_MaleAceByAD_AllCellTypes` |
| **MaleAncovaAD** (Model 3) | per AD var `v ∈ {niareagansc, tangsqrt, braaksc, amylsqrt}`: `+ v` | `_MaleAncovaAD_AllCellTypes` |

`broad_Exc` is excluded (101 GB > R sparse-matrix limit); the 7 Ex subtypes provide
excitatory coverage. Each arm prints a donor-level collinearity diagnostic
(`cor(ACE, AD var)` + VIF of ACE) to stdout.

Launch the NEW arms (MaleNiaReagan + MaleAncovaAD) efficiently — one SLURM job per
cell type, all concurrent — then auto-run the cross-model comparison:

```bash
bash run_male_ad_models.sh                 # full run
SMOKE_FLAG=--smoke bash run_male_ad_models.sh   # dry run (load + pseudobulk dims)
RERUN_EXISTING=1 bash run_male_ad_models.sh      # also re-fire the 4 existing arms
```

Cross-model comparison: `report/summarize_ACE_AD_models.sh` →
`results_derived_batch_AD_comparison/` (`deg_counts.csv`, `lfc_correlations.csv`,
`ace_attenuation.csv`, `ace_vs_ad_beta_concordance.csv`, top-N tables, scatter PNGs).
