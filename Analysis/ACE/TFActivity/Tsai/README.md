# ACE Transcription Factor & Pathway Activity Analysis (Tsai Cohort)

Infers TF activity and pathway activity from pseudobulked snRNA-seq data
using [decoupler](https://decoupler-py.readthedocs.io/), then tests
associations with adverse childhood experience (ACE) phenotypes via OLS
regression with demographic covariates.

## Method

### Overview

1. **Input**: per-cell-type h5ad files (raw counts) produced by the DEG
   preprocessing step.
2. **Pseudobulk**: raw counts are summed per patient to create one
   observation per patient per cell type.
3. **Normalization**: CPM (counts per million) followed by log1p.
4. **Activity inference**: decoupler's multivariate linear model (MLM)
   estimates TF/pathway activity scores for each patient.
5. **Statistical testing**: OLS regression of each activity score against
   the ACE phenotype, with covariates (age_death, pmi, niareagansc).
6. **Multiple testing**: Benjamini-Hochberg FDR correction across all
   TFs/pathways within each cell type and sex stratum.

### Databases

| Database | Type | Description |
|----------|------|-------------|
| **DoRothEA** (A, B, C) | TF regulons | Curated TF-target interactions at confidence levels A-C. High confidence, moderate coverage. |
| **CollecTRI** | TF regulons | Community-collected TF-target interactions. Broader coverage than DoRothEA. |
| **PROGENy** (top 500) | Pathway footprints | Pathway-responsive genes weighted by consistency across perturbation experiments. 14 canonical signaling pathways. |

### Sex stratification

All analyses are run separately for males (msex=1) and females (msex=0).
Cell type/sex strata with fewer than 5 patients are skipped.

### Regression model

```
activity ~ phenotype + age_death + pmi + niareagansc
```

Where `phenotype` is one of: `tot_adverse_exp`, `early_hh_ses`,
`ace_aggregate`.

## Interpretation Guide

### Coefficients

- **Positive coefficient**: higher ACE score is associated with **higher**
  TF/pathway activity. For TFs, this means the TF's target gene program is
  more active in patients with more adversity.
- **Negative coefficient**: higher ACE score is associated with **lower**
  TF/pathway activity.

### Epigenetic TFs

The visualization script highlights known epigenetic regulators:

| Category | TFs |
|----------|-----|
| Histone deacetylases | HDAC1-11 |
| DNA methyltransferases | DNMT1, DNMT3A, DNMT3B |
| Polycomb | EZH2, SUZ12 |
| Chromatin remodelers | REST, SIN3A, CTCF |
| Histone demethylases | KDM5A, KDM5B |
| Histone acetyltransferases | KAT2A, KAT2B, EP300, CREBBP |

These are of particular interest because ACE-related epigenetic
reprogramming (e.g., altered DNA methylation, histone modifications) is a
hypothesized mechanism linking early adversity to late-life
neurodegeneration.

### Hub TFs

The convergence analysis identifies TFs that are significant (FDR < 0.05)
in 3 or more cell types with a consistent direction of effect. These
"hub" regulators suggest cell-type-general transcriptional programs
affected by ACE.

## Entry Point

```bash
# Submit all three phenotype jobs
bash aceTfActT.sh              # default integration = derived_batch
bash aceTfActT.sh projid       # alternative integration
```

Individual steps:

```bash
# Single phenotype via SLURM
sbatch run_tf_activity.sh derived_batch tot_adverse_exp

# Local smoke test (no SLURM)
bash smoke_test.sh

# Post-hoc convergence (interactive, no SLURM)
python tf_convergence_analysis.py \
  --results-dir ${ACE_OUTPUT_ROOT}/TFActivity/Tsai/results_derived_batch/tot_adverse_exp \
  --output-dir  ${ACE_OUTPUT_ROOT}/TFActivity/Tsai/convergence_derived_batch/tot_adverse_exp \
  --phenotype tot_adverse_exp
```

## Output Structure

```text
${ACE_OUTPUT_ROOT}/TFActivity/Tsai/
  logs/
    <jobid>_<phenotype>_<integration>.out
    <jobid>_<phenotype>_<integration>.err
  results_<integration>/<phenotype>/
    tf_dorothea_<celltype>_<sex>.csv
    tf_collectri_<celltype>_<sex>.csv
    pathway_progeny_<celltype>_<sex>.csv
    tf_summary.csv           # all TF results combined
    pathway_summary.csv      # all pathway results combined
  figures_<integration>/<phenotype>/
    tf_heatmap_<sex>.pdf
    tf_volcano_<sex>.pdf
    epigenetic_tf_<sex>.pdf
    hub_tf_barplot_<sex>.pdf
    pathway_heatmap_<sex>.pdf
  convergence_<integration>/<phenotype>/
    hub_tf_table.csv
    cross_database_validation.csv
    convergence_heatmap_<sex>.pdf
    scenic_overlap.csv       # if SCENIC results are available
```

### Result CSV Columns

| Column | Description |
|--------|-------------|
| `name` | TF or pathway name |
| `coef` | OLS coefficient for the ACE phenotype term |
| `stderr` | Standard error of the coefficient |
| `tstat` | t-statistic |
| `pvalue` | Raw p-value |
| `padj` | BH-corrected FDR |
| `n_targets` | Number of target genes in the regulon/pathway |
| `n_obs` | Number of patients in the regression |
| `source_db` | Database (DoRothEA, CollecTRI, or PROGENy) |
| `cell_type` | Cell type label |
| `sex` | Male or Female |
| `phenotype` | ACE phenotype tested |

## Environment

Uses the `decoupler_compat` conda environment:
- decoupler 1.8.0
- scanpy, anndata, pandas, statsmodels, matplotlib, seaborn

## Files

| File | Purpose |
|------|---------|
| `tf_activity_analysis.py` | Core analysis: pseudobulk, activity inference, OLS testing |
| `tf_activity_visualize.py` | Heatmaps, volcanos, epigenetic panel, hub barplot |
| `tf_convergence_analysis.py` | Post-hoc hub TF identification and cross-validation |
| `aceTfActT.sh` | Launcher: submits one SLURM job per phenotype |
| `run_tf_activity.sh` | SLURM job script: runs analysis + visualization |
| `smoke_test.sh` | Quick smoke test with fixture data |
