# Phenotype DEG Comparison — DeJager

Compares 5 batch correction strategies (`library_id`, `patient_id`, `pool_batch`, `derived_batch`, `sequencing_date`) by running pseudobulk DESeq2 for two well-known phenotypes (sex, AD-vs-control) and measuring effect sizes on curated truth-set vs negative-control gene lists.

The batch correction that **amplifies known biology** while **suppressing noise** has the highest signal-to-noise ratio (SNR).

## Inputs

- `Processing_Outputs/03_Integrated*/dejager_annotated.h5ad` — produced by Stage 3 / 3b / 3c / 3d / 3e
- `Data/Phenotypes/TSAI_DEJAGER_all_patients_wACEscores.csv` — has `msex`, `cogdx`, `age_death`, `pmi`

## Pipeline

```
Phase D1: prep_splits.sh  per integration  (5 jobs, ~2-4h each)
   ↓ produces celltype_splits_<integration>/{Exc,Inh,Ast,Mic,Oli,OPC}.h5ad
Phase D2: phenoDeg.sh     per (integration, phenotype, celltype) (60 jobs, ~30 min each)
   ↓ produces results_<integration>/<phenotype>/deseq_<phenotype>_<celltype>.csv
Phase D3: aggregate_effects.py     (1 job, ~10 min)
   ↓ produces comparison/{effect_size_summary.csv, effect_size_ranking.csv, *.png}
```

## Run

```bash
bash run_pipeline.sh
```

## Truth Sets (in `truth_genes.py`)

- **Sex (msex):** XIST, RPS4Y1, DDX3Y, KDM5D, UTY, ... (Y-chromosome + X-inactivation)
- **AD (cogdx_binary):** APOE, BIN1, MS4A6A, TREM2, CR1, ... (top GWAS hits + canonical markers)
- **Negative control:** ACTB, GAPDH, B2M, RPL13A, ... (housekeeping, should show zero effect)

## Interpretation

- **High `truth_mean_abs_lfc`** = batch correction preserves real biology
- **Low `neg_ctrl_mean_abs_lfc`** = batch correction doesn't inflate noise
- **High SNR = truth/noise** = best batch correction

The `truth_vs_noise_scatter.png` plot makes the trade-off visible: points well above the diagonal indicate real biological signal.
