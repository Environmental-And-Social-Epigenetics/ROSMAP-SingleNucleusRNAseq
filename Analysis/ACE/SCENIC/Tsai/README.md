# SCENIC -- ACE -- Tsai

pySCENIC gene regulatory network inference for the ACE (Adverse Childhood
Experiences) phenotype on the Tsai dataset.

## Method

1. **Load** cell-type-specific h5ad (from DEG preprocessing splits)
2. **Filter by sex** using metadata (msex: 1=Male, 0=Female) -- derived from
   phenotype CSV, not hardcoded patient lists
3. **Micropool** cells (~50 per pool per patient) to reduce single-cell
   sparsity while preserving patient-level variability
4. **CPM normalize** micropools (counts per million)
5. **GRNBoost2** -- infer gene regulatory network (TF-gene adjacencies)
6. **cisTarget motif pruning** -- filter adjacencies by cis-regulatory motif
   enrichment using two ranking databases (500bp upstream/100bp downstream and
   10kbp flanking)
7. **Regulon extraction** -- convert pruned motif enrichment results to
   regulons (TF + target gene sets)
8. **AUCell** -- score regulon activity per micropool
9. **OLS regression** -- test association between regulon activity and ACE
   phenotype (`tot_adverse_exp` by default) controlling for `age_death`,
   `pmi`, and `niareagansc`
10. **FDR correction** (Benjamini-Hochberg) across regulons

## Prerequisites

- `${SCENIC_ANALYSIS_ENV}` conda environment with pySCENIC, arboreto, scanpy
- pySCENIC ranking databases in `${SCENIC_RANKING_DIR}`:
  - `hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather`
  - `hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather`
  - `motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl`
  - `hg.txt` (TF names list)
- Download databases from: https://resources.aertslab.org/cistarget/
- Cell-type h5ad splits produced by the DEG preprocessing step

## Entry Point

```bash
# From the repo root:
source config/paths.sh
bash Analysis/ACE/SCENIC/Tsai/aceScenicT.sh [INTEGRATION] [PHENOTYPE]

# Defaults: INTEGRATION=derived_batch, PHENOTYPE=tot_adverse_exp
```

This submits one SLURM job per cell-type x sex combination via `run_scenic.sh`.

## Files

| File | Description |
|---|---|
| `aceScenicT.sh` | Main launcher -- submits SLURM jobs |
| `run_scenic.sh` | SLURM batch script for one cell-type x sex. **Shared by both cohorts**: takes an optional 5th `COHORT` arg (`tsai`\|`dejager`, default `tsai`) that selects the cohort-specific DEG input dir and SCENIC output dir. |
| `scenic_analysis.py` | Core pySCENIC pipeline (GRNBoost2 + cisTarget + AUCell + OLS). Cohort-aware via `--cohort {tsai,dejager}`; `--pool-size` (default 50) is a single arg shared by both cohorts. |
| `scenic_associate.py` | Per-arm OLS over cached AUCell (cohort-agnostic; `--cohort` for provenance). |
| `scenic_visualize.py` | Post-hoc visualization (volcano, heatmap, bar plots); cohort-agnostic, runs on either cohort's results dir. |
| `smoke_test.sh` | Quick validation with synthetic data |
| `legacy/aceScenic.py` | Archived monolithic script (DEPRECATED; kept for reference) |
| `legacy/aceScenic.sh` | Archived SLURM wrapper (DEPRECATED; kept for reference) |

### Cohort parity

The DeJager cohort runs **this same pipeline**. `Analysis/ACE/SCENIC/DeJager/`
contains only thin launchers (`aceScenicDJ.sh`, `run_scenic.sh`) that delegate
to the files here with `--cohort dejager`. There is no separate DeJager analysis
code, and the micropool `pool_size` (default 50) is **identical** for both
cohorts and both sexes — see `../README.md` "Unified pipeline".

## Output Structure

```
${ACE_OUTPUT_ROOT}/SCENIC/Tsai/
  results_${INTEGRATION}/${PHENOTYPE}/
    Male_Mic/
      adjacencies.csv.gz       # GRNBoost2 TF-gene edges (compressed)
      regulons.pkl             # Regulon objects (pickle)
      regulons_list.csv        # Human-readable regulon names + target count
      auc_matrix.csv           # Micropool x regulon AUCell scores
      regression_results.csv   # Regulon, coef, stderr, pvalue, padj, n_targets
    Female_Ast/
      ...
  logs/
    %j_Male_Mic.out
    %j_Male_Mic.err
    ...
```

## Visualization

After all jobs complete, generate summary figures:

```bash
source config/paths.sh
activate_env "${SCENIC_ANALYSIS_ENV}"

python Analysis/ACE/SCENIC/Tsai/scenic_visualize.py \
  --results-dir "${ACE_OUTPUT_ROOT}/SCENIC/Tsai/results_derived_batch/tot_adverse_exp" \
  --output-dir "${ACE_OUTPUT_ROOT}/SCENIC/Tsai/figures" \
  --phenotype tot_adverse_exp
```

Produces:
- `scenic_summary_tot_adverse_exp.pdf` -- volcano plots, heatmap, bar chart
- `all_regression_results_tot_adverse_exp.csv` -- combined results table

## Interpreting Results

- **regression_results.csv**: Each row is a regulon (TF). A positive
  coefficient means higher regulon activity is associated with higher ACE
  exposure. `padj < 0.05` indicates significance after FDR correction.
- **Volcano plot**: x-axis = effect size (coefficient), y-axis = -log10(FDR).
  Regulons above the dashed line are significant.
- **Heatmap**: Shows which regulons are significant across cell types,
  revealing shared vs cell-type-specific regulatory programs.

## Resource Requirements

- **Memory**: 400 GB per job (GRNBoost2 + ranking database loading)
- **CPUs**: 20 cores per job
- **Time**: Up to 48 hours per cell type (varies by cell count)
- **Disk**: ~2-5 GB per cell type (adjacencies are the largest output)

## Smoke Test

```bash
bash Analysis/ACE/SCENIC/Tsai/smoke_test.sh
```

Generates synthetic data and runs the pipeline in `--smoke` mode (500 genes,
skip motif pruning) to verify the code executes without errors.
