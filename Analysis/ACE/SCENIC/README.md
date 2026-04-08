# ACE SCENIC — Regulatory Network Inference

pySCENIC regulatory network analysis for the ACE phenotype, adapted from the
SocIsl implementation in `Analysis/SocIsl/SCENIC/`.

## Method

1. Load integrated annotated H5AD, subset by cell type and sex
2. Micropool cells (~50 per pool) to reduce sparsity
3. Run GRNBoost2 for gene regulatory network inference
4. Prune with motif enrichment (cisTarget ranking databases)
5. Score regulon activity per sample with AUCell
6. Test association between regulon activity and ACE phenotype using OLS
   regression with covariates (age, PMI, NIA, batch)
7. FDR correction (Benjamini-Hochberg)

## Prerequisites

- `${SCENIC_ANALYSIS_ENV}` conda environment
- pySCENIC ranking databases in `${SCENIC_RANKING_DIR}`:
  - `hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather`
  - `hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather`
  - `motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl`
  - `hg.txt` (TF names list)
- Download from: https://resources.aertslab.org/cistarget/

## Running

```bash
source config/paths.sh
cd Analysis/ACE/SCENIC/Tsai
sbatch aceScenic.sh
```

## Outputs

Under `${ACE_OUTPUT_ROOT}/SCENIC/Tsai/`:
- `{sex}_regulonsFULL_Tsai{celltype}.csv` — regulon list
- `{sex}_auc_mtxFULL_Tsai{celltype}.csv` — AUCell activity matrix
- `{sex}_pVals{celltype}.csv` — statistical results with FDR
- `{celltype}{sex}Volcano.png` — volcano plots

## Phenotype

Replaces `social_isolation_avg` with ACE phenotype column (configurable via
`--phenotype`, default: `ace_aggregate`). Covariates: `age_death`, `pmi`,
`niareagansc`, `patient_id` (batch).

## Subdirectories

- `DeJager/` — SCENIC analysis on DeJager dataset
- `Tsai/` — SCENIC analysis on Tsai dataset
