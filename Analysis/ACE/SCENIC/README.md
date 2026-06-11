# ACE SCENIC -- Regulatory Network Inference

pySCENIC regulatory network analysis for the ACE (Adverse Childhood
Experiences) phenotype. Infers gene regulatory networks and transcription
factor regulons, then tests their association with ACE exposure.

## Method

1. Load integrated annotated h5ad, subset by cell type and sex
2. Micropool cells (~50 per pool) to reduce sparsity
3. CPM normalize micropools
4. Run GRNBoost2 for gene regulatory network inference
5. Prune with motif enrichment (cisTarget ranking databases)
6. Extract regulons (TF + target gene sets)
7. Score regulon activity per sample with AUCell
8. Test association between regulon activity and ACE phenotype using OLS
   regression with covariates (age_death, pmi, niareagansc)
9. FDR correction (Benjamini-Hochberg)

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

# Tsai dataset
bash Analysis/ACE/SCENIC/Tsai/aceScenicT.sh [INTEGRATION] [PHENOTYPE]
```

Defaults: INTEGRATION=derived_batch, PHENOTYPE=tot_adverse_exp

## Phenotype

Default phenotype is `tot_adverse_exp` (total adverse experiences count).
Other options: `early_hh_ses`, `ace_aggregate`. The phenotype column must
exist in `${ACE_SCORES_CSV}`.

Covariates: `age_death`, `pmi`, `niareagansc`.

Sex stratification: msex=1 (Male), msex=0 (Female). Patient lists are derived
from the phenotype metadata at runtime, not hardcoded.

## Subdirectories

- `Tsai/` -- SCENIC analysis on Tsai dataset
- `DeJager/` -- SCENIC analysis on DeJager dataset

## Outputs

Under `${ACE_OUTPUT_ROOT}/SCENIC/{cohort}/results_{integration}/{phenotype}/`:

```
{Sex}_{CellType}/
  adjacencies.csv.gz       # GRNBoost2 TF-gene edges
  regulons.pkl             # Regulon objects
  regulons_list.csv        # Regulon names and target counts
  auc_matrix.csv           # Micropool x regulon AUCell scores
  regression_results.csv   # Statistical results with FDR
```

## Resource Requirements

- Memory: 400 GB per job
- CPUs: 20 cores per job
- Time: up to 48 hours per cell type
