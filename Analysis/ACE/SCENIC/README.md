# ACE SCENIC -- Regulatory Network Inference

pySCENIC regulatory network analysis for the ACE (Adverse Childhood
Experiences) phenotype. Infers gene regulatory networks and transcription
factor regulons, then tests their association with ACE exposure.

## Unified pipeline (both cohorts)

Tsai and DeJager run the **same modular pipeline**
(`Tsai/scenic_analysis.py`). The only cohort difference is the input DEG
celltype-split directory and the output directory; the METHOD is identical.

- **Single `pool_size` for both cohorts and both sexes** — the micropool size
  is one CLI arg (`--pool-size`, default 50) shared by Tsai and DeJager. The
  old DeJager monolith used an inconsistent 30 (male) / 39 (female); that is
  gone (the script is archived under `DeJager/legacy/`).
- **Cohort flag** — `scenic_analysis.py --cohort {tsai,dejager}` only selects
  which obs column is used as the cell->patient key. Both cohorts ultimately
  key cells by ROSMAP `projid` and merge against the same pooled phenotype CSV,
  so micropooling, CPM, GRNBoost2, cisTarget, AUCell and the OLS association are
  bit-for-bit the same code.
- **Single source of truth** — `DeJager/run_scenic.sh` delegates to
  `Tsai/run_scenic.sh` with `dejager`; there is no duplicated DeJager analysis
  code.

## Method

1. Load integrated annotated h5ad, subset by cell type and sex
2. Micropool cells (`--pool-size`, default 50) per patient to reduce sparsity
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
# Defaults: INTEGRATION=derived_batch, PHENOTYPE=tot_adverse_exp

# DeJager dataset (same pipeline, --cohort dejager, same pool_size)
bash Analysis/ACE/SCENIC/DeJager/aceScenicDJ.sh [INTEGRATION] [PHENOTYPE]
# Defaults: INTEGRATION=library_id, PHENOTYPE=tot_adverse_exp
```

## Phenotype

Default phenotype is `tot_adverse_exp` (total adverse experiences count).
Other options: `early_hh_ses`, `ace_aggregate`. The phenotype column must
exist in `${ACE_SCORES_CSV}`.

Covariates: `age_death`, `pmi`, `niareagansc`.

Sex stratification: msex=1 (Male), msex=0 (Female). Patient lists are derived
from the phenotype metadata at runtime, not hardcoded.

## Subdirectories

- `Tsai/` -- modular SCENIC pipeline (canonical method: `scenic_analysis.py`,
  `scenic_associate.py`, `scenic_visualize.py`, `run_scenic.sh`, `aceScenicT.sh`)
- `DeJager/` -- DeJager launchers (`aceScenicDJ.sh`, `run_scenic.sh`) that reuse
  the Tsai modular pipeline with `--cohort dejager`
- `{Tsai,DeJager}/legacy/` -- archived monolithic `aceScenic.py`/`aceScenic.sh`
  (deprecated; do not run)

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
