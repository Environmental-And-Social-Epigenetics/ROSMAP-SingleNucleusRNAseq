# SCENIC — ACE — DeJager

SCENIC gene regulatory network inference for the ACE (Adverse Childhood
Experiences) phenotype on the DeJager cohort.

DeJager runs the **same modular pipeline as Tsai** — there is no separate
DeJager analysis code. The launchers here delegate to the canonical Tsai
modular Python (`../Tsai/scenic_analysis.py`) with `--cohort dejager`, so
micropooling, CPM, GRNBoost2, cisTarget pruning, AUCell, and the OLS
association are **bit-for-bit the same method**.

## Identical method, single pool_size

- The micropool `pool_size` is a single CLI arg (`--pool-size`, default 50),
  used for **both cohorts and both sexes**. This replaces the deprecated
  monolith's inconsistent 30 (male) / 39 (female).
- `--cohort dejager` only changes which `obs` column is used as the
  cell->patient key. Both cohorts key cells by ROSMAP `projid` (the DeJager DEG
  prep step sets `obs['projid'] = patient_id`, already a ROSMAP projid) and
  merge against the same pooled `${ACE_SCORES_CSV}`
  (`TSAI_DEJAGER_all_patients_wACEscores.csv`). The DeJager ID map
  (`${DEJAGER_ID_MAP}`, projid <-> individualID) is therefore not needed by the
  SCENIC step — the join key is already a projid.

## Entry point

```bash
source config/paths.sh
bash Analysis/ACE/SCENIC/DeJager/aceScenicDJ.sh [INTEGRATION] [PHENOTYPE]

# Defaults: INTEGRATION=library_id, PHENOTYPE=tot_adverse_exp
```

This submits one independent SLURM job per cell-type x sex via
`DeJager/run_scenic.sh`, which forwards to `../Tsai/run_scenic.sh ... dejager`.

## Files

| File | Description |
|---|---|
| `aceScenicDJ.sh` | Launcher — submits one SLURM job per cell-type x sex (mirrors `../Tsai/aceScenicT.sh`, DeJager paths). |
| `run_scenic.sh` | Thin wrapper — forwards to `../Tsai/run_scenic.sh` with `cohort=dejager`. No analysis logic of its own. |
| `legacy/aceScenic.py` | Archived monolithic DeJager script (DEPRECATED; inconsistent pool_size 30/39). |
| `legacy/aceScenic.sh` | Archived SLURM wrapper (DEPRECATED). |

## Inputs / Outputs

```
inputs : ${ACE_OUTPUT_ROOT}/DEG/DeJager/celltype_splits_${INTEGRATION}/${CELL_TYPE}.h5ad
phenos : ${ACE_SCORES_CSV}   (pooled Tsai+DeJager, keyed by projid)
outputs: ${ACE_OUTPUT_ROOT}/SCENIC/DeJager/results_${INTEGRATION}/${PHENOTYPE}/${SEX}_${CELL_TYPE}/
           adjacencies.csv.gz
           regulons.pkl
           regulons_list.csv
           auc_matrix.csv
           regression_results.csv
logs   : ${ACE_OUTPUT_ROOT}/SCENIC/DeJager/logs/
```

## Visualization

The Tsai `scenic_visualize.py` is cohort-agnostic; point it at the DeJager
results dir:

```bash
source config/paths.sh
activate_env "${SCENIC_ANALYSIS_ENV}"

python Analysis/ACE/SCENIC/Tsai/scenic_visualize.py \
  --results-dir "${ACE_OUTPUT_ROOT}/SCENIC/DeJager/results_library_id/tot_adverse_exp" \
  --output-dir  "${ACE_OUTPUT_ROOT}/SCENIC/DeJager/figures" \
  --phenotype tot_adverse_exp
```

## Prerequisites

Same as Tsai (see `../Tsai/README.md`): `${SCENIC_ANALYSIS_ENV}` conda
environment, the pySCENIC ranking databases in `${SCENIC_RANKING_DIR}`, and the
DeJager DEG celltype splits produced by `Analysis/ACE/DEG/DeJager/aceDegDJ.sh`.
