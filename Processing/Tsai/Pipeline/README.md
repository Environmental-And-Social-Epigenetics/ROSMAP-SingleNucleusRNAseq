# Tsai Pipeline

This directory contains a self-contained three-stage processing pipeline for the Tsai ROSMAP snRNA-seq dataset, starting from CellBender-filtered matrices and ending with an integrated, cell type annotated `AnnData` object.

## Inputs

- CellBender outputs: `Tsai_Data/Cellbender_Outputs/{projid}/processed_feature_bc_matrix_filtered.h5`
- Sample metadata: `Preprocessing/Tsai/02_Cellranger_Counts/Tracking/patient_metadata.csv`
- Marker reference copied into this repo: `Processing/Tsai/Pipeline/Resources/Brain_Human_PFC_Markers_Mohammadi2020.rds`

## Outputs

The pipeline writes data outside the repo into `Tsai_Data/Processing_Outputs/`:

```text
Tsai_Data/Processing_Outputs/
  01_QC_Filtered/         per-sample QC-filtered h5ad files + qc_summary.csv
  02_Doublet_Removed/     per-sample singlet-only h5ad files + doublet_summary.csv
  03_Integrated/          integrated and annotated h5ad files
  Logs/                   SLURM stdout/stderr
```

`qc_summary.csv` and `doublet_summary.csv` are automatically generated and
track per-sample cell counts before/after filtering across all array tasks.

## Pipeline Stages

### Stage 1: QC Filtering

- Script: `01_qc_filter.py`
- Wrapper: `01_qc_filter.sh`
- Environment: `QC_ENV` (see `config/paths.sh`; env spec in `envs/stage1_qc.yml`)
- Input: CellBender `processed_feature_bc_matrix_filtered.h5`
- Output: `{projid}_qc.h5ad`

Filtering rules (percentile-based):

| Metric | Threshold |
|--------|-----------|
| `log1p_total_counts` | Below 4.5th or above 96th percentile |
| `log1p_n_genes_by_counts` | Below 5th percentile |
| `pct_counts_mt` | Above 10% |

### Stage 2: Doublet Removal

- Script: `02_doublet_removal.Rscript`
- Wrapper: `02_doublet_removal.sh`
- Environment: `SINGLECELL_ENV` (see `config/paths.sh`; env spec in `envs/stage2_doublets.yml`)
- Method: `scDblFinder` (seed=123, 4 workers)
- Input: `{projid}_qc.h5ad`
- Output: `{projid}_singlets.h5ad`

If the active R environment is missing a required Bioconductor package such as `GenomeInfoDbData`, the script bootstraps it into `Processing/Tsai/Pipeline/R_libs/` on first run.

### Stage 3: Integration and Annotation

- Script: `03_integration_annotation.py`
- Wrapper: `03_integration_annotation.sh`
- Environment: `BATCHCORR_ENV` (see `config/paths.sh`; env spec in `envs/stage3_integration.yml`)
- Input: all Stage 2 singlet files
- Outputs:
  - `tsai_integrated.h5ad`
  - `tsai_annotated.h5ad`
  - `cluster_annotation_rankings.csv`
  - `cluster_annotation_top3.csv`
  - UMAP figures in `03_Integrated/figures/`

Stage 3 processing steps:

| Step | Function | Key parameters |
|------|----------|----------------|
| Normalization | `sc.pp.normalize_total` | Per-cell median-based normalization |
| Log transform | `sc.pp.log1p` | Natural log(1 + x) |
| HVG selection | `sc.pp.highly_variable_genes` | `flavor="seurat_v3"`, `n_top_genes=3000`, `layer="counts"` |
| PCA | `sc.tl.pca` | `n_comps=30`, `svd_solver="arpack"`, `use_highly_variable=True` |
| Batch correction | `harmonypy.run_harmony` | Input: `X_pca`; batch variable: `derived_batch` (flowcell-based, ~41 groups); configurable via `--harmony-batch-key`; `--harmony-theta 2.0`; `--skip-harmony` for no correction |
| Neighbors | `sc.pp.neighbors` | `n_neighbors=30`, `n_pcs=30`, `metric="cosine"`, `use_rep="X_harmony"` (or `"X_pca"` if `--skip-harmony`) |
| Clustering | `sc.tl.leiden` | Resolutions: **0.2**, **0.5**, **1.0** |
| UMAP | `sc.tl.umap` | `min_dist=0.15`, `random_state=0` |
| Annotation | `decoupler.run_ora` | Markers: Mohammadi 2020 PFC reference; `use_raw=False`; top-1 cell type per `leiden_res0_5` cluster |

## Quick Start

Run the entire pipeline with automatic dependency chaining:

```bash
cd <REPO_ROOT>/Processing/Tsai/Pipeline
./submit_pipeline.sh all
```

Submit individual stages:

```bash
./submit_pipeline.sh 1        # Stage 1 only
./submit_pipeline.sh 2 3      # Stages 2 and 3 (chained)
```

For testing on a small subset (runs locally, no SLURM):

```bash
python 01_qc_filter.py --sample-ids 10100574,10100862
Rscript 02_doublet_removal.Rscript --sample-ids 10100574,10100862
python 03_integration_annotation.py --sample-ids 10100574,10100862
```

## Sample Discovery

The per-sample stages auto-discover complete inputs and preserve metadata order from `patient_metadata.csv`.

The SLURM wrappers default to `#SBATCH --array=1-478%32` for the current
CellBender-complete dataset. `patient_metadata.csv` lists 480 Tsai samples; the
currently missing CellBender outputs are `11467746` and `20834164`. If that
count changes, override the array range at submit time:

```bash
sbatch --array=1-$(python 01_qc_filter.py --list-samples | wc -l) 01_qc_filter.sh
```

## Configuration

All shell wrappers source `config/paths.sh` at the repo root, which is the
single source of truth for conda environments, data directories, and SLURM log
paths.  To adapt the pipeline for a new cluster, edit `config/paths.sh` only —
no changes to the pipeline scripts are needed.

Conda environment specs are in `envs/` for recreating environments from scratch.
