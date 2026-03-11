# DeJager Pipeline

This directory contains a self-contained three-stage processing pipeline for the DeJager ROSMAP snRNA-seq dataset, starting from CellBender-filtered matrices and ending with an integrated, cell type annotated `AnnData` object.

The DeJager pipeline uses **identical processing methods** to the Tsai pipeline, with the only data-specific differences being:

1. **Patient ID assignment** from a barcode-to-patient mapping CSV (DeJager libraries are multiplexed)
2. **Harmony batch correction on `patient_id`** (not `projid` or `batch`)

## Inputs

- CellBender outputs: `{DEJAGER_PREPROCESSED}/{library}/processed_feature_bc_matrix_filtered.h5`
- Patient barcode mapping: `/om/scratch/Mon/shared_folder/WGS/cell_to_patient_assignmentsFinal0.csv`
- Patient ID overrides for "alone" libraries: `Resources/patient_id_overrides.json`
- Marker reference: `Resources/Brain_Human_PFC_Markers_Mohammadi2020.rds` (symlink to Tsai's copy)

## Outputs

The pipeline writes data into `DeJager_Data/Processing_Outputs/`:

```text
DeJager_Data/Processing_Outputs/
  01_QC_Filtered/         per-library QC-filtered h5ad files + qc_summary.csv
  02_Doublet_Removed/     per-library singlet-only h5ad files + doublet_summary.csv
  03_Integrated/          integrated and annotated h5ad files
  Logs/                   SLURM stdout/stderr
```

## Pipeline Stages

### Stage 1: QC Filtering

- Script: `01_qc_filter.py`
- Wrapper: `01_qc_filter.sh`
- Environment: `QC_ENV` (env spec in `envs/stage1_qc.yml`)
- Input: CellBender `processed_feature_bc_matrix_filtered.h5`
- Output: `{library}_qc.h5ad`

**Patient ID assignment** (DeJager-specific):
- For multiplexed libraries: cell barcodes are mapped to patient IDs via `cell_to_patient_assignmentsFinal0.csv`
- For "alone" libraries (suffix `-alone`): the R-number is extracted from the library name and mapped to a patient ID via `patient_id_overrides.json`
- Cells without a patient ID mapping are dropped

**QC filtering** (percentile-based, identical to Tsai):

| Metric | Threshold |
|--------|-----------|
| `log1p_total_counts` | Below 4.5th or above 96th percentile |
| `log1p_n_genes_by_counts` | Below 5th percentile |
| `pct_counts_mt` | Above 10% |

### Stage 2: Doublet Removal

- Script: `02_doublet_removal.Rscript`
- Wrapper: `02_doublet_removal.sh`
- Environment: `SINGLECELL_ENV` (env spec in `envs/stage2_doublets.yml`)
- Method: `scDblFinder` (seed=123, 4 workers)
- Input: `{library}_qc.h5ad`
- Output: `{library}_singlets.h5ad`

If the active R environment is missing a required Bioconductor package such as `GenomeInfoDbData`, the script bootstraps it into `Processing/DeJager/Pipeline/R_libs/` on first run.

### Stage 3: Integration and Annotation

- Script: `03_integration_annotation.py`
- Wrapper: `03_integration_annotation.sh`
- Environment: `BATCHCORR_ENV` (env spec in `envs/stage3_integration.yml`)
- Input: all Stage 2 singlet files
- Outputs:
  - `dejager_integrated.h5ad`
  - `dejager_annotated.h5ad`
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
| Batch correction | `harmonypy.run_harmony` | Input: `X_pca`; batch variable: **`patient_id`** |
| Neighbors | `sc.pp.neighbors` | `n_neighbors=30`, `n_pcs=30`, `metric="cosine"`, `use_rep="X_harmony"` |
| Clustering | `sc.tl.leiden` | Resolutions: **0.2**, **0.5**, **1.0** |
| UMAP | `sc.tl.umap` | `min_dist=0.15`, `random_state=0` |
| Annotation | `decoupler.run_ora` | Markers: Mohammadi 2020 PFC reference; `use_raw=False`; top-1 cell type per `leiden_res0_5` cluster |

## Quick Start

Run the entire pipeline with automatic dependency chaining:

```bash
cd <REPO_ROOT>/Processing/DeJager/Pipeline
./submit_pipeline.sh all
```

Submit individual stages:

```bash
./submit_pipeline.sh 1        # Stage 1 only
./submit_pipeline.sh 2 3      # Stages 2 and 3 (chained)
```

For testing on a small subset (runs locally, no SLURM):

```bash
python 01_qc_filter.py --sample-ids LIB5001_R1234567-alone
Rscript 02_doublet_removal.Rscript --sample-ids LIB5001_R1234567-alone
python 03_integration_annotation.py --sample-ids LIB5001_R1234567-alone
```

## Sample Discovery

Stage 1 discovers libraries by scanning `--input-dir` for subdirectories containing `processed_feature_bc_matrix_filtered.h5`. Stages 2-3 discover samples via `*_qc.h5ad` / `*_singlets.h5ad` file patterns.

The SLURM wrappers default to `#SBATCH --array=1-200%32`. Array tasks beyond the actual sample count exit gracefully. To override:

```bash
sbatch --array=1-$(python 01_qc_filter.py --list-samples | wc -l) 01_qc_filter.sh
```

## Configuration

All shell wrappers source `config/paths.sh` at the repo root, which is the
single source of truth for conda environments, data directories, and SLURM log
paths. Key DeJager-specific variables:

| Variable | Purpose |
|----------|---------|
| `DEJAGER_PREPROCESSED` | CellBender output directory |
| `DEJAGER_QC_FILTERED` | Stage 1 output |
| `DEJAGER_DOUBLET_REMOVED` | Stage 2 output |
| `DEJAGER_INTEGRATED` | Stage 3 output |
| `DEJAGER_PATIENT_MAP` | Barcode-to-patient CSV |
| `DEJAGER_PATIENT_ID_OVERRIDES` | R-number overrides JSON |

Conda environment specs are in `envs/` for recreating environments from scratch.

## Resources

- `patient_id_overrides.json` — Maps R-numbers to patient IDs for "alone" libraries where the barcode mapping CSV does not apply.
- `Brain_Human_PFC_Markers_Mohammadi2020.rds` — Symlink to the shared marker reference in `Processing/Tsai/Pipeline/Resources/`.
