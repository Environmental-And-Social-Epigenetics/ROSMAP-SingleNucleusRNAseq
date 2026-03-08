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

Filtering rules:

- Mitochondrial percentage > 10%
- `log1p_total_counts` outside the 4.5th to 96th percentile
- `log1p_n_genes_by_counts` outside the 5th to 100th percentile
- `pct_counts_in_top_20_genes` beyond 4 MADs from median (disable with `--top20-nmads 0`)

### Stage 2: Doublet Removal

- Script: `02_doublet_removal.Rscript`
- Wrapper: `02_doublet_removal.sh`
- Environment: `SINGLECELL_ENV` (see `config/paths.sh`; env spec in `envs/stage2_doublets.yml`)
- Method: `scDblFinder`
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

Stage 3 performs:

1. Metadata join on `projid`
2. Raw counts preservation in `layers["counts"]`
3. Normalization and log transform
4. HVG selection with `seurat_v3`
5. PCA
6. Harmony batch correction
7. Leiden clustering at 0.2, 0.5, and 1.0
8. UMAP
9. ORA-based annotation with the Mohammadi 2020 marker list

## Batch Variable

Harmony uses the metadata column `batch`, which is the sequencing batch identifier from `patient_metadata.csv`.

This is the correct technical covariate for correction because each Tsai sample corresponds to a single individual. Correcting on `projid` would remove biological signal rather than batch effects.

## Sample Discovery

The per-sample stages auto-discover complete inputs and preserve metadata order from `patient_metadata.csv`.

Current data snapshot:

- 474 sample directories in `Tsai_Data/Cellbender_Outputs/`
- 447 complete CellBender outputs with `processed_feature_bc_matrix_filtered.h5`

The SLURM wrappers default to `#SBATCH --array=1-447%32` for the current dataset. If that count changes, override the array range at submit time instead of editing the scripts:

```bash
cd <REPO_ROOT>/Processing/Tsai/Pipeline
sbatch --array=1-$(python 01_qc_filter.py --list-samples | wc -l) 01_qc_filter.sh
sbatch --array=1-$(Rscript 02_doublet_removal.Rscript --list-samples | wc -l) 02_doublet_removal.sh
```

## Suggested Run Order

```bash
cd <REPO_ROOT>/Processing/Tsai/Pipeline
sbatch 01_qc_filter.sh
sbatch 02_doublet_removal.sh
sbatch 03_integration_annotation.sh
```

For testing on a small subset:

```bash
python 01_qc_filter.py --sample-ids 10100574,10100862
Rscript 02_doublet_removal.Rscript --sample-ids 10100574,10100862
python 03_integration_annotation.py --sample-ids 10100574,10100862
```

## Configuration

All shell wrappers source `config/paths.sh` at the repo root, which is the
single source of truth for conda environments, data directories, and SLURM log
paths.  To adapt the pipeline for a new cluster, edit `config/paths.sh` only —
no changes to the pipeline scripts are needed.

Conda environment specs are in `envs/` for recreating environments from scratch.

## Notes

- The pipeline is self-contained within the repo for scripts and secondary resources.
- It reuses the tracked Tsai metadata file already present in the repository.
- The Stage 3 marker reference is now stored in `Resources/` so annotation does not depend on external paths.
- Stages 1 and 2 produce per-sample summary CSVs (`qc_summary.csv`, `doublet_summary.csv`) for tracking cell counts across the dataset.
