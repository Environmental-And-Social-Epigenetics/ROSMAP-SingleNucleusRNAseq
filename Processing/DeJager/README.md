# DeJager Processing

This directory contains the processing pipeline for the DeJager ROSMAP snRNA-seq dataset.

## Pipeline

The active pipeline is in `Pipeline/`. See [Pipeline/README.md](Pipeline/README.md) for full documentation.

The pipeline consists of three SLURM-automated stages:

1. **QC Filtering** — Percentile-based outlier removal + patient ID assignment from barcode mapping
2. **Doublet Removal** — scDblFinder singlet detection
3. **Integration & Annotation** — Harmony batch correction on patient_id, Leiden clustering, ORA cell type annotation

All processing parameters are identical to the Tsai pipeline (see `Processing/Tsai/Pipeline/`).

## Legacy Scripts

Original scripts are archived in `_legacy/` for reference. These are superseded by the `Pipeline/` implementation.
