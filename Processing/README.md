# Processing

This directory contains scripts for quality control, normalization, batch correction, and cell type annotation of snRNA-seq data.

## Overview

The Processing phase takes CellBender-corrected count matrices and prepares them for downstream analysis through:

1. Quality control (doublet removal, outlier filtering)
2. Normalization and feature selection
3. Batch correction (Harmony)
4. Dimensionality reduction and clustering
5. Cell type annotation

## Directory Structure

```
Processing/
├── DeJager/
│   ├── firstStageTsaiPipeline.py    # Initial QC and processing
│   ├── secondStageTsaiPipeline.Rscript  # R-based processing
│   ├── thirdStageTsaiPipeline.py    # Final processing steps
│   ├── doublets.Rscript             # DoubletFinder
│   ├── batchCorrect.py              # Batch correction
│   └── onePatientScript.py          # Per-patient processing
└── Tsai/
    ├── Pipeline/                    # Primary 3-stage pipeline (current)
    │   ├── 01_qc_filter.*           # Stage 1: QC filtering
    │   ├── 02_doublet_removal.*     # Stage 2: scDblFinder doublet removal
    │   ├── 03_integration_annotation.*  # Stage 3: Harmony + annotation
    │   ├── Resources/               # Marker gene references
    │   └── envs/                    # Conda environment specs
    ├── QC/                          # Legacy per-sample scripts
    ├── Batch_Correction/            # Legacy batch correction scripts
    └── *.ipynb                      # Prototyping notebooks
```

## Tsai Pipeline (Primary Workflow)

The current Tsai processing pipeline is in `Tsai/Pipeline/`.  See
[Tsai/Pipeline/README.md](Tsai/Pipeline/README.md) for full documentation.

### Stages

| Stage | Script | Method |
|-------|--------|--------|
| 1 — QC Filtering | `01_qc_filter.py` | Percentile-based outlier removal, MT% filter |
| 2 — Doublet Removal | `02_doublet_removal.Rscript` | scDblFinder |
| 3 — Integration & Annotation | `03_integration_annotation.py` | Harmony, Leiden, ORA annotation |

### Resource Requirements

| Stage | Cores | Memory | Time |
|-------|-------|--------|------|
| 1 | 4 | 32GB | ~1h/sample |
| 2 | 4 | 32GB | ~1h/sample |
| 3 | 32 | 256GB | 4-8h |

### Outputs

| Stage | Output |
|-------|--------|
| 1 | `{projid}_qc.h5ad`, `qc_summary.csv` |
| 2 | `{projid}_singlets.h5ad`, `doublet_summary.csv` |
| 3 | `tsai_integrated.h5ad`, `tsai_annotated.h5ad`, annotation CSVs, UMAP figures |

## DeJager Pipeline

The DeJager scripts in `DeJager/` follow a similar workflow but are named
`*TsaiPipeline.*` because the pipeline was originally developed on Tsai data.

## Next Steps

After processing, proceed to the `Analysis/` directory for differential
expression and downstream analyses.

