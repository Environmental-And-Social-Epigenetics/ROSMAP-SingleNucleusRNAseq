# Processing

This directory contains scripts for quality control, normalization, batch correction, and cell type annotation of snRNA-seq data.

## Overview

The Processing phase takes CellBender-corrected count matrices and prepares them for downstream analysis through:

1. Quality control (percentile-based outlier filtering)
2. Doublet removal (scDblFinder)
3. Normalization and feature selection (seurat_v3 HVGs)
4. Batch correction (Harmony)
5. Dimensionality reduction and clustering (PCA, UMAP, Leiden)
6. Cell type annotation (ORA with Mohammadi 2020 markers)

## Directory Structure

```
Processing/
├── DeJager/
│   ├── Pipeline/                    # Primary 3-stage pipeline (SLURM-automated)
│   │   ├── 01_qc_filter.*          # Stage 1: QC filtering + patient ID assignment
│   │   ├── 02_doublet_removal.*    # Stage 2: scDblFinder doublet removal
│   │   ├── 03_integration_annotation.*  # Stage 3: Harmony + annotation
│   │   ├── submit_pipeline.sh      # SLURM orchestrator
│   │   ├── Resources/              # Patient ID overrides, marker references
│   │   └── envs/                   # Conda environment specs
│   └── _legacy/                    # Archived original scripts
└── Tsai/
    ├── Pipeline/                    # Primary 3-stage pipeline (SLURM-automated)
    │   ├── 01_qc_filter.*          # Stage 1: QC filtering
    │   ├── 02_doublet_removal.*    # Stage 2: scDblFinder doublet removal
    │   ├── 03_integration_annotation.*  # Stage 3: Harmony + annotation
    │   ├── submit_pipeline.sh      # SLURM orchestrator
    │   ├── Resources/              # Marker gene references
    │   └── envs/                   # Conda environment specs
    ├── QC/                          # Legacy per-sample scripts
    ├── Batch_Correction/            # Legacy batch correction scripts
    └── *.ipynb                      # Prototyping notebooks
```

## Unified Processing Methods

Both pipelines use **identical processing methods**. The only data-specific differences are:

| Parameter | Tsai | DeJager |
|-----------|------|---------|
| Harmony batch variable | `projid` | `patient_id` |
| Patient ID source | From `patient_metadata.csv` | From barcode mapping CSV + overrides JSON |
| Output file prefix | `tsai_` | `dejager_` |

### Shared Parameters

#### Stage 1: QC Filtering

| Metric | Threshold |
|--------|-----------|
| `log1p_total_counts` | Below 4.5th or above 96th percentile |
| `log1p_n_genes_by_counts` | Below 5th percentile |
| `pct_counts_mt` | Above 10% |

#### Stage 2: Doublet Removal

| Parameter | Value |
|-----------|-------|
| Method | `scDblFinder` |
| Random seed | 123 |
| Parallelism | MulticoreParam, 4 workers |
| Filter | Retain `scDblFinder.class == "singlet"` |

#### Stage 3: Integration & Annotation

| Step | Function | Key Parameters |
|------|----------|----------------|
| Normalization | `sc.pp.normalize_total` + `sc.pp.log1p` | Median-based, no scaling |
| HVG selection | `sc.pp.highly_variable_genes` | `flavor="seurat_v3"`, `n_top_genes=3000`, `layer="counts"` |
| PCA | `sc.tl.pca` | `n_comps=30`, `svd_solver="arpack"` |
| Batch correction | `harmonypy.run_harmony` | Input: `X_pca` |
| Neighbors | `sc.pp.neighbors` | `n_neighbors=30`, `n_pcs=30`, `metric="cosine"` |
| Clustering | `sc.tl.leiden` | Resolutions: 0.2, 0.5, 1.0 |
| UMAP | `sc.tl.umap` | `min_dist=0.15`, `random_state=0` |
| Annotation | `decoupler.run_ora` | Mohammadi 2020 PFC markers, `use_raw=False` |

### Resource Requirements

| Stage | Cores | Memory | Time | Job Type |
|-------|-------|--------|------|----------|
| 1 | 4 | 32GB | ~1h/sample | SLURM array |
| 2 | 4 | 32GB | ~1h/sample | SLURM array |
| 3 | 32 | 500GB | 4-8h | Single job |

## Quick Start

```bash
# Run the Tsai pipeline
cd <REPO_ROOT>/Processing/Tsai/Pipeline
./submit_pipeline.sh all

# Run the DeJager pipeline
cd <REPO_ROOT>/Processing/DeJager/Pipeline
./submit_pipeline.sh all
```

See each pipeline's README for detailed instructions:
- [Tsai/Pipeline/README.md](Tsai/Pipeline/README.md)
- [DeJager/Pipeline/README.md](DeJager/Pipeline/README.md)

## Next Steps

After processing, proceed to the `Analysis/` directory for differential
expression and downstream analyses.
