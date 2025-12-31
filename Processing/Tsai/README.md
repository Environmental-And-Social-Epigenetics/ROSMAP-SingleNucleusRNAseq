# Tsai Processing

Processing scripts for the Tsai ROSMAP snRNA-seq dataset.

## Overview

This directory contains the full processing pipeline for Tsai snRNA-seq data, including QC, normalization, batch correction, and cell type annotation. The scripts are organized into subdirectories for different processing stages.

## Directory Structure

```
Tsai/
├── QC/
│   ├── Doublets/           # Doublet detection
│   └── Outliers/           # Outlier removal
├── Batch_Correction/        # Harmony batch correction
├── Norm_FeatSel_DimRed_Cluster.ipynb  # Main processing notebook
└── Tsai_Sample_SingleCell_Pipeline.ipynb  # Sample pipeline
```

## QC

### Doublets

Uses DoubletFinder to identify and remove doublets.

**Scripts:**
- `Master_Remove_Doublets.py`: Generates per-sample scripts
- `Master_Remove_Doublets.r`: R implementation for all samples
- `example_doublet.r`: Example per-sample script

**Usage:**
```bash
# Generate per-sample scripts
python Master_Remove_Doublets.py

# Or run master R script
Rscript Master_Remove_Doublets.r
```

### Outliers

Removes cells with extreme QC metrics.

**Scripts:**
- `Master_Remove_Outliers.py`: Generates per-sample scripts
- `example_outlier.py`: Example per-sample script

**Filtering Criteria:**
- Minimum genes: 200
- Maximum mitochondrial %: 20%
- UMI count thresholds (sample-specific)

## Batch Correction

Uses Harmony to correct for batch effects across samples.

**Scripts:**

| Script | Description |
|--------|-------------|
| `Batch_Correction.py` | Main batch correction |
| `Batch_Correction2.py` | Alternative implementation |
| `Batch_Correction.sh` | SLURM wrapper |
| `featureSelOnward.py` | Feature selection and PCA |
| `featureSelOnward2.py` | Updated implementation |
| `batchCluster.py` | Clustering post-correction |
| `anno.py` | Cell type annotation |

**Usage:**
```bash
# Submit batch correction job
sbatch Batch_Correction.sh

# Or run directly
python Batch_Correction.py --input <h5ad> --output <output_dir>
```

## Main Notebooks

### Norm_FeatSel_DimRed_Cluster.ipynb

Complete processing pipeline notebook:
1. Load QC'd data
2. Normalization (scanpy `normalize_total`, `log1p`)
3. Feature selection (highly variable genes)
4. Dimensionality reduction (PCA, UMAP)
5. Clustering (Leiden)

### Tsai_Sample_SingleCell_Pipeline.ipynb

Sample/test pipeline for prototyping new analyses.

## Workflow

```
CellBender outputs
    ↓
Doublet Detection (DoubletFinder)
    ↓
Outlier Removal (QC metrics)
    ↓
Normalization + HVG Selection
    ↓
PCA
    ↓
Harmony Batch Correction
    ↓
UMAP + Clustering
    ↓
Cell Type Annotation
```

## Resource Requirements

| Step | Cores | Memory | Time |
|------|-------|--------|------|
| Doublet Detection | 4 | 32GB | 2h per sample |
| Outlier Removal | 4 | 16GB | 30m per sample |
| Batch Correction | 32 | 256GB | 4-8h |
| Clustering | 16 | 128GB | 2h |

## Outputs

| Step | Output |
|------|--------|
| QC | `*_qc_filtered.h5ad` |
| Batch Correction | `batch_corrected.h5ad` |
| Clustering | `postClustering*.h5ad` |
| Annotation | `postAnnotation*.h5ad` |

## Known Issues

> **⚠️ Memory Requirements**: The full dataset requires 256GB+ RAM for batch correction. Use chunked processing if limited.

> **⚠️ Multiple Cohorts**: Ensure you're using the correct input paths for ACE, Resilient, or SocIsl cohorts.

> **⚠️ Hardcoded Paths**: Many scripts have hardcoded paths. Update for your environment.

## Dependencies

### Python
```
scanpy>=1.9.0
anndata>=0.8.0
harmony-pytorch
scrublet
```

### R
```
Seurat>=4.0
DoubletFinder
```

