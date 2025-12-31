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
    ├── QC/
    │   ├── Doublets/                # Doublet detection scripts
    │   └── Outliers/                # Outlier removal scripts
    ├── Batch_Correction/            # Harmony batch correction
    ├── Norm_FeatSel_DimRed_Cluster.ipynb  # Main processing notebook
    └── Tsai_Sample_SingleCell_Pipeline.ipynb  # Sample pipeline
```

## Workflow

### 1. Quality Control

#### Doublet Detection

Uses DoubletFinder (R) to identify and remove doublets.

**Scripts:**
- `Master_Remove_Doublets.py`: Generate per-sample doublet scripts
- `Master_Remove_Doublets.r`: R implementation
- `example_doublet.r`: Example per-sample script

**Method:**
1. Preprocess with Seurat
2. Run DoubletFinder with expected doublet rate
3. Filter cells classified as doublets

#### Outlier Removal

Remove cells with extreme QC metrics.

**Scripts:**
- `Master_Remove_Outliers.py`: Generate per-sample outlier scripts
- `example_outlier.py`: Example per-sample script

**Criteria:**
- Low gene count (< 200 genes)
- High mitochondrial percentage (> 20%)
- Extreme UMI counts

### 2. Normalization and Feature Selection

**Script:** `Norm_FeatSel_DimRed_Cluster.ipynb`

**Steps:**
1. Size factor normalization
2. Log transformation
3. Highly variable gene (HVG) selection
4. Scaling

### 3. Batch Correction

Uses Harmony to integrate samples and remove batch effects.

**Scripts:**
- `Batch_Correction.py`: Main batch correction script
- `Batch_Correction.sh`: SLURM wrapper
- `Batch_Correction2.py`: Alternative implementation
- `featureSelOnward.py`/`featureSelOnward2.py`: Post-correction processing

**Method:**
1. Run PCA on normalized data
2. Apply Harmony on PC space
3. Use corrected embeddings for downstream analysis

### 4. Dimensionality Reduction & Clustering

**Scripts:**
- `batchCluster.py`: Clustering after batch correction
- `makeUMAPS.py`: Generate UMAP visualizations

**Steps:**
1. PCA on HVG
2. Harmony batch correction
3. UMAP visualization
4. Leiden/Louvain clustering

### 5. Cell Type Annotation

**Scripts:**
- `anno.py`: Annotation script
- `anno.Rscript`: R-based annotation

**Method:**
- Marker-based annotation using known brain cell type markers
- Reference markers in `/Analysis/Tsai/Gene_Markers/`

## Key Notebooks

| Notebook | Purpose |
|----------|---------|
| `Norm_FeatSel_DimRed_Cluster.ipynb` | Main processing pipeline |
| `Tsai_Sample_SingleCell_Pipeline.ipynb` | Sample/test pipeline |
| `BatchCorrection_Test.ipynb` | Batch correction testing |

## Resource Requirements

| Step | Cores | Memory | Time |
|------|-------|--------|------|
| Doublet Detection | 4 | 32GB | 2h |
| Outlier Removal | 4 | 32GB | 1h |
| Batch Correction | 32 | 256GB | 4-8h |

## Known Issues

> **⚠️ Pipeline Naming**: The DeJager scripts are named "TsaiPipeline" because the workflow was developed on Tsai data first. The pipeline works for both datasets.

> **⚠️ Resource Intensive**: Batch correction on the full dataset requires substantial memory (256GB+). Consider processing in chunks if memory is limited.

## Outputs

| File | Description |
|------|-------------|
| `*_doublet_filtered.h5ad` | Doublet-removed data |
| `*_qc_filtered.h5ad` | QC-filtered data |
| `batch_corrected.h5ad` | Harmony-corrected integrated data |
| `postClustering*.h5ad` | Clustered data |
| `postAnnotation*.h5ad` | Annotated data |

## Next Steps

After processing, proceed to the `Analysis/` directory for differential expression and downstream analyses.

