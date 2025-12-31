# DeJager Processing

Processing scripts for the DeJager ROSMAP snRNA-seq dataset.

## Overview

This directory contains scripts for QC, normalization, batch correction, and cell type annotation of the DeJager snRNA-seq data. The workflow was initially developed for the Tsai dataset and adapted for DeJager.

## Scripts

### Pipeline Scripts

| Script | Description |
|--------|-------------|
| `firstStageTsaiPipeline.py` | Initial QC and preprocessing |
| `secondStageTsaiPipeline.Rscript` | R-based processing (Seurat) |
| `thirdStageTsaiPipeline.py` | Final processing and integration |

### QC Scripts

| Script | Description |
|--------|-------------|
| `doublets.Rscript` | DoubletFinder for doublet detection |
| `onePatientScript.py` | Per-patient processing workflow |

### Batch Correction

| Script | Description |
|--------|-------------|
| `batchCorrect.py` | Harmony batch correction |

## Workflow

### Stage 1: Initial Processing

**Script:** `firstStageTsaiPipeline.py`

1. Load CellBender outputs
2. Add patient IDs from demuxlet results
3. Basic QC filtering
4. Merge samples

### Stage 2: R Processing

**Script:** `secondStageTsaiPipeline.Rscript`

1. Create Seurat object
2. Run DoubletFinder
3. Seurat normalization and HVG selection
4. PCA and clustering

### Stage 3: Integration

**Script:** `thirdStageTsaiPipeline.py`

1. Load processed objects
2. Integrate across samples
3. Final clustering and UMAP
4. Cell type annotation

## Individual Analyses

The `Individual_Analyses/` directory contains per-patient analyses.

Each patient directory (e.g., `R1015854/`) contains:
- `seurat_object_*.rds`: Processed Seurat objects
- QC plots (e.g., `FeatureScatterQC_noDF.pdf`, `PercentMito_noDF.pdf`)

## Usage

### Running the Pipeline

```bash
# Stage 1
python firstStageTsaiPipeline.py

# Stage 2
Rscript secondStageTsaiPipeline.Rscript

# Stage 3
python thirdStageTsaiPipeline.py
```

### Per-Patient Analysis

```bash
python onePatientScript.py --patient R1015854
```

## Known Issues

> **⚠️ Script Naming**: Scripts are named "TsaiPipeline" for historical reasons but work with DeJager data.

> **⚠️ Demuxlet Integration**: The pipeline expects demuxlet results to be available. Ensure preprocessing is complete before running.

> **⚠️ Memory Requirements**: Processing all libraries together requires significant memory. Consider batching if resources are limited.

## Dependencies

### Python
- scanpy
- anndata
- pandas
- numpy

### R
- Seurat
- DoubletFinder
- dplyr

## Output Locations

| Output | Path |
|--------|------|
| Per-patient analyses | `Individual_Analyses/{patientID}/` |
| Integrated data | (Define based on your workflow) |

