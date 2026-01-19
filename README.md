# ROSMAP Single Nucleus RNA Sequencing Pipeline

This repository contains the single nucleus RNA sequencing (snRNA-seq) analysis pipeline for the Religious Orders Study and Memory and Aging Project (ROSMAP) dataset. The pipeline processes dorsolateral prefrontal cortex (DLPFC) samples from two primary data sources:

- **DeJager Dataset**: Data downloaded from Synapse, requiring sample demultiplexing via Demuxlet/Freemuxlet
- **Tsai Dataset**: Data located on MIT Engaging cluster, with known patient assignments

## Repository Structure

```
ROSMAP-SingleNucleusRNAseq/
├── Preprocessing/           # Raw data processing (FASTQs → Cellbender outputs)
│   ├── DeJager/            # DeJager-specific preprocessing
│   │   ├── 01_FASTQ_Download/
│   │   ├── 02_Cellranger_Counts/
│   │   ├── 03_Cellbender/
│   │   └── 04_Demuxlet_Freemuxlet/
│   └── Tsai/               # Tsai-specific preprocessing (480 patients)
│       ├── 01_FASTQ_Location/    # FASTQ discovery and indexing
│       ├── 02_Cellranger_Counts/ # Automated batch pipeline (CR + CellBender)
│       └── 03_Cellbender/        # Legacy cohort notebooks
├── Processing/              # QC and cell type annotation
│   ├── DeJager/
│   └── Tsai/
│       ├── QC/
│       │   ├── Doublets/
│       │   └── Outliers/
│       └── Batch_Correction/
└── Analysis/                # Downstream analysis
    ├── DeJager/
    │   ├── DEG_Analysis/
    │   ├── SCENIC/
    │   └── Transcription_Factors/
    └── Tsai/
        ├── DEG_Analysis/
        ├── SCENIC/
        └── Transcription_Factors/
```

## Pipeline Overview

### Phase 1: Preprocessing

Converts raw FASTQ files into ambient RNA-corrected count matrices.

| Step | DeJager | Tsai |
|------|---------|------|
| 1. Data Acquisition | Download FASTQs from Synapse | Locate FASTQs on Engaging (indexed in CSV) |
| 2. Alignment & Counting | Cell Ranger `count` | Cell Ranger `count` (batched, 30 patients/batch) |
| 3. Ambient RNA Removal | CellBender | CellBender (GPU-accelerated) |
| 4. Sample Assignment | Demuxlet/Freemuxlet (WGS-based) | Known from metadata |

### Phase 2: Processing

Quality control and preprocessing of the count matrices.

1. **Doublet Detection**: DoubletFinder (R) to identify and remove doublets
2. **Outlier Removal**: Filter cells based on QC metrics (gene count, UMI count, mitochondrial %)
3. **Normalization**: Size factor normalization and log transformation
4. **Feature Selection**: Highly variable gene selection
5. **Batch Correction**: Harmony integration across samples
6. **Dimensionality Reduction**: PCA and UMAP
7. **Clustering**: Leiden/Louvain clustering
8. **Cell Type Annotation**: Marker-based annotation

### Phase 3: Analysis

Downstream biological analyses.

1. **Differential Expression**: DEG analysis between conditions/cell types
2. **SCENIC**: Single-cell regulatory network inference
3. **Transcription Factor Analysis**: TF activity and regulatory analysis

## Prerequisites

### Software Requirements

- **Cell Ranger** v8.0.0
- **CellBender** (GPU-accelerated)
- **Python** 3.8+ with:
  - scanpy
  - anndata
  - harmony-pytorch
  - scrublet
- **R** 4.0+ with:
  - Seurat
  - DoubletFinder
- **Singularity** (for Demuxafy container)

### Path Configuration

Before running the pipeline, configure paths in `config/paths.sh`:

```bash
# Source the configuration
source config/paths.sh

# Verify paths are correct
check_paths
```

**Important**: Update `SCRATCH_ROOT` in `config/paths.sh` to point to your cluster's scratch filesystem.

### Conda Environments

All environments are located at `/orcd/data/lhtsai/001/om2/mabdel03/conda_envs/`:

```bash
# Initialize conda
source /orcd/data/lhtsai/001/om2/mabdel03/miniforge3/etc/profile.d/conda.sh

# CellBender environment
conda activate /orcd/data/lhtsai/001/om2/mabdel03/conda_envs/Cellbender_env

# Synapse download environment
conda activate /orcd/data/lhtsai/001/om2/mabdel03/conda_envs/synapse_env
```

### Reference Files

- Human reference genome: `refdata-gex-GRCh38-2020-A`
- Location: `/orcd/data/lhtsai/001/om2/mabdel03/yard/references/human/refdata-gex-GRCh38-2020-A`

## SLURM Resources

### DeJager Pipeline

| Step | Partition | Cores | Memory | Time | GPU |
|------|-----------|-------|--------|------|-----|
| Cell Ranger | mit_normal | 32 | 128GB | 47h | - |
| CellBender | mit_normal_gpu | 32 | 128-500GB | 47h | A100 |
| Demuxlet | mit_normal | 80 | 400GB | 48h | - |

### Tsai Pipeline (Updated)

| Step | Partition | Cores | Memory | Time | GPU |
|------|-----------|-------|--------|------|-----|
| Cell Ranger | mit_preemptable | 16 | 64GB | 2 days | - |
| CellBender | mit_normal_gpu | 4 | 64GB | 4h | 1 |

## Data Locations

Data paths are configured in `config/paths.sh`. Update `SCRATCH_ROOT` for your cluster.

### DeJager

| Data Type | Variable | Default Path |
|-----------|----------|--------------|
| FASTQs | `DEJAGER_FASTQS` | `${SCRATCH_ROOT}/FASTQs/` |
| Counts | `DEJAGER_COUNTS` | `${SCRATCH_ROOT}/Counts/` |
| Preprocessed | `DEJAGER_PREPROCESSED` | `${DATA_ROOT}/Data/DeJager/Preprocessed_Counts/` |

### Tsai

| Data Type | Variable | Default Path |
|-----------|----------|--------------|
| FASTQs | (from CSV) | Located on Engaging filesystem |
| FASTQ Index | `TSAI_FASTQS_CSV` | `${REPO_ROOT}/Data/Tsai/All_ROSMAP_FASTQs.csv` |
| Cell Ranger (temp) | `TSAI_CELLRANGER_SCRATCH` | `${SCRATCH_ROOT}/Tsai/Cellranger_Counts/` |
| CellBender (temp) | `TSAI_CELLBENDER_SCRATCH` | `${SCRATCH_ROOT}/Tsai/Cellbender_Output/` |
| Preprocessed | `TSAI_PREPROCESSED` | `${DATA_ROOT}/Data/Tsai/Preprocessed_Counts/` |

**Dataset:** 480 patients, 5,197 FASTQ files, processed in 16 batches of 30 patients each.

## Known Issues

> **Warning**: Several scripts contain hardcoded paths that may need to be updated for your environment. See individual directory READMEs for specific issues.

1. **Hardcoded paths**: Many scripts reference `/orcd/data/lhtsai/001/om2/mabdel03/` or `/om/scratch/Mon/mabdel03/` paths
2. **Demuxlet scripts**: Contains extensive commented-out parameter tuning code that should be cleaned up
3. **Scratch paths**: Some scripts reference temporary scratch directories that may not persist

## Contributing

When adding new scripts or analyses:

1. Place scripts in the appropriate phase/dataset directory
2. Update the relevant README.md with usage instructions
3. Use relative paths or configuration files where possible
4. Document any new dependencies

## Contact

For questions about this pipeline, contact the Environmental and Social Epigenetics lab.
