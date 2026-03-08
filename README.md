# ROSMAP Single Nucleus RNA Sequencing Pipeline

This repository contains the single nucleus RNA sequencing (snRNA-seq) analysis pipeline for the Religious Orders Study and Memory and Aging Project (ROSMAP) dataset. The pipeline processes dorsolateral prefrontal cortex (DLPFC) samples from two primary data sources:

- **DeJager Dataset**: Data downloaded from Synapse, requiring sample demultiplexing via Demuxlet/Freemuxlet
- **Tsai Dataset**: Data located on MIT Engaging cluster, with known patient assignments

## Repository Structure

```
ROSMAP-SingleNucleusRNAseq/
├── config/                  # Central path configuration (paths.sh)
├── setup/                   # First-time setup guide and install scripts
├── Preprocessing/           # Raw data processing (FASTQs → CellBender outputs)
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
│       ├── Pipeline/              # Primary 3-stage pipeline (current)
│       │   ├── 01_qc_filter.*
│       │   ├── 02_doublet_removal.*
│       │   ├── 03_integration_annotation.*
│       │   ├── submit_pipeline.sh # SLURM submission wrapper
│       │   ├── Resources/         # Marker gene references
│       │   └── envs/              # Conda environment specs
│       └── archive/               # Superseded legacy scripts
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

### Phase 2: Processing (Tsai Pipeline)

The primary processing pipeline lives in `Processing/Tsai/Pipeline/` and runs in three stages:

1. **Stage 1 — QC Filtering** (`01_qc_filter.py`): Percentile-based outlier removal and mitochondrial % filtering
2. **Stage 2 — Doublet Removal** (`02_doublet_removal.Rscript`): scDblFinder-based doublet detection
3. **Stage 3 — Integration & Annotation** (`03_integration_annotation.py`): Normalization, HVG selection, PCA, Harmony batch correction, Leiden clustering, UMAP, and ORA-based cell type annotation with Mohammadi 2020 markers

See `Processing/Tsai/Pipeline/README.md` for full details.

### Phase 3: Analysis

Downstream biological analyses.

1. **Differential Expression**: DEG analysis between conditions/cell types
2. **SCENIC**: Single-cell regulatory network inference
3. **Transcription Factor Analysis**: TF activity and regulatory analysis

## Prerequisites

### Software Requirements

- **Cell Ranger** v8.0.0
- **CellBender** (GPU-accelerated)
- **Python** 3.10+ with: scanpy, anndata, harmonypy, decoupler
- **R** 4.2+ with: scDblFinder, zellkonverter, BiocParallel
- **Singularity** (for Demuxafy container, DeJager only)

Conda environment specs are provided in `Processing/Tsai/Pipeline/envs/`.

### Path Configuration

Before running the pipeline, configure paths in `config/paths.sh` (or create
`config/paths.local.sh` for per-user overrides):

```bash
source config/paths.sh
check_paths
```

See `setup/README.md` for a complete first-time setup guide.

### Conda Environments

Create the processing pipeline environments using the automated installer:

```bash
source config/paths.sh
bash setup/install_envs.sh
```

Or manually from the YAML specs in `Processing/Tsai/Pipeline/envs/`.
After creation, use the variables from `config/paths.sh`:

```bash
source config/paths.sh
init_conda
conda activate "${QC_ENV}"          # Stage 1
conda activate "${SINGLECELL_ENV}"  # Stage 2
conda activate "${BATCHCORR_ENV}"   # Stage 3
```

### Reference Files

- Human reference genome: `refdata-gex-GRCh38-2020-A`
- Set `CELLRANGER_REF` in `config/paths.sh` to the reference directory

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

See `KNOWN_ISSUES.md` for a detailed tracker. Key remaining items:

1. **Demuxlet scripts**: Contains commented-out parameter tuning code that should be cleaned up
2. **Scratch dependencies**: Intermediate data on scratch may be cleaned up before processing completes
3. **DeJager processing scripts**: Named "TsaiPipeline" (historical artifact) — rename pending

## Contributing

When adding new scripts or analyses:

1. Place scripts in the appropriate phase/dataset directory
2. Update the relevant README.md with usage instructions
3. Use relative paths or configuration files where possible
4. Document any new dependencies

## Contact

For questions about this pipeline, contact the Environmental and Social Epigenetics lab.
