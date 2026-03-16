# ROSMAP Single Nucleus RNA Sequencing Pipeline

This repository contains the single nucleus RNA sequencing (snRNA-seq) analysis pipeline for the Religious Orders Study and Memory and Aging Project (ROSMAP) dataset. The pipeline processes dorsolateral prefrontal cortex (DLPFC) samples from two primary data sources:

- **DeJager Dataset**: Data downloaded from Synapse, requiring sample demultiplexing via Demuxlet/Freemuxlet
- **Tsai Dataset**: Data located on MIT Engaging cluster, with known patient assignments

## Getting Started

1. **Understand the project** — Read [BACKGROUND.md](BACKGROUND.md) for scientific context (ROSMAP, snRNA-seq, and why each pipeline step matters)
2. **Set up your environment** — Follow [setup/README.md](setup/README.md) for first-time cluster setup
3. **Run Preprocessing** — See [Preprocessing/](Preprocessing/README.md) (choose DeJager or Tsai pathway)
4. **Run Processing** — Submit all three stages with one command:
   ```bash
   cd Processing/Tsai/Pipeline
   ./submit_pipeline.sh all
   ```
5. **Downstream Analysis** — See [Analysis/](Analysis/README.md)

### For New Lab Members

1. Get MIT Engaging cluster access: `ssh <Kerberos ID>@orcd-login003.mit.edu`
2. Clone this repo to your workspace
3. Copy `config/paths.local.sh.template` to `config/paths.local.sh` and fill in your paths
4. Run `bash setup/install_envs.sh` to create conda environments
5. Download data from the Tsai Lab NAS (see [Data_Access/](Data_Access/README.md))
6. Run the pipeline: `cd Processing/Tsai/Pipeline && ./submit_pipeline.sh all`

## Data Provenance (March 2026)

The processed data from Mahmoud Abdelmoneum, Nina Khera, and Ravikiran Raju's March 2026 analysis is available at these locations. **Openmind was decommissioned March 17, 2026 — all data now lives on Engaging.**

| Data | Primary Location | Backup | Access Method |
|------|-----------------|--------|---------------|
| Tsai raw FASTQs | Engaging cluster | NAS | Already on Engaging (CSV mapping in repo) |
| DeJager raw FASTQs | Synapse ([syn21438684](https://www.synapse.org/#!Synapse:syn21438684)) | Engaging, NAS | Synapse download or Engaging |
| Tsai CellBender outputs (478 samples) | Engaging | NAS | Already on Engaging |
| DeJager CellBender outputs (131 dirs) | Engaging | NAS | Already on Engaging |
| Tsai cell-type annotated object | Engaging (`tsai_annotated.h5ad`, 83GB) | NAS | Already on Engaging |
| **DeJager processing outputs** | **Not yet generated** | — | **Run pipeline on Engaging** |
| Phenotype CSVs | In repo (`Data/Phenotypes/`) | Engaging, NAS | Already included |
| DeJager WGS VCF + Demuxafy SIF | Engaging (`/home/nkhera/orcd/pool/WGS/`) | — | Already on Engaging |

See [Data_Access/README.md](Data_Access/README.md) for the full Engaging data layout, transfer scripts, and Globus audit trail.

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
├── Processing/              # QC, doublet removal, batch correction, and annotation
│   ├── DeJager/
│   │   ├── Pipeline/              # 3-stage pipeline (mirrors Tsai)
│   │   │   ├── 01_qc_filter.*
│   │   │   ├── 02_doublet_removal.*
│   │   │   ├── 03_integration_annotation.*
│   │   │   ├── submit_pipeline.sh
│   │   │   ├── Resources/
│   │   │   └── envs/
│   │   └── _legacy/               # Archived original scripts
│   └── Tsai/
│       ├── Pipeline/              # 3-stage pipeline (primary development target)
│       │   ├── 01_qc_filter.*
│       │   ├── 02_doublet_removal.*
│       │   ├── 03_integration_annotation.*
│       │   ├── submit_pipeline.sh # SLURM submission wrapper
│       │   ├── Resources/         # Marker gene references
│       │   └── envs/              # Conda environment specs
│       └── archive/               # Superseded legacy scripts
├── Data/                    # Pipeline data (phenotypes tracked, transcriptomics gitignored)
│   ├── Phenotypes/          # ROSMAP clinical, ACE scores, ID maps (tracked in git)
│   └── Transcriptomics/     # FASTQs, CellRanger, CellBender outputs (populate via Data_Access/)
├── Data_Access/             # Scripts to download/transfer data (NAS, Globus, Synapse)
└── Analysis/                # Downstream analysis (organized by phenotype)
    ├── _template/           # Template for adding new phenotype analyses
    ├── ACE/                 # Adverse Childhood Experiences
    │   ├── DEG/  (DeJager/, Tsai/)
    │   ├── TF/   (DeJager/, Tsai/)
    │   └── SCENIC/ (DeJager/, Tsai/)
    ├── Resilient/           # Cognitive Resilience
    │   └── (same structure)
    └── SocIsl/              # Social Isolation
        └── (same structure)
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

Both datasets use an identical three-stage pipeline (`Processing/{Dataset}/Pipeline/`):

1. **Stage 1 — QC Filtering** (`01_qc_filter.py`): Percentile-based outlier removal and mitochondrial % filtering
2. **Stage 2 — Doublet Removal** (`02_doublet_removal.Rscript`): scDblFinder-based doublet detection
3. **Stage 3 — Integration & Annotation** (`03_integration_annotation.py`): Normalization, HVG selection, PCA, Harmony batch correction, Leiden clustering, UMAP, and ORA-based cell type annotation with Mohammadi 2020 markers

Submit all stages with dependency chaining:

```bash
cd Processing/Tsai/Pipeline && ./submit_pipeline.sh all     # Tsai (478 samples with CellBender outputs)
cd Processing/DeJager/Pipeline && ./submit_pipeline.sh all   # DeJager
```

See `Processing/Tsai/Pipeline/README.md` or `Processing/DeJager/Pipeline/README.md` for full details.

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

Before running the pipeline, create your local path overrides:

```bash
cp config/paths.local.sh.template config/paths.local.sh
# Edit config/paths.local.sh with your cluster-specific paths
source config/paths.sh
check_paths
```

Validate that all pipeline prerequisites (conda environments, input files, references) are in place:

```bash
bash config/preflight.sh
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
| CellBender | mit_normal_gpu | 32 | 128GB | 47h | A100 |
| Demuxlet | mit_normal | 80 | 400GB | 48h | - |

### Tsai Pipeline (Updated)

| Step | Partition | Cores | Memory | Time | GPU |
|------|-----------|-------|--------|------|-----|
| Cell Ranger | mit_preemptable | 16 | 64GB | 2 days | - |
| CellBender | mit_normal_gpu | 4 | 64GB | 4h | 1 |

### Processing Pipeline (Both Datasets)

| Stage | Cores | Memory | Time | Notes |
|-------|-------|--------|------|-------|
| 1 — QC Filtering | 4 | 32GB | 12h | Array job (Tsai: 478 tasks, DeJager: 200 tasks, 32 concurrent) |
| 2 — Doublet Removal | 4 | 32GB | 12h | Array job (same dimensions as Stage 1) |
| 3 — Integration & Annotation | 32 | 500GB | 48h | Single job (loads all samples) |

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
| FASTQs | `TSAI_FASTQS_DIR` | `${SCRATCH_ROOT}/Tsai_Data/FASTQs` |
| FASTQ Index | `TSAI_FASTQS_CSV` | `${DATA_ROOT}/Data/Tsai/Preprocessing/FASTQ_Transfer/New/CSVs/All_ROSMAP_FASTQs.csv` |
| Cell Ranger | `TSAI_CELLRANGER_OUTPUT` | `${SCRATCH_ROOT}/Tsai_Data/Cellranger_Outputs` |
| CellBender (temp) | `TSAI_CELLBENDER_SCRATCH` | `${SCRATCH_ROOT}/Tsai/Cellbender_Output` |
| Preprocessed | `TSAI_PREPROCESSED` | `${WORKSPACE_ROOT}/Tsai_Data/Cellbender_Outputs` |

**Dataset:** 480 patients, 5,197 FASTQ files, processed in 16 batches of 30 patients each.

## Known Issues

See `KNOWN_ISSUES.md` for a detailed tracker. Key remaining items:

1. **Scratch dependencies**: Intermediate Cell Ranger/CellBender data on scratch may be cleaned before processing completes
2. **Pipeline naming**: Legacy scripts in `Processing/DeJager/_legacy/` are named "TsaiPipeline" (historical artifact)
3. **SocIsl preprocessing**: SocIsl cohort has batch scripts but no dedicated notebooks for Cell Ranger/CellBender script generation

## Contributing

When adding new scripts or analyses:

1. Place scripts in the appropriate phase/dataset directory
2. Update the relevant README.md with usage instructions
3. Use relative paths or configuration files where possible
4. Document any new dependencies

## Contact

For questions about this pipeline, contact:

- Mahmoud Abdelmoneum (mabdel03@mit.edu)
- Nina Khera (nkhera@college.harvard.edu)
- Ravikiran Raju (rraju@mit.edu)

Environmental and Social Epigenetics Lab, Tsai Lab, MIT.
