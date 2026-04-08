# ROSMAP Single Nucleus RNA Sequencing Pipeline

This repository packages the ROSMAP single-nucleus RNA-seq workflows for two
sequencing sets:

- **Tsai**: FASTQs live on Engaging and libraries already map to individual donors
- **DeJager**: FASTQs come from Synapse and multiplexed libraries require donor assignment

The intended flow is:

1. raw FASTQs
2. preprocessing to CellBender-corrected counts
3. a shared three-stage processing pipeline
4. downstream phenotype analysis

The current cleanup pass standardizes environments, removes hardcoded user paths,
keeps generated outputs outside the repo, and formalizes ACE analysis for both
cohorts.

## Supported Workflows

| Area | Status | Notes |
|------|--------|-------|
| Tsai preprocessing | `official` | FASTQ discovery, Cell Ranger batching, CellBender |
| DeJager preprocessing | `official` | Synapse download, Cell Ranger, scripted CellBender, Demuxlet/Freemuxlet |
| Tsai processing | `official` | Canonical Stage 3 Harmony key: `derived_batch` |
| DeJager processing | `official` | Canonical Stage 3 Harmony key: `library_id`; not yet run (see KNOWN_ISSUES #19) |
| ACE DEG | `official` | Tsai and DeJager; outputs under `ANALYSIS_OUTPUT_ROOT` |
| ACE cell-type proportion | `official` | sccomp workflow for Tsai and DeJager |
| ACE SCENIC | `implemented` | Ported from SocIsl; requires SCENIC ranking databases |
| ACE TF/COMPASS | `implemented` | Ported from SocIsl; requires IBM CPLEX license |
| ACE GSEA | `implemented` | WebGestaltR pathway enrichment |
| Resilient DEG | `implemented` | Ported from ACE with resilience group derivation |
| Resilient cell-type proportion | `implemented` | sccomp with resilience groups |
| SocIsl DEG | `migrated` | Paths migrated to config/paths.sh |
| SocIsl SCENIC | `migrated` | Paths migrated; uses `${SCENIC_RANKING_DIR}` |
| SocIsl TF/COMPASS | `migrated` | Paths migrated; uses `${CPLEX_DIR}` |
| SocIsl GSEA | `migrated` | Paths migrated |
| SocIsl cell-type proportion | `implemented` | Ported from ACE |

## Analysis Status Matrix

| | DEG | SCENIC | TF/COMPASS | GSEA | CellTypeProportion |
|---|---|---|---|---|---|
| **ACE** | official | implemented | implemented | implemented | official |
| **SocIsl** | migrated | migrated | migrated | migrated | implemented |
| **Resilient** | implemented | scaffold | scaffold | scaffold | implemented |

## Quick Start

1. Read [setup/README.md](setup/README.md).
2. Copy `config/paths.local.sh.template` to `config/paths.local.sh`.
3. Run:
   ```bash
   source config/paths.sh
   check_paths
   bash setup/install_envs.sh --all
   ```
4. Run preprocessing for the cohort you need:
   - [Preprocessing/Tsai/README.md](Preprocessing/Tsai/README.md)
   - [Preprocessing/DeJager/README.md](Preprocessing/DeJager/README.md)
5. Run processing:
   ```bash
   cd Processing/Tsai/Pipeline && ./submit_pipeline.sh all
   cd Processing/DeJager/Pipeline && ./submit_pipeline.sh all
   ```
6. Run downstream analysis:
   - [Analysis/README.md](Analysis/README.md)
   - [Analysis/ACE/README.md](Analysis/ACE/README.md)

## Repository Layout

```text
Transcriptomics/
├── config/            central path and environment configuration
├── setup/             first-time setup and environment installation
├── Preprocessing/     FASTQ -> CellBender workflows
├── Processing/        QC, doublet removal, Harmony integration, annotation
├── Analysis/          phenotype analyses
├── Data/Phenotypes/   tracked phenotype and ID-map inputs
└── Data_Access/       download / transfer helpers
```

Large generated outputs are expected to live **outside** the repo tree. By
default:

- processing outputs go under the configured `*_PROCESSING_OUTPUTS` roots
- downstream analysis outputs go under `ANALYSIS_OUTPUT_ROOT`

## Pipeline Summary

### Preprocessing

| Step | Tsai | DeJager |
|------|------|---------|
| FASTQ access | indexed on Engaging | download from Synapse |
| Alignment | Cell Ranger `count` | Cell Ranger `count` |
| Ambient RNA removal | CellBender | CellBender |
| Donor assignment | known from metadata | Demuxlet/Freemuxlet + overrides |

### Processing

Both cohorts use the same three-stage structure:

1. **Stage 1**: percentile-based QC filtering
2. **Stage 2**: `scDblFinder` doublet removal
3. **Stage 3**: normalization, HVG selection, PCA, Harmony, clustering, ORA annotation

Canonical batch variables:

- **Tsai**: `derived_batch`
- **DeJager**: `library_id`

The ORA annotation step uses the shared Mohammadi 2020 PFC marker reference and
the raw gene space in both cohorts.

### Analysis

Three phenotype analyses are implemented:

- **ACE**: pseudobulk DEG, cell-type proportion (sccomp), SCENIC, TF/COMPASS, GSEA
- **Resilient**: pseudobulk DEG, cell-type proportion (sccomp)
- **SocIsl**: DEG, SCENIC, TF/COMPASS, GSEA, cell-type proportion (migrated from legacy)

See [Analysis/ACE/README.md](Analysis/ACE/README.md), [Analysis/Resilient/README.md](Analysis/Resilient/README.md),
and [Analysis/SocIsl/README.md](Analysis/SocIsl/README.md) for phenotype contracts.

## Environment Setup

Every official environment now lives in its own directory and includes:

- `environment.yml`
- `requirements.txt`
- `README.md`

Use the installer:

```bash
source config/paths.sh
bash setup/install_envs.sh --all
```

Or choose an install method explicitly:

```bash
bash setup/install_envs.sh --analysis --method=conda
bash setup/install_envs.sh --analysis --method=requirements
```

`requirements.txt` is a companion artifact for every official env. For hybrid
Python/R or system-dependent environments, the env README explains the required
conda bootstrap before `pip install -r requirements.txt`.

## Data Notes

- `Preprocessing/Tsai/02_Cellranger_Counts/Tracking/patient_metadata.csv` lists
  **480** Tsai samples.
- The current CellBender-complete Tsai preprocessing set contains **478**
  sample directories.
- Missing Tsai CellBender outputs: `11467746`, `20834164`.

## Configuration

The shared path contract lives in [config/paths.sh](config/paths.sh). Local
machine- or user-specific overrides belong in `config/paths.local.sh`.

Important variables:

- `TSAI_INTEGRATED`, `TSAI_INTEGRATED_PROJID`
- `DEJAGER_INTEGRATED`, `DEJAGER_INTEGRATED_PATIENT_ID`, `DEJAGER_INTEGRATED_POOL_BATCH`, `DEJAGER_INTEGRATED_DERIVED_BATCH`
- `ACE_SCORES_CSV`
- `ANALYSIS_OUTPUT_ROOT`
- `NEBULA_ENV`, `SCCOMP_ENV`

## Next Reads

- [BACKGROUND.md](BACKGROUND.md)
- [setup/README.md](setup/README.md)
- [Processing/README.md](Processing/README.md)
- [Analysis/README.md](Analysis/README.md)
