# DeJager Preprocessing

Preprocessing pipeline for the DeJager ROSMAP snRNA-seq dataset.

## Overview

The DeJager dataset consists of multiplexed snRNA-seq libraries downloaded from Synapse. Each library contains cells from multiple patients, requiring genotype-based demultiplexing.

## Pipeline Steps

```
Synapse → FASTQs → Cell Ranger → CellBender → Demuxlet → Patient-assigned cells
```

### 01_FASTQ_Download

Download FASTQ files from Synapse using the Python synapse client.

**Scripts:**
- `Download_FASTQs.py`: Main download script using synapse SDK
- `Download_FASTQs.sh`: SLURM batch wrapper
- `Download_FASTQs_Rerun.py/sh`: Retry failed downloads

**Prerequisites:**
```bash
source config/paths.sh
init_conda
conda activate "${SYNAPSE_ENV}"
```

**Input:**
- `Synapse_FASTQ_IDs.csv`: CSV with Synapse IDs for each FASTQ file

**Output:**
- FASTQs organized by library ID in `${DEJAGER_FASTQS}` (see `config/paths.sh`)

### 02_Cellranger_Counts

Run Cell Ranger `count` on each library.

**Scripts:**
- `Count_DeJager.py`: Generates batch scripts for all libraries and submits them
- `Count_DeJager.sh`: SLURM wrapper for the Python script
- `example_count.sh`: Example batch script for a single library

**Key Parameters:**
- Cores: 32
- Memory: 128GB
- Time: 47 hours
- `--create-bam true`: Generate BAM (required for Demuxlet)
- `--include-introns true`: Include intronic reads

**Output:**
- Counts in `${DEJAGER_COUNTS}/{LibraryID}/` (see `config/paths.sh`)

### 03_Cellbender

Remove ambient RNA from Cell Ranger outputs.

**Scripts:**
- `DeJager_Cellbender.ipynb`: Jupyter notebook to generate batch scripts
- `example_cellbender.sh`: Example batch script

**Key Parameters:**
- Cores: 32
- Memory: 128GB
- Time: 47 hours
- GPU: A100
- `--fpr 0`: False positive rate (stringent)

**Output:**
- `processed_feature_bc_matrix.h5` in `${DEJAGER_PREPROCESSED}` (see `config/paths.sh`)

### 04_Demuxlet_Freemuxlet

Assign cells to patients using WGS genotype data.

**Workflow:** BAM filtering -> Pileup generation -> Demuxlet -> Aggregate assignments

**Scripts:**
- `Demuxlet_DeJager.py`: Generates batch scripts for all libraries and submits them
- `Demuxlet_DeJager.sh`: SLURM wrapper for the Python script
- `example_demuxlet.sh`: Example batch script for a single library
- `postprocess_assignments.py`: Aggregate results into patient assignment CSV
- `validate_demuxlet.py`: Validate demuxlet results against annotations

**Tools:**
- Uses Demuxafy Singularity container (`${DEMUXAFY_SIF}`)
- Popscle toolkit for pileup and demuxlet
- Bundled `popscle_helper_tools/` for BAM filtering

**Prerequisites:**
```bash
source config/paths.sh
# Verify: ${DEMUXAFY_SIF}, ${DEJAGER_DEMUX_VCF}, ${DEJAGER_PATIENT_IDS_DIR}
```

**Input:**
- Possorted BAM from Cell Ranger (Step 02)
- Cell barcodes from CellBender (Step 03)
- WGS VCF and patient ID lists (external data)

**Output:**
- Per-library `.best` files with cell-patient assignments
- Aggregated CSV at `${DEJAGER_PATIENT_MAP}`

**Documentation:**
- `docs/PARAMETER_TUNING.md`: Parameter sweep results and rationale
- `docs/VCF_PREPARATION.md`: How the WGS VCF was prepared

## Known Issues

> **Path Configuration**: All paths are configured in `config/paths.sh`. Run `source config/paths.sh && check_paths` to verify your setup.

## Data Files

| File | Description |
|------|-------------|
| `Synapse_FASTQ_IDs.csv` | Maps Synapse IDs to library IDs |
| `individPat*.txt` | Patient IDs for each library (in `${DEJAGER_PATIENT_IDS_DIR}`) |
| `snp_fixedconcatenated_liftedROSMAP.vcf.gz` | SNP VCF for demuxlet (`${DEJAGER_DEMUX_VCF}`) |

## Resource Requirements

| Step | Cores | Memory | Time | GPU |
|------|-------|--------|------|-----|
| FASTQ Download | 64 | 512GB | 47h | - |
| Cell Ranger | 32 | 128GB | 47h | - |
| CellBender | 32 | 128GB | 47h | A100 |
| Demuxlet (BAM filter) | 45 | 400GB | 3h | - |
| Demuxlet (pileup+demux) | 10 | 500GB | 36h | - |

