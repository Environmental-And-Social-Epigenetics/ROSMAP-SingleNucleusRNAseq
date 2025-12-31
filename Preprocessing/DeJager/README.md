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
conda activate /orcd/data/lhtsai/001/om2/mabdel03/conda_envs/synapse_env
```

**Input:**
- `Synapse_FASTQ_IDs.csv`: CSV with Synapse IDs for each FASTQ file

**Output:**
- FASTQs organized by library ID in `/om/scratch/Mon/mabdel03/FASTQs/`

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
- Counts in `/om/scratch/Mon/mabdel03/Counts/{LibraryID}/`

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
- `processed_feature_bc_matrix.h5` in `/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/DeJager/Preprocessed_Counts/`

### 04_Demuxlet_Freemuxlet

Assign cells to patients using WGS genotype data.

**Scripts:**
- `demuxTest.sh`: Main demuxlet workflow (extensive parameter testing)
- `demux.sh`: Simplified demuxlet script
- `freemux.sh`: Freemuxlet alternative (no reference genotypes)
- `pileup.sh`: Generate pileup files
- `demuxValCode.py`: Validate demuxlet results
- `freemuxletValidateCode.py`: Validate freemuxlet results

**Tools:**
- Uses Demuxafy Singularity container
- Popscle toolkit for pileup and demuxlet

**Input:**
- Possorted BAM from Cell Ranger
- Filtered VCF with patient genotypes
- Cell barcodes from CellBender

**Output:**
- `.best` file with cell-patient assignments

## Known Issues

> **⚠️ Hardcoded Paths**: Scripts reference `/orcd/data/lhtsai/001/om2/mabdel03/` and `/om/scratch/Mon/mabdel03/`. Update these for your environment.

> **⚠️ demuxTest.sh**: Contains extensive commented-out parameter tuning experiments. The active commands are at lines 53-54 and 69. Consider cleaning up for production use.

## Data Files

| File | Description |
|------|-------------|
| `Synapse_FASTQ_IDs.csv` | Maps Synapse IDs to library IDs |
| `individualsTest*.txt` | Patient IDs for each library |
| `filteredSNPFreq2.vcf.gz` | Filtered VCF for demuxlet |

## Resource Requirements

| Step | Cores | Memory | Time | GPU |
|------|-------|--------|------|-----|
| FASTQ Download | 64 | 512GB | 47h | - |
| Cell Ranger | 32 | 128GB | 47h | - |
| CellBender | 32 | 128GB | 47h | A100 |
| Demuxlet | 80 | 400GB | 48h | - |

