# Tsai Preprocessing

Preprocessing pipeline for the Tsai ROSMAP snRNA-seq dataset.

## Overview

The Tsai dataset contains snRNA-seq data from the Tsai Lab, already located on the MIT Engaging cluster. Unlike the DeJager dataset, patient assignments are known from sequencing metadata, so no demultiplexing is required.

## Cohorts

The Tsai data is organized into three cohorts:

| Cohort | Description |
|--------|-------------|
| **ACE** | Adverse Childhood Experiences study |
| **Resilient** | Cognitive resilience study |
| **SocIsl** | Social Isolation study |

## Pipeline Steps

```
Locate FASTQs → Cell Ranger → CellBender → Patient-assigned cells
```

### 01_FASTQ_Location

Locate and index FASTQ files across the Engaging filesystem.

**Directory Structure:**
```
01_FASTQ_Location/
├── 01_Build_Master_CSV/     # Build comprehensive FASTQ index
├── 02_Organize_FASTQs/      # Organize FASTQs by patient
└── Config/                   # Configuration files
```

**Scripts:**
- `run_build_csv.sh`: SLURM script for parallel FASTQ discovery
- `list_fastqs_for_path.sh`: Worker script to enumerate FASTQs
- `merge_and_finalize_csv.py`: Merge and validate results
- `organize_patient.sh`: Organize FASTQs by patient ID
- `config.sh`: Configuration parameters

**Output:**
- `All_ROSMAP_FASTQs.csv`: Master index with columns:
  - `projid`: Patient identifier
  - `Library_ID`: Sequencing library ID
  - `source_dir`: Directory containing FASTQ
  - `fastq_filename`: Filename
  - `full_path`: Complete path
  - `file_size`: Size in bytes

### 02_Cellranger_Counts

Run Cell Ranger `count` on each sample.

**Scripts:**
- `Count_TsaiACE.ipynb`: Generate batch scripts for ACE cohort
- `Count_TsaiResilient.ipynb`: Generate batch scripts for Resilient cohort
- `example_count_Resilient.sh`: Example batch script (Resilient)
- `example_count_SocIsl.sh`: Example batch script (SocIsl)

**Key Differences from DeJager:**
- `--create-bam false`: No BAM needed (no demuxlet)
- FASTQs may span multiple directories (comma-separated in `--fastqs`)

**Key Parameters:**
- Cores: 32
- Memory: 128GB
- Time: 47 hours
- `--include-introns true`: Include intronic reads

**Output:**
- Counts in `/om/scratch/Mon/mabdel03/Tsai/{cohort}/Counts/{projid}/`

### 03_Cellbender

Remove ambient RNA from Cell Ranger outputs.

**Scripts:**
- `ACE_Cellbender.ipynb`: Generate batch scripts for ACE cohort
- `Resilient_Cellbender.ipynb`: Generate batch scripts for Resilient cohort
- `example_cellbender_Resilient.sh`: Example batch script

**Key Parameters:**
- Cores: 32
- Memory: 500GB
- Time: 47 hours
- GPU: A100
- `--fpr 0`: False positive rate (stringent)

**Output:**
- `processed_feature_bc_matrix.h5` in `/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/Tsai/Preprocessing/Preprocessed_Counts/{cohort}/{projid}/`

## Known Issues

> **⚠️ Hardcoded Paths**: Scripts reference `/orcd/data/lhtsai/001/om2/mabdel03/` and `/om/scratch/Mon/mabdel03/`. Update these for your environment.

> **⚠️ Multiple FASTQ Directories**: Some samples have FASTQs split across multiple directories. The `--fastqs` parameter accepts comma-separated paths.

> **⚠️ Cohort Organization**: Each cohort (ACE, Resilient, SocIsl) has its own batch scripts and output directories. Ensure you're using the correct paths.

## Resource Requirements

| Step | Cores | Memory | Time | GPU |
|------|-------|--------|------|-----|
| FASTQ Location | 4 | 8GB | 1h | - |
| Cell Ranger | 32 | 128GB | 47h | - |
| CellBender | 32 | 500GB | 47h | A100 |

## Configuration

Edit `01_FASTQ_Location/Config/config.sh` to customize:

- Input CSV locations
- Output directories
- Parallelization settings
- SLURM parameters

## Prerequisites

- GNU `parallel` for FASTQ discovery
- Python 3 with `pandas`
- Access to `/nfs/picower*` filesystems

