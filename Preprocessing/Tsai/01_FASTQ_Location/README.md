# 01_FASTQ_Location

This directory contains scripts to locate, index, organize, and transfer FASTQ files for the Tsai ROSMAP snRNA-seq dataset.

## Overview

The pipeline discovers FASTQs scattered across the Engaging filesystem, builds a comprehensive index, and transfers them to Openmind for Cell Ranger processing.

## Directory Structure

```
01_FASTQ_Location/
├── 01_Build_Master_CSV/     # Discover and index FASTQs
├── 02_Organize_FASTQs/      # Create organized symlink structure
├── 03_Globus_Transfer/      # Transfer to Openmind via Globus
├── Config/                   # Configuration files
└── README.md                 # This file
```

## Workflow

### Step 1: Build Master CSV

Discover all FASTQs and create the master index.

```bash
cd 01_Build_Master_CSV
sbatch run_build_csv.sh
```

**Output:** `All_ROSMAP_FASTQs.csv` in `ROSMAP-SingleNucleusRNAseq/Data/Tsai/`

### Step 2: Organize FASTQs (Optional)

Create symlink structure for local use on Engaging.

```bash
cd 02_Organize_FASTQs
sbatch run_organize_fastqs.sh
```

**Output:** Organized symlinks in `FASTQs_By_Patient/`

### Step 3: Transfer to Openmind

Transfer FASTQs to Openmind for Cell Ranger processing.

```bash
cd 03_Globus_Transfer

# Generate batch file
python Scripts/generate_globus_batch.py

# Submit transfer
./Scripts/submit_globus_transfer.sh

# Verify transfer (run on Openmind after completion)
python Scripts/verify_transfer.py
```

**Output:** FASTQs in `/om/scratch/Mon/mabdel03/Tsai_Data/FASTQs/<projid>/<Library_ID>/`

## Data Summary

| Metric | Value |
|--------|-------|
| Total FASTQ files | 5,197 |
| Unique patients (projids) | 480 |
| Unique Library_IDs | 492 |
| Total data size | ~9 TB |
| Multi-source samples | 272 (sequenced on multiple flow cells) |

## Destination Structure on Openmind

```
/om/scratch/Mon/mabdel03/Tsai_Data/FASTQs/
├── <projid>/
│   └── <Library_ID>/
│       ├── [run_id/]           # For multi-source samples (e.g., 10x-4182G)
│       │   └── *.fastq.gz
│       └── *.fastq.gz          # For single-source samples
```

## Globus Endpoints

| Cluster | Endpoint ID | Name |
|---------|-------------|------|
| Engaging | `ec54b570-cac5-47f7-b2a1-100c2078686f` | MIT ORCD Engaging Collection |
| Openmind | `cbc6f8da-d37e-11eb-bde9-5111456017d9` | mithpc#openmind |

## Configuration

Edit `Config/config.sh` to customize paths and parameters.

## Prerequisites

- GNU `parallel` for FASTQ discovery
- Python 3 with `pandas`
- Globus CLI (`/home/mabdel03/conda_envs/globus_env`)
- Access to `/nfs/picower*` filesystems
