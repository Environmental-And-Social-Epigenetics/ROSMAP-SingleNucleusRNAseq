# Tsai Preprocessing

Preprocessing pipeline for the Tsai ROSMAP snRNA-seq dataset.

## Overview

The Tsai dataset contains snRNA-seq data from the Tsai Lab, already located on the MIT Engaging cluster. Unlike the DeJager dataset, patient assignments are known from sequencing metadata, so no demultiplexing is required.

## Dataset Summary

| Metric | Value |
|--------|-------|
| Total Patients | 480 |
| Total FASTQ Files | 5,197 |
| Processing Batches | 16 (30 patients each) |

## Pipeline Steps

```
Locate FASTQs -> Cell Ranger -> CellBender -> Patient-assigned cells
```

### 01_FASTQ_Location

Locate and index FASTQ files across the Engaging filesystem.

**Key Output:**
- `All_ROSMAP_FASTQs.csv`: Master index of all 5,197 FASTQ files

**Columns:**

| Column | Description |
|--------|-------------|
| `projid` | Patient identifier |
| `Library_ID` | Sequencing library ID |
| `source_dir` | Directory containing FASTQ |
| `fastq_filename` | Filename |
| `full_path` | Complete path |
| `file_size` | Size in bytes |

### 02_Cellranger_Counts

Run Cell Ranger `count` and CellBender on all 480 patients.

**Key Features:**
- Automated batched processing (16 batches of 30 patients)
- GPU-accelerated CellBender
- Automatic scratch space management
- Failure tracking and retry capability

**Quick Start:**

```bash
cd 02_Cellranger_Counts

# Generate batch scripts
python Scripts/generate_batch_scripts.py

# Apply fixes for problematic patients
python Scripts/fix_220311Tsa_patients.py

# Run the full pipeline
sbatch Scripts/pipeline_slurm_wrapper.sh
```

See `02_Cellranger_Counts/README.md` for detailed documentation.

### 03_Cellbender

CellBender is now integrated into the `02_Cellranger_Counts` pipeline. This directory contains legacy cohort-specific notebooks.

## Resource Requirements

| Step | Partition | Cores | Memory | Time | GPU |
|------|-----------|-------|--------|------|-----|
| FASTQ Location | mit_normal | 4 | 8GB | 1h | - |
| Cell Ranger | mit_preemptable | 16 | 64GB | 2 days | - |
| CellBender | mit_normal_gpu | 4 | 64GB | 4h | 1 |

## Data Locations

### Input FASTQs

FASTQs are located across multiple directories on the Engaging filesystem. The master index is stored at:

```
ROSMAP-SingleNucleusRNAseq/Data/Tsai/All_ROSMAP_FASTQs.csv
```

### Scratch Space (Temporary)

```
/home/mabdel03/orcd/scratch/Tsai/
    Cellranger_Counts/    # Cell Ranger outputs (temporary)
    Cellbender_Output/    # CellBender outputs (temporary)
```

### Permanent Storage

```
/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/Tsai/Preprocessed_Counts/
    {projid}/
        cellbender_output.h5
        cellbender_output_filtered.h5
```

## Known Issues

### Hardcoded Paths

Scripts reference `/orcd/data/lhtsai/001/om2/mabdel03/` paths. Update `config/paths.sh` for your environment.

### Problematic FASTQ Directories

- **220311Tsa patients** (13 patients): FASTQs in flat directory structure that Cell Ranger cannot scan. Fixed via symlinks in `Data/Tsai/FASTQ_Symlinks/`.
- **Patient 10490993**: Invalid FASTQ directory in source CSV. Manually fixed in batch script.

### Scratch Space Limit

The 1TB scratch quota requires careful batch management. The pipeline automatically cleans scratch between batches.

## Monitoring Progress

```bash
# Check job queue
squeue -u $USER

# Check completion stats
cd 02_Cellranger_Counts
wc -l Tracking/cellranger_completed.txt Tracking/cellbender_completed.txt

# Check failures
cat Tracking/cellranger_failed.txt
cat Tracking/cellbender_failed.txt

# View pipeline log
tail -f Logs/Outs/pipeline_master_*.out
```

## Prerequisites

- Python 3.8+ with `pandas`
- Cell Ranger 8.0.0
- CellBender conda environment
- Access to `/nfs/picower*` filesystems

### Conda Environments

```bash
# Initialize conda
source /orcd/data/lhtsai/001/om2/mabdel03/miniforge3/etc/profile.d/conda.sh

# Python data analysis (for scripts)
conda activate /orcd/data/lhtsai/001/om2/mabdel03/conda_envs/python_data_analysis

# CellBender
conda activate /orcd/data/lhtsai/001/om2/mabdel03/conda_envs/Cellbender_env
```
