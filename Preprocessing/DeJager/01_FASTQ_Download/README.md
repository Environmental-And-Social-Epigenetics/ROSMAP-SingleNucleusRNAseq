# 01_FASTQ_Download

Download FASTQ files from Synapse for the DeJager ROSMAP snRNA-seq dataset.

## Overview

This step downloads raw FASTQ files from Synapse (syn21438684) using the Synapse Python client.

## Scripts

| Script | Description |
|--------|-------------|
| `Download_FASTQs.py` | Main download script |
| `Download_FASTQs.sh` | SLURM batch wrapper |
| `Download_FASTQs_Rerun.py` | Retry failed downloads |
| `Download_FASTQs_Rerun.sh` | SLURM wrapper for rerun |

## Prerequisites

### Synapse Account

1. Create account at [synapse.org](https://www.synapse.org)
2. Request access to syn21438684
3. Generate personal access token

### Environment

```bash
conda activate /orcd/data/lhtsai/001/om2/mabdel03/conda_envs/synapse_env
```

## Usage

### 1. Configure Synapse credentials

```bash
synapse login --auth-token <YOUR_TOKEN>
```

### 2. Submit download job

```bash
sbatch Download_FASTQs.sh
```

### 3. Monitor progress

```bash
squeue -u $USER
```

## Input

- `Synapse_FASTQ_IDs.csv`: CSV mapping Synapse IDs to library IDs
  - Located at: `/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/DeJager/FASTQs_Download/FASTQ_Download_CSVs/`

## Output

FASTQs organized by library ID:
```
/om/scratch/Mon/mabdel03/FASTQs/
├── 190403-B4-A/
│   ├── 190403-B4-A_Broad_S1_L001_R1_001.fastq.gz
│   ├── 190403-B4-A_Broad_S1_L001_R2_001.fastq.gz
│   └── ...
├── 190403-B4-B/
└── ...
```

## Resource Requirements

| Parameter | Value |
|-----------|-------|
| Cores | 64 |
| Memory | 512GB |
| Time | 47 hours |

## Known Issues

> **⚠️ Large Downloads**: The full dataset is several terabytes. Ensure sufficient scratch space.

> **⚠️ Network Timeouts**: Use the rerun script if downloads fail.

> **⚠️ Hardcoded Paths**: Update paths in scripts for your environment.

