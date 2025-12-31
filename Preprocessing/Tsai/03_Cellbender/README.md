# 03_Cellbender

Remove ambient RNA from Cell Ranger outputs using CellBender.

## Overview

CellBender uses a deep generative model to distinguish true cell-associated RNA from ambient background contamination.

## Scripts

| Script | Description |
|--------|-------------|
| `ACE_Cellbender.ipynb` | Generate batch scripts for ACE cohort |
| `Resilient_Cellbender.ipynb` | Generate batch scripts for Resilient cohort |
| `example_cellbender_Resilient.sh` | Example batch script |

## Prerequisites

### CellBender Environment

```bash
source /orcd/data/lhtsai/001/om2/mabdel03/miniforge3/etc/profile.d/conda.sh
conda activate /orcd/data/lhtsai/001/om2/mabdel03/conda_envs/Cellbender_env
```

### GPU Access

CellBender requires GPU acceleration:
```bash
#SBATCH --gres=gpu:a100:1
```

## Usage

### 1. Generate batch scripts

Open and run the appropriate Jupyter notebook:
- `ACE_Cellbender.ipynb` for ACE cohort
- `Resilient_Cellbender.ipynb` for Resilient cohort

### 2. Submit jobs

```bash
sbatch example_cellbender_Resilient.sh
```

## Key Parameters

```bash
cellbender remove-background \
    --cuda \                      # Use GPU
    --input <RAW_MATRIX.h5> \     # Cell Ranger raw output
    --fpr 0 \                     # False positive rate (stringent)
    --output <OUTPUT.h5>          # Output path
```

## Input

Cell Ranger raw counts:
```
/om/scratch/Mon/mabdel03/Tsai/{cohort}/Counts/{projid}/outs/raw_feature_bc_matrix.h5
```

## Output

CellBender-corrected counts:
```
/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/Tsai/Preprocessing/Preprocessed_Counts/{cohort}/{projid}/
├── processed_feature_bc_matrix.h5           # Corrected counts
├── processed_feature_bc_matrix_filtered.h5  # Filtered version
└── processed_feature_bc_matrix.pdf          # QC plots
```

## Cohort Paths

| Cohort | Output Path |
|--------|-------------|
| ACE | `Preprocessed_Counts/ACE/` |
| Resilient | `Preprocessed_Counts/Resilient/` |
| SocIsl | `Preprocessed_Counts/SocIsl/` |

## Resource Requirements

| Parameter | Value |
|-----------|-------|
| Cores | 32 |
| Memory | 500GB |
| Time | 47 hours |
| GPU | A100 |

> **Note**: Tsai samples may require more memory (500GB) compared to DeJager (128GB).

## Known Issues

> **⚠️ GPU Memory**: Large samples may require more GPU memory.

> **⚠️ Output Directory**: The script changes to output directory before running. Ensure it exists.

> **⚠️ Hardcoded Paths**: Update input/output paths in notebooks for your environment.

