# 03_Cellbender

Remove ambient RNA from Cell Ranger outputs using CellBender.

## Overview

CellBender uses a deep generative model to distinguish true cell-associated RNA from ambient background contamination.

## Scripts

| Script | Description |
|--------|-------------|
| `DeJager_Cellbender.ipynb` | Jupyter notebook to generate batch scripts |
| `example_cellbender.sh` | Example batch script for a single library |

## Prerequisites

### CellBender Environment

```bash
source config/paths.sh
init_conda
conda activate "${CELLBENDER_ENV}"
```

### GPU Access

CellBender requires GPU acceleration:
```bash
#SBATCH --gres=gpu:a100:1
```

## Usage

### 1. Generate batch scripts

Open and run `DeJager_Cellbender.ipynb` to generate per-library scripts.

### 2. Submit jobs

```bash
sbatch example_cellbender.sh
```

## Key Parameters

```bash
cellbender remove-background \
    --cuda \                      # Use GPU
    --input <RAW_MATRIX.h5> \     # Cell Ranger raw output
    --fpr 0 \                     # False positive rate (stringent)
    --output <OUTPUT.h5>          # Output path
```

### False Positive Rate (FPR)

- `--fpr 0`: Most stringent, removes all detected ambient RNA
- `--fpr 0.01`: Allows 1% false positive rate
- Higher FPR retains more signal but also more noise

## Input

Cell Ranger raw counts:
```
${DEJAGER_COUNTS}/{LibraryID}/outs/raw_feature_bc_matrix.h5
```

## Output

CellBender-corrected counts:
```
${DEJAGER_PREPROCESSED}/{LibraryID}/
├── processed_feature_bc_matrix.h5           # Corrected counts
├── processed_feature_bc_matrix_filtered.h5  # Filtered version
└── processed_feature_bc_matrix.pdf          # QC plots
```

## Resource Requirements

| Parameter | Value |
|-----------|-------|
| Cores | 32 |
| Memory | 128GB |
| Time | 47 hours |
| GPU | A100 |

## Known Issues

> **⚠️ GPU Memory**: Large libraries may require more GPU memory. A100 is recommended.

> **⚠️ Output Directory**: The script changes to the output directory before running. Ensure directory exists.

> **⚠️ Path Configuration**: All paths are configured in `config/paths.sh`. Run `source config/paths.sh && check_paths` to verify.

## Troubleshooting

### CUDA out of memory

Reduce `--total-droplets-included` or use a GPU with more memory.

### Convergence issues

Try adjusting `--epochs` or `--learning-rate`.

