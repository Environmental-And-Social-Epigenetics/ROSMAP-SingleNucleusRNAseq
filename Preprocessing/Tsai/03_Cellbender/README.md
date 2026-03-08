# 03_Cellbender

Standalone CellBender pipeline that runs on CellRanger outputs produced by
`Preprocessing/Tsai/02_Cellranger_Counts`.

## Overview

This directory provides a Python-based batch script generator plus helper
scripts to submit and track CellBender jobs. It reads the patient list and
batch assignments from:

`Preprocessing/Tsai/02_Cellranger_Counts/Tracking/batch_assignments.csv`

## Scripts

| Script | Description |
|--------|-------------|
| `Scripts/generate_batch_scripts.py` | Generate SLURM scripts for each patient |
| `Scripts/run_batch.sh` | Submit all CellBender jobs for one batch |
| `Scripts/run_all_batches.sh` | Submit all batches in sequence |
| `Scripts/check_status.sh` | Show progress from tracking files |

## Configuration

All defaults are defined in:

`Config/cellbender_config.sh`

This file sources `config/paths.sh` and defines:
- Input/output paths
- CellBender parameters
- SLURM resource settings

Override any value by exporting an environment variable before running the
generator (e.g., `CB_FPR=0.02`).

## Usage

```bash
# Generate scripts from batch assignments
python Scripts/generate_batch_scripts.py

# Optional: only generate scripts for completed CellRanger samples
python Scripts/generate_batch_scripts.py --only-completed

# Submit a single batch
./Scripts/run_batch.sh 1

# Submit all batches
./Scripts/run_all_batches.sh

# Check progress
./Scripts/check_status.sh
```

## Inputs

CellRanger raw counts (per patient):

```
${TSAI_CELLRANGER_OUTPUT}/{projid}/outs/raw_feature_bc_matrix.h5
```

## Outputs

CellBender output (scratch):

```
${TSAI_CELLBENDER_SCRATCH}/{projid}/cellbender_output.h5
```

Final outputs (permanent):

```
${TSAI_PREPROCESSED}/{projid}/
├── cellbender_output.h5
└── cellbender_output_filtered.h5
```

## Tracking

Tracking files are stored in:

`Preprocessing/Tsai/03_Cellbender/Tracking/`

- `cellbender_completed.txt`
- `cellbender_failed.txt`

