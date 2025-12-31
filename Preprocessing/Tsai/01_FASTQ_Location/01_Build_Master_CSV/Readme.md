# 01_Build_Master_CSV

This directory contains scripts to build a master CSV indexing all FASTQ files for ROSMAP patients on Engaging.

## Overview

These scripts read the existing `Tsai_To_Openmind.csv` mapping (797 rows, 483 unique projids) and enumerate all FASTQ files at each path, creating a comprehensive index.

## Scripts

### `run_build_csv.sh`
SLURM batch script that orchestrates the parallel FASTQ discovery process.

**Usage:**
```bash
sbatch run_build_csv.sh
```

**What it does:**
1. Reads the master mapping CSV (`Tsai_To_Openmind.csv`)
2. Loads backup paths from `Reformatted_Backups.csv`
3. Splits work across parallel jobs using GNU `parallel`
4. Calls `list_fastqs_for_path.sh` for each unique path
5. Merges results and runs `merge_and_finalize_csv.py`

### `list_fastqs_for_path.sh`
Worker script that lists all FASTQ files for a given path and Library_ID.

**Usage:**
```bash
./list_fastqs_for_path.sh <row_index> <projid> <library_id> <tsai_path> <backup_path> <output_file>
```

**Logic:**
- Searches the primary `tsai_path` for `*.fastq.gz` files matching the `library_id`
- Falls back to `backup_path` if primary is empty or missing
- Outputs CSV rows: `projid,Library_ID,source_dir,fastq_filename,full_path,file_size`

### `merge_and_finalize_csv.py`
Python script that merges worker outputs and creates the final master CSV.

**Usage:**
```bash
python3 merge_and_finalize_csv.py --input combined.csv --output All_ROSMAP_FASTQs.csv --master-csv Tsai_To_Openmind.csv
```

**What it does:**
1. Reads combined worker outputs
2. Deduplicates entries
3. Validates all file paths exist
4. Compares with master CSV to identify missing projids
5. Outputs `All_ROSMAP_FASTQs.csv`

## Output

The main output is `New/CSVs/All_ROSMAP_FASTQs.csv` with columns:
- `projid`: Patient identifier
- `Library_ID`: Sequencing library ID (e.g., D19-4111)
- `source_dir`: Directory containing the FASTQ
- `fastq_filename`: Name of the FASTQ file
- `full_path`: Complete path to the file
- `file_size`: File size in bytes

## Configuration

All paths and parameters are defined in `../Config/config.sh`. Edit that file to customize:
- Input CSV locations
- Output directories
- Parallelization settings
- SLURM parameters

## Prerequisites

- GNU `parallel` must be available
- Python 3 with `pandas` installed
- Access to the Engaging `/nfs/picower*` filesystems

