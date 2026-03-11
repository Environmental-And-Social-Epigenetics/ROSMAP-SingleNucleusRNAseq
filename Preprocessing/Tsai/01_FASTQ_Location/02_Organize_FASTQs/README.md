# 02_Organize_FASTQs

This directory contains scripts to organize FASTQ files by patient (projid), creating a directory structure with symbolic links.

## Overview

These scripts read the master FASTQ CSV (`All_ROSMAP_FASTQs.csv`) and create a organized directory structure where each patient has their own subdirectory containing symlinks to all their FASTQ files.

## Scripts

### `run_organize_fastqs.sh`
SLURM batch script that orchestrates the parallel organization process.

**Usage:**
```bash
sbatch run_organize_fastqs.sh
```

**Prerequisites:**
- Must run `01_Build_Master_CSV/run_build_csv.sh` first to generate `All_ROSMAP_FASTQs.csv`

**What it does:**
1. Reads the master FASTQ CSV
2. Preprocesses to assign occurrence numbers (handles multiple sequencing runs per Library_ID)
3. Distributes work across parallel jobs by projid
4. Calls `organize_patient.sh` for each projid
5. Runs validation via `validate_organization.py`

### `organize_patient.sh`
Worker script that creates the directory structure and symlinks for a single patient.

**Usage:**
```bash
./organize_patient.sh <projid> <fastq_csv> <output_dir>
```

**What it does:**
- Creates `{output_dir}/{projid}/`
- For each Library_ID belonging to this projid:
  - Creates subdirectory `{Library_ID}_N/` (N = occurrence number, 1, 2, etc.)
  - Creates symbolic links to all FASTQs for that Library_ID

### `validate_organization.py`
Python script that validates the organized directory structure.

**Usage:**
```bash
python3 validate_organization.py --fastq-csv All_ROSMAP_FASTQs.csv --output-dir FASTQs_By_Patient --report validation_report.txt
```

**What it checks:**
1. All expected patient directories exist
2. Symlinks are created and not broken
3. File counts match the master CSV
4. Generates a detailed validation report

## Output Directory Structure

```
FASTQs_By_Patient/
├── 20244554/
│   └── D19-10914_1_1/
│       ├── D19-10914_1_S1_L001_I1_001.fastq.gz -> /nfs/...
│       ├── D19-10914_1_S1_L001_R1_001.fastq.gz -> /nfs/...
│       └── D19-10914_1_S1_L001_R2_001.fastq.gz -> /nfs/...
├── 68745332/
│   ├── D19-4111_1/
│   │   └── *.fastq.gz -> /nfs/...
│   └── D19-4111_2/
│       └── *.fastq.gz -> /nfs/...
└── ...
```

The naming convention `{Library_ID}_N` matches the original notebook convention:
- `_1` = first occurrence/sequencing run
- `_2` = second occurrence/sequencing run
- etc.

## Validation Report

The validation report (`validation_report.txt`) includes:
- Total expected vs actual file counts
- Per-patient validation status
- List of any issues found
- Overall pass/fail verdict

**Validation criteria:**
- PASSED: ≥95% coverage, no broken symlinks
- PASSED WITH WARNINGS: ≥80% coverage
- FAILED: <80% coverage or broken symlinks

## Configuration

All paths and parameters are defined in `../Config/config.sh`.

## Troubleshooting

**Missing symlinks:**
- Check if source files exist at the original paths
- Verify backup paths are accessible
- Review the `All_ROSMAP_FASTQs_missing.csv` for files that couldn't be found

**Broken symlinks:**
- Source files may have been moved or deleted
- Check filesystem access permissions
- Re-run `01_Build_Master_CSV` to update the index

