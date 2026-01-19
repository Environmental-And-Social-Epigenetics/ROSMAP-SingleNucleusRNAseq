# 02_Cellranger_Counts

Run Cell Ranger `count` and CellBender on all 480 Tsai patients in batches.

## Overview

This pipeline processes all Tsai patients through:
1. **Cell Ranger count**: Align reads and generate count matrices
2. **CellBender**: Remove ambient RNA contamination (GPU-accelerated)

Due to the 1TB scratch space limit, patients are processed in **16 batches of 30 patients** each.

## Directory Structure

```
02_Cellranger_Counts/
├── Config/
│   └── cellranger_config.sh      # Pipeline configuration
├── Scripts/
│   ├── generate_batch_scripts.py # Generate all batch scripts
│   ├── run_batch.sh              # Run a single batch
│   ├── run_all_batches.sh        # Master orchestrator for all batches
│   ├── pipeline_slurm_wrapper.sh # SLURM wrapper for full automation
│   ├── cleanup_batch.sh          # Cleanup scratch after batch
│   ├── retry_failed.sh           # Retry failed jobs
│   └── fix_220311Tsa_patients.py # Fix problematic FASTQ directories
├── Batch_Scripts/                # Generated per-patient scripts
│   ├── batch_1/
│   │   ├── cellranger/           # 30 Cell Ranger scripts
│   │   └── cellbender/           # 30 CellBender scripts
│   ├── batch_2/
│   └── ... (batch_3 through batch_16)
├── Logs/
│   ├── Outs/                     # SLURM stdout
│   └── Errs/                     # SLURM stderr
├── Tracking/
│   ├── batch_assignments.csv     # Patient -> batch mapping
│   ├── patient_metadata.csv      # Full patient metadata
│   ├── cellranger_completed.txt  # Completed Cell Ranger jobs
│   ├── cellbender_completed.txt  # Completed CellBender jobs
│   ├── cellranger_failed.txt     # Failed Cell Ranger jobs
│   └── cellbender_failed.txt     # Failed CellBender jobs
└── Old/                          # Legacy notebooks
```

## Quick Start

### 1. Generate batch scripts

```bash
cd /orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/ROSMAP-SingleNucleusRNAseq/Preprocessing/Tsai/02_Cellranger_Counts

# Activate Python environment
source /orcd/data/lhtsai/001/om2/mabdel03/miniforge3/etc/profile.d/conda.sh
conda activate /orcd/data/lhtsai/001/om2/mabdel03/conda_envs/python_data_analysis

# Generate scripts (dry-run first)
python Scripts/generate_batch_scripts.py --dry-run

# Generate scripts for real
python Scripts/generate_batch_scripts.py

# Apply fixes for problematic patients
python Scripts/fix_220311Tsa_patients.py
```

### 2. Run the full pipeline (recommended)

The easiest way to run all 480 patients is using the automated orchestrator:

```bash
# Submit the master pipeline job
sbatch Scripts/pipeline_slurm_wrapper.sh
```

This will:
- Automatically process all 16 batches
- Skip Cell Ranger for batches where it's already complete
- Run CellBender on GPU nodes
- Clean up scratch space between batches
- Handle preemption with automatic requeue

### 3. Run batches manually (alternative)

```bash
# Run batch 1 (Cell Ranger + CellBender + cleanup)
./Scripts/run_batch.sh 1

# Or run steps separately:
./Scripts/run_batch.sh 1 --cellranger-only
./Scripts/run_batch.sh 1 --cellbender-only
./Scripts/run_batch.sh 1 --cleanup-only
```

### 4. Monitor progress

```bash
# Check running jobs
squeue -u $USER

# Check completed patients
wc -l Tracking/cellranger_completed.txt
wc -l Tracking/cellbender_completed.txt

# Check failures
cat Tracking/cellranger_failed.txt
cat Tracking/cellbender_failed.txt

# Check scratch usage
du -sh /home/mabdel03/orcd/scratch/Tsai/

# View pipeline master log (if using automated orchestrator)
tail -f Logs/Outs/pipeline_master_*.out
```

## Batch Workflow

Each batch follows this workflow:

```
┌─────────────────────────────────────────────────────────────┐
│ Batch N                                                     │
├─────────────────────────────────────────────────────────────┤
│ 1. Submit 30 Cell Ranger jobs (parallel, mit_preemptable)   │
│    └── Output: /home/mabdel03/orcd/scratch/Tsai/Cellranger/ │
│                                                             │
│ 2. Wait for all Cell Ranger jobs to complete                │
│                                                             │
│ 3. Submit 30 CellBender jobs (parallel, mit_normal_gpu)     │
│    └── Input: Cell Ranger raw_feature_bc_matrix.h5          │
│    └── Output: /home/mabdel03/orcd/scratch/Tsai/Cellbender/ │
│                                                             │
│ 4. Wait for all CellBender jobs to complete                 │
│                                                             │
│ 5. Copy final .h5 files to permanent storage                │
│    └── /orcd/data/.../Data/Tsai/Preprocessed_Counts/        │
│                                                             │
│ 6. Delete scratch files for this batch                      │
│    └── Frees ~750GB for next batch                          │
└─────────────────────────────────────────────────────────────┘
```

## Scripts Reference

| Script | Description |
|--------|-------------|
| `generate_batch_scripts.py` | Generate all Cell Ranger and CellBender scripts |
| `run_batch.sh` | Run a single batch (CR + CB + cleanup) |
| `run_all_batches.sh` | Master orchestrator for multiple batches |
| `pipeline_slurm_wrapper.sh` | SLURM wrapper with auto-resume capability |
| `cleanup_batch.sh` | Cleanup scratch after a batch (supports `--force`) |
| `retry_failed.sh` | Retry failed jobs from tracking files |
| `fix_220311Tsa_patients.py` | Create symlinks for problematic 220311Tsa patients |

## Cell Ranger Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `--create-bam` | false | No BAM needed (patient assignment known) |
| `--include-introns` | true | Include intronic reads (nuclear RNA) |
| `--nosecondary` | true | Skip secondary analysis |
| `--r1-length` | 26 | Trim R1 to 26bp |
| `--localcores` | 16 | Cores for Cell Ranger |
| `--localmem` | 64 | Memory (GB) for Cell Ranger |

## CellBender Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `--expected-cells` | 5000 | Expected number of cells |
| `--total-droplets-included` | 20000 | Total droplets to analyze |
| `--fpr` | 0.01 | False positive rate |
| `--epochs` | 150 | Training epochs |
| `--cuda` | (flag) | Use GPU acceleration |

## Resource Requirements

### Cell Ranger (per job)
- Partition: `mit_preemptable`
- Time: 2 days
- CPUs: 16
- Memory: 64GB
- Disk: ~5-10GB per patient

### CellBender (per job)
- Partition: `mit_normal_gpu`
- Time: 4 hours
- CPUs: 4
- Memory: 64GB
- GPU: 1
- Disk: ~1-2GB per patient

### Total Estimates
- Total patients: 480
- Batches: 16 (30 patients each)
- Time per batch: ~12-24h (Cell Ranger) + ~4h (CellBender)
- Total time: ~4-10 days (with queue wait times)
- Final storage: ~240GB (CellBender outputs only)

## Input Data

### On Engaging (Source)
FASTQ locations are indexed in:
```
ROSMAP-SingleNucleusRNAseq/Data/Tsai/All_ROSMAP_FASTQs.csv
```

This CSV contains 5,197 FASTQ files across 480 patients.

### On Openmind (Processing)
FASTQs are transferred via Globus to:
```
/om/scratch/Mon/mabdel03/Tsai_Data/FASTQs/<projid>/<Library_ID>/
```

**Note:** For Cell Ranger on Openmind, point `--fastqs` to the patient directory. Cell Ranger will automatically discover all FASTQs in subdirectories.

See `01_FASTQ_Location/03_Globus_Transfer/README.md` for transfer details.

## Output

Final CellBender outputs are stored in:
```
/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/Tsai/Preprocessed_Counts/
├── {projid_1}/
│   ├── cellbender_output.h5
│   └── cellbender_output_filtered.h5
├── {projid_2}/
│   └── ...
└── ...
```

## Known Patient Fixes

### 220311Tsa Patients (13 patients)
These patients have FASTQs in a flat directory structure that Cell Ranger cannot scan properly. The `fix_220311Tsa_patients.py` script creates symlink directories for these patients.

Affected patients: 2518573, 10310236, 21000054, 22396591, 27586957, 34542628, 38264019, 42099069, 50111971, 62574447, 64336939, 65002324, 65499271

### Patient 10490993
This patient has an invalid FASTQ directory in the source CSV. The Cell Ranger script was manually fixed to exclude the invalid directory.

## Troubleshooting

### Job failures

Check the error logs:
```bash
cat Logs/Errs/cellranger_{projid}_*.err
cat Logs/Errs/cellbender_{projid}_*.err
```

### Retry failed patients

```bash
# Dry-run first
./Scripts/retry_failed.sh --dry-run

# Actually retry
./Scripts/retry_failed.sh

# Retry only Cell Ranger failures
./Scripts/retry_failed.sh --cellranger-only

# Retry only CellBender failures
./Scripts/retry_failed.sh --cellbender-only
```

### Scratch space full

If scratch fills up unexpectedly:
```bash
# Check what's using space
du -sh /home/mabdel03/orcd/scratch/Tsai/*

# Manually run cleanup for completed batches
./Scripts/cleanup_batch.sh 1 --force
```

### Pipeline preempted

If the master pipeline job gets preempted, it will automatically requeue and resume from where it left off using the tracking files.

## Legacy Files

Old cohort-specific notebooks are preserved in `Old/`:
- `Count_TsaiACE.ipynb`
- `Count_TsaiResilient.ipynb`
- `example_count_*.sh`
