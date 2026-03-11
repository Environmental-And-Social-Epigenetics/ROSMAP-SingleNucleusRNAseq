# Scripts

Pipeline scripts for Cell Ranger counting and CellBender ambient RNA removal across all 480 Tsai patients. All scripts are designed to be called from the parent directory (`02_Cellranger_Counts/`), not from within this directory.

## Script Inventory

### Core Pipeline

| Script | Description |
|--------|-------------|
| `generate_batch_scripts.py` | Reads `Tracking/patient_metadata.csv` and generates per-patient SLURM scripts for both Cell Ranger and CellBender, organized into batch directories under `Batch_Scripts/`. Must source `Config/cellranger_config.sh` first so environment variables are set. |
| `run_batch.sh` | Runs a single batch end-to-end: submits Cell Ranger jobs, waits for completion, submits CellBender jobs, waits, then cleans up scratch. Supports `--cellranger-only`, `--cellbender-only`, and `--cleanup-only` flags. |
| `run_all_batches.sh` | Master orchestrator that iterates over batches 1-16 sequentially. Accepts `--start-batch`, `--end-batch`, and `--skip-cellranger-batch` options to resume interrupted runs. |
| `run_all_at_once.sh` | Alternative to `run_all_batches.sh` that submits all Cell Ranger scripts first, waits for all to finish, then submits all CellBender scripts. No per-batch scratch management. |
| `pipeline_slurm_wrapper.sh` | SLURM job that wraps `run_all_batches.sh` so the entire pipeline can run as a single submitted job. Automatically requeues on preemption and resumes from the last incomplete batch using tracking files. |
| `cleanup_batch.sh` | Verifies CellBender outputs are copied to permanent storage, then deletes Cell Ranger and CellBender temporary files from scratch. Supports `--force` for automated runs. |

### Recovery and Fixes

| Script | Description |
|--------|-------------|
| `retry_failed.sh` | Reads failed patient IDs from `Tracking/cellranger_failed.txt` and `Tracking/cellbender_failed.txt`, locates their SLURM scripts, and resubmits them. Supports `--cellranger-only`, `--cellbender-only`, and `--dry-run`. |
| `resubmit_fixed.sh` | One-off script that resubmits specific patients whose Cell Ranger scripts were regenerated after fixes (e.g., chemistry conflicts, path corrections). |
| `fix_220311Tsa_patients.py` | Creates symlink directories for patients whose FASTQs reference the problematic `220311Tsa` source path, then regenerates their Cell Ranger scripts with corrected paths. |
| `fix_chemistry_conflicts.sh` | SLURM script that handles patients with mixed 10x chemistry versions (SC3Pv2 + SC3Pv3) by running each library separately with `cellranger count`, then aggregating with `cellranger aggr`. |
| `fix_20151388.sh` | Patient-specific SLURM fix for patient 20151388, which has libraries with conflicting chemistries (D19-2458 as SC3Pv2, D19-4144 as SC3Pv3). Runs each library independently. |

## Typical Workflow

```bash
# 1. Source configuration
source Config/cellranger_config.sh

# 2. Generate batch scripts
python Scripts/generate_batch_scripts.py

# 3. Run all batches (or submit as SLURM job via pipeline_slurm_wrapper.sh)
bash Scripts/run_all_batches.sh
```
