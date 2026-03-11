# Scripts

Scripts for the standalone CellBender ambient RNA removal pipeline. These operate on Cell Ranger outputs that have already been produced by the `02_Cellranger_Counts` pipeline.

## Script Inventory

| Script | Description |
|--------|-------------|
| `generate_batch_scripts.py` | Reads batch assignments from `02_Cellranger_Counts/Tracking/batch_assignments.csv` and generates per-patient CellBender SLURM scripts under `Batch_Scripts/`. Sources `Config/cellbender_config.sh` for environment variables. |
| `run_batch.sh` | Submits all CellBender scripts for a single batch. Usage: `./run_batch.sh <batch_num>`. |
| `run_all_batches.sh` | Iterates over all batch directories and submits CellBender jobs for each. |
| `run_array.sh` | SLURM array job alternative that processes patients needing reruns as a single array submission (max 10 concurrent tasks). Reads the list of patients from `Tracking/needs_rerun.txt`. |
| `check_status.sh` | Reports pipeline progress: total patients, completed, failed, and remaining counts based on the tracking files. |

## Typical Workflow

```bash
# 1. Source configuration
source Config/cellbender_config.sh

# 2. Generate batch scripts (if not already generated)
python Scripts/generate_batch_scripts.py

# 3. Run all batches
bash Scripts/run_all_batches.sh

# 4. Check progress
bash Scripts/check_status.sh
```
