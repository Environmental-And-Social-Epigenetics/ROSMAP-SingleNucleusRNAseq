# Scripts

Scripts for transferring Tsai FASTQ files between HPC clusters via Globus.

## Script Inventory

| Script | Description |
|--------|-------------|
| `generate_globus_batch.py` | Reads `All_ROSMAP_FASTQs.csv` and generates a Globus batch transfer file mapping each FASTQ from its source location on Engaging to the destination on Openmind, organized by `projid/Library_ID/run_id`. Handles samples with duplicate filenames from different sequencing runs by creating unique subdirectories per source. Output goes to `Batch_Files/`. |
| `submit_globus_transfer.sh` | Activates the Globus CLI conda environment and submits the batch transfer between the Engaging and Openmind endpoints. Supports `--check-only` to validate login status and batch file without submitting. |
| `verify_transfer.py` | Compares expected files (from `All_ROSMAP_FASTQs.csv`) against actual files at the destination after transfer completes. Reports missing and extra files with summary statistics. Supports `--generate-manifest` for pre-transfer validation and `--dest-dir` to check a specific directory. |

## Typical Workflow

```bash
# 1. Generate the batch transfer file
python Scripts/generate_globus_batch.py

# 2. Submit the transfer
bash Scripts/submit_globus_transfer.sh

# 3. After transfer completes, verify
python Scripts/verify_transfer.py
```
