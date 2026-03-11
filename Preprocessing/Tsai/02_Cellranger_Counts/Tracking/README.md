# Tracking

Metadata and completion-tracking files used by the Cell Ranger + CellBender batch pipeline.

## Files

### Master Data

| File | Description |
|------|-------------|
| `patient_metadata.csv` | Master patient list (480 rows + header). Columns: `projid`, `library_ids`, `fastq_dirs`, `n_source_dirs`, `n_fastqs`, `batch`, `index`. This file is the single source of truth for which patients to process and where their FASTQs are located. |
| `batch_assignments.csv` | Patient-to-batch mapping extracted from the metadata. Columns: `projid`, `library_ids`, `batch`, `index`. Used by `Scripts/run_batch.sh` and `Scripts/generate_batch_scripts.py` to determine which patients belong to each batch. |

### Completion Tracking

These plain-text files contain one `projid` per line and are updated by the SLURM job scripts as each patient completes or fails. The pipeline uses them to skip already-processed patients on restart.

| File | Description |
|------|-------------|
| `cellranger_completed.txt` | Patient IDs that finished Cell Ranger successfully. |
| `cellranger_failed.txt` | Patient IDs whose Cell Ranger jobs failed. Fed to `Scripts/retry_failed.sh`. |
| `cellbender_completed.txt` | Patient IDs that finished CellBender successfully. |
| `cellbender_failed.txt` | Patient IDs whose CellBender jobs failed. Fed to `Scripts/retry_failed.sh`. |

## Notes

- The tracking files are append-only during normal operation.
- `pipeline_slurm_wrapper.sh` reads these files to find the next incomplete batch after preemption.
- To reprocess a patient, remove its `projid` from the relevant completed/failed file before resubmitting.
