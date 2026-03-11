# Tracking

Completion-tracking files for the standalone CellBender pipeline.

## Files

| File | Description |
|------|-------------|
| `cellbender_completed.txt` | Patient IDs (`projid`, one per line) that finished CellBender successfully. Appended by each SLURM job on completion. |
| `cellbender_completed.txt.bak` | Backup of the completed list, preserved before bulk edits or reruns. |
| `cellbender_failed.txt` | Patient IDs whose CellBender jobs failed. Used by recovery scripts. |
| `needs_rerun.txt` | Patient IDs queued for reprocessing, consumed by `Scripts/run_array.sh` as a SLURM array job manifest. |

## Notes

- These files are append-only during normal operation.
- To reprocess a patient, add its `projid` to `needs_rerun.txt` and remove it from `cellbender_completed.txt`.
- `Scripts/check_status.sh` reads these files to report pipeline progress.
