# Known Issues and Flagged Problems

This document consolidates known issues identified during the organization of this repository. These issues should be addressed before production use.

## Critical Issues

### 1. Scratch Space Paths Need Configuration

**Severity**: High
**Status**: Resolved

All scripts now source `config/paths.sh` and derive paths from configured
variables.  Preprocessing scripts (DeJager and Tsai) have been migrated.
Hardcoded Synapse auth tokens in the DeJager download scripts have been
removed (use `synapse login --rememberMe` or `SYNAPSE_AUTH_TOKEN` env var
instead).

---

### 2. Demuxlet Script Cleanup Needed

**Severity**: Medium
**Status**: Resolved

The experimental scripts have been archived to `archive/` and replaced with
clean, production-ready scripts following the pattern of Steps 01-03:
- `Demuxlet_DeJager.py`: Batch script generator (like `Count_DeJager.py`)
- `Demuxlet_DeJager.sh`: SLURM wrapper
- `example_demuxlet.sh`: Single-library example
- `postprocess_assignments.py`: Aggregates demux results into CSV
- `validate_demuxlet.py`: Parameterized validation script

Parameter tuning experiments extracted to `docs/PARAMETER_TUNING.md`.
VCF preparation documented in `docs/VCF_PREPARATION.md`.
BAM/VCF helper tools bundled in `popscle_helper_tools/`.
Demuxlet-specific path variables added to `config/paths.sh`.

---

### 3. Scratch Space Dependencies

**Severity**: Medium  
**Affected Scripts**: Cell Ranger and CellBender scripts

Scripts reference temporary scratch directories:
- `/om/scratch/Mon/mabdel03/` (weekly cleanup)
- `/om/scratch/Sun/mabdel03/` (weekly cleanup)

**Risk**: Data loss if scratch is cleaned before processing completes.

**Recommendation**: 
- Move intermediate outputs to persistent storage
- Or document scratch cleanup schedule

---

## Moderate Issues

### 4. Pipeline Naming Confusion

**Severity**: Low  
**Affected Files**: `Processing/DeJager/` scripts

Scripts are named "TsaiPipeline" (e.g., `firstStageTsaiPipeline.py`) even though they work for both datasets.

**Reason**: Pipeline was developed on Tsai data first.

**Recommendation**: Consider renaming for clarity (e.g., `stage1_processing.py`).

---

### 5. Conda Environment Paths - RESOLVED

**Severity**: Low
**Status**: Fixed

All conda environment paths are now centralized in `config/paths.sh`.
Processing pipeline shell wrappers source this config.  Environment specs for
recreating envs from scratch are in `Processing/Tsai/Pipeline/envs/`.

---

### 6. Missing SocIsl Processing Scripts

**Severity**: Medium  
**Location**: `Preprocessing/Tsai/`

SocIsl cohort has batch scripts but no dedicated notebooks for script generation.

**Affected Steps**:
- 02_Cellranger_Counts: Uses Resilient notebook pattern
- 03_Cellbender: No SocIsl notebook

**Recommendation**: Create dedicated notebooks or document the workflow for SocIsl.

---

## Informational Notes

### 7. Large Resource Requirements

Some steps require significant resources:
- Stage 3 Integration: 500GB RAM
- Demuxlet: 400GB RAM, 80 cores
- Batch Correction: 256GB RAM

Ensure adequate cluster allocation.

### 8. GPU Requirements

CellBender requires GPU:
- Recommended: A100
- Minimum: GPU with 16GB VRAM

### 9. Dependencies on External Data

- WGS VCF files for Demuxlet (not in repo)
- Reference genome (not in repo)
- Synapse credentials (not in repo)

---

### 10. `.gitignore` Blocked Essential Pipeline Files — RESOLVED

**Severity**: High
**Status**: Fixed

The blanket `*.csv` and `*.rds` rules in `.gitignore` prevented
`patient_metadata.csv` and `Brain_Human_PFC_Markers_Mohammadi2020.rds` from
being tracked.  Whitelist exceptions have been added.

---

### 11. Hemoglobin Gene Regex Bug — RESOLVED

**Severity**: Low
**Status**: Fixed

`r"^HB[^(P)]"` in `01_qc_filter.py` and `03_integration_annotation.py` used a
character class that excluded literal `(`, `)`, and `P`.  Corrected to
`r"^HB[^P]"`.

---

### 12. No QC Summary Tracking — RESOLVED

**Severity**: Medium
**Status**: Fixed

Stages 1 and 2 now produce `qc_summary.csv` and `doublet_summary.csv`
respectively, tracking per-sample cell counts before/after filtering.

---

### 13. Stale Parent READMEs — RESOLVED

**Severity**: Medium
**Status**: Fixed

`Processing/README.md`, `Processing/Tsai/README.md`, and the root `README.md`
referenced DoubletFinder/Scrublet and outdated directory structures.  Updated
to reflect the current Pipeline.

---

### 14. README Hardcoded Paths — RESOLVED

**Severity**: Medium
**Status**: Fixed

README documentation files contained hardcoded paths to
`/orcd/data/lhtsai/001/om2/mabdel03/` and `/om/scratch/Mon/mabdel03/`
despite the migration to `config/paths.sh`. All READMEs have been updated
to reference config variables instead.

---

### 15. Harmony Batch Variable Was Patient ID — RESOLVED

**Severity**: High
**Status**: Fixed

`03_integration_annotation.py` used `projid` (patient ID, 480 levels) as the
Harmony batch variable. This was too aggressive — it forced cells from all 480
patients to overlap in PC space, removing genuine inter-individual biological
variation needed for downstream pseudobulk DEG analysis. Additionally, the
`batch` column in `patient_metadata.csv` (values 1-16) was a computational
convenience grouping for SLURM jobs, not a real sequencing/wet-lab batch.

**Fix**: Created `derive_batches.py` which extracts actual Illumina flowcell IDs
from FASTQ headers and groups samples by shared flowcells (~41 natural batch
groups). The integration script now defaults to `--harmony-batch-key derived_batch`
but supports switching to `projid` or any other obs column via CLI flags.

---

### 16. Data_Access/ and Key Pipeline Files Were Untracked — RESOLVED

**Severity**: High
**Status**: Fixed

The entire `Data_Access/` directory (29 transfer scripts), all 4 phenotype CSVs
in `Data/Phenotypes/`, `derive_batches.py`, `03_integration_annotation_projid.sh`,
`03b_evaluate_correction.py`, `derived_batches.csv`, and `batch_assignments.csv`
were not tracked in git. A fresh clone was missing data download capability,
phenotype data, and batch derivation tooling.

**Fix**: Added `.gitignore` whitelist entries for phenotype CSVs and pipeline
resources. `git add`-ed all untracked files. Created `Data/Transcriptomics/README.md`
documenting the canonical data layout.

---

### 17. Hardcoded mabdel03 Paths in Production Scripts — RESOLVED

**Severity**: High
**Status**: Fixed

All defaults in `config/paths.sh` and fallback paths in Preprocessing/Processing
Python and shell scripts were hardcoded to mabdel03's Openmind environment. Globus
and NAS transfer scripts also had hardcoded conda and SFTP paths.

**Fix**: Replaced all user-specific defaults with `${HOME}`-based auto-detection,
`__UNCONFIGURED__` sentinels for must-configure values, and `${VAR:-default}`
patterns that respect environment overrides. Created `config/paths.local.sh.template`.
Archive/legacy scripts are marked as non-portable with header comments.

---

### 18. Archive Scripts Are Not Portable (Intentional)

**Severity**: Informational
**Status**: By design

Scripts under `Processing/Tsai/archive/`, `Processing/DeJager/_legacy/`,
`Preprocessing/Tsai/02_Cellranger_Counts/Old/`, and generated `Batch_Scripts/`
directories contain hardcoded paths from the original development environment.
These are preserved for historical reference only and are not intended to be
portable. The current production pipeline is in `Processing/{Dataset}/Pipeline/`.

---

## Resolution Tracking

| Issue | Status | Date Fixed |
|-------|--------|------------|
| 1. Hardcoded paths (pipeline) | Resolved | 2026-03-08 |
| 1. Hardcoded paths (preprocessing) | Resolved | 2026-03-08 |
| 2. Demuxlet cleanup | Resolved | 2026-03-11 |
| 3. Scratch dependencies | Open | - |
| 4. Pipeline naming | Open | - |
| 5. Conda paths | Resolved | 2026-03-08 |
| 6. SocIsl scripts | Open | - |
| 10. .gitignore blocking files | Resolved | 2026-03-08 |
| 11. Hemoglobin regex bug | Resolved | 2026-03-08 |
| 12. No QC summary tracking | Resolved | 2026-03-08 |
| 13. Stale parent READMEs | Resolved | 2026-03-08 |
| 14. README hardcoded paths | Resolved | 2026-03-10 |
| 15. Harmony batch variable | Resolved | 2026-03-11 |
| 16. Untracked files | Resolved | 2026-03-14 |
| 17. Hardcoded mabdel03 paths | Resolved | 2026-03-14 |
| 18. Archive scripts non-portable | By design | - |

---

## How to Report New Issues

When you encounter a new issue:

1. Add it to this document with:
   - Severity (High/Medium/Low)
   - Affected files
   - Description
   - Recommendation

2. Update the Resolution Tracking table

3. Consider creating a GitHub Issue for tracking

