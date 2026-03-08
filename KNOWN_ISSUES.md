# Known Issues and Flagged Problems

This document consolidates known issues identified during the organization of this repository. These issues should be addressed before production use.

## Critical Issues

### 1. Scratch Space Paths Need Configuration

**Severity**: High
**Status**: Partially resolved
**Affected Scripts**: Preprocessing scripts that read/write intermediate data

The **Processing/Tsai/Pipeline/** shell wrappers now source `config/paths.sh` and derive all paths from it.  Preprocessing scripts still reference old OpenMind scratch paths and need updating.

**Remaining files to update:**

| Directory | File | Lines to Check |
|-----------|------|----------------|
| `Preprocessing/DeJager/02_Cellranger_Counts/` | `Count_DeJager.py` | Lines 10-19 |
| `Preprocessing/DeJager/02_Cellranger_Counts/` | `example_count.sh` | Lines 5-6, 9-10 |
| `Preprocessing/DeJager/03_Cellbender/` | `example_cellbender.sh` | Lines 4-5, 15, 17 |
| `Preprocessing/Tsai/02_Cellranger_Counts/` | All example scripts | Throughout |
| `Preprocessing/Tsai/03_Cellbender/` | All example scripts | Throughout |

**Recommendation**: Migrate these to source `config/paths.sh` as done for the processing pipeline.

---

### 2. Demuxlet Script Cleanup Needed

**Severity**: Medium  
**Affected File**: `Preprocessing/DeJager/04_Demuxlet_Freemuxlet/demuxTest.sh`

This script contains extensive commented-out parameter tuning experiments spanning ~120 lines. The active commands are only at:
- Lines 53-54 (pileup generation)
- Line 69 (demuxlet execution)

**Issues:**
- Lines 1-52: Commented setup and early experiments
- Lines 70-151: Commented parameter variations with accuracy notes

**Recommendation**: 
1. Extract active commands to a clean production script
2. Move parameter experiments to a separate documentation file
3. Keep the original as a reference in an `archive/` directory

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
- CellBender (Tsai): 500GB RAM
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

## Resolution Tracking

| Issue | Status | Date Fixed |
|-------|--------|------------|
| 1. Hardcoded paths (pipeline) | Resolved | 2026-03-08 |
| 1. Hardcoded paths (preprocessing) | Open | - |
| 2. Demuxlet cleanup | Open | - |
| 3. Scratch dependencies | Open | - |
| 4. Pipeline naming | Open | - |
| 5. Conda paths | Resolved | 2026-03-08 |
| 6. SocIsl scripts | Open | - |
| 10. .gitignore blocking files | Resolved | 2026-03-08 |
| 11. Hemoglobin regex bug | Resolved | 2026-03-08 |
| 12. No QC summary tracking | Resolved | 2026-03-08 |
| 13. Stale parent READMEs | Resolved | 2026-03-08 |

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

