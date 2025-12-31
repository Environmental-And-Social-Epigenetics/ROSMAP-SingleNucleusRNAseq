# Known Issues and Flagged Problems

This document consolidates known issues identified during the organization of this repository. These issues should be addressed before production use.

## Critical Issues

### 1. Scratch Space Paths Need Configuration

**Severity**: High  
**Affected Scripts**: Preprocessing scripts that read/write intermediate data

The scripts reference old OpenMind scratch paths (`/om/scratch/Mon/mabdel03/`, `/om/scratch/Sun/mabdel03/`) that need to be updated for your cluster.

**Solution**: 
1. Edit `config/paths.sh` and set `SCRATCH_ROOT` to your cluster's scratch filesystem
2. Update individual scripts to use the configuration file, or manually update paths

**Files to update:**

| Directory | File | Lines to Check |
|-----------|------|----------------|
| `Preprocessing/DeJager/02_Cellranger_Counts/` | `Count_DeJager.py` | Lines 10-19 |
| `Preprocessing/DeJager/02_Cellranger_Counts/` | `example_count.sh` | Lines 5-6, 9-10 |
| `Preprocessing/DeJager/03_Cellbender/` | `example_cellbender.sh` | Lines 4-5, 15, 17 |
| `Preprocessing/Tsai/02_Cellranger_Counts/` | All example scripts | Throughout |
| `Preprocessing/Tsai/03_Cellbender/` | All example scripts | Throughout |

**Recommendation**: Create a configuration file with paths as variables.

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

All conda environment paths have been updated to use:
- Init: `source /orcd/data/lhtsai/001/om2/mabdel03/miniforge3/etc/profile.d/conda.sh`
- Envs: `/orcd/data/lhtsai/001/om2/mabdel03/conda_envs/<env_name>`

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

## Resolution Tracking

| Issue | Status | Assigned | Date Fixed |
|-------|--------|----------|------------|
| Hardcoded paths | Open | - | - |
| Demuxlet cleanup | Open | - | - |
| Scratch dependencies | Open | - | - |
| Pipeline naming | Open | - | - |
| Conda paths | Open | - | - |
| SocIsl scripts | Open | - | - |

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

