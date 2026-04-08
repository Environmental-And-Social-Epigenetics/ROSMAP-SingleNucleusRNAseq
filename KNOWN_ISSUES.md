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
**Status**: Resolved

Scripts previously referenced temporary scratch directories on Openmind:
- `/om/scratch/Mon/mabdel03/` (weekly cleanup)
- `/om/scratch/Sun/mabdel03/` (weekly cleanup)

**Update (March 2026)**: Openmind was decommissioned. All data has been transferred
to Engaging, which uses a different storage tier model. Users should configure
`SCRATCH_ROOT` in `paths.local.sh` to an appropriate Engaging scratch location.

**Update (March 2026 audit)**: `config/paths.local.sh` created with `SCRATCH_ROOT`
pointing to `/orcd/data/lhtsai/001/mabdel03/ROSMAP_Data`. All pipeline path variables
now resolve correctly to actual data in `ROSMAP_Data/Single_Nucleus/`.

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

### 6. Missing SocIsl Processing Scripts — PARTIALLY RESOLVED

**Severity**: Medium
**Location**: `Preprocessing/Tsai/`, `Analysis/SocIsl/`
**Status**: Partially resolved

SocIsl cohort has batch scripts but no dedicated notebooks for Cell Ranger/CellBender
script generation in `Preprocessing/`.

**Update (March 2026)**: Legacy SocIsl analysis scripts (DEG, SCENIC, COMPASS, GSEA)
from the Openmind working directory have been copied into `Analysis/SocIsl/`. These
scripts contain hardcoded Openmind paths and will need updates to run on Engaging.
The preprocessing gap for SocIsl remains open — the current unified pipeline
(`Processing/Tsai/Pipeline/submit_pipeline.sh`) handles all Tsai patients including
SocIsl cohort members, so dedicated SocIsl preprocessing scripts are not needed.

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

### 19. DeJager Processing Pipeline Not Yet Run

**Severity**: High
**Status**: Open

The DeJager processing pipeline (QC filtering, doublet removal, integration &
annotation) has never been executed. CellBender outputs (131 directories) and
CellRanger outputs (47 libraries) are on Engaging and ready to process.

**Recommendation**: Run the pipeline on Engaging:
```bash
cd Processing/DeJager/Pipeline && ./submit_pipeline.sh all
```

Prerequisites:
- Configure `DEJAGER_PREPROCESSED` to point to the CellBender outputs on Engaging
- Obtain `cell_to_patient_assignmentsFinal1.csv` from `/home/nkhera/orcd/pool/WGS/`
- Create conda environments via `setup/install_envs.sh`

---

### 20. No Environment Specs for Preprocessing — RESOLVED

**Severity**: Medium
**Status**: Resolved

The preprocessing steps (CellBender, Synapse download, bcftools, Globus) had no
conda YAML specifications in the repo, making it impossible to recreate the
environments from scratch.

**Fix**: Created `Preprocessing/envs/` with YAML specs for all four preprocessing
environments (cellbender.yml, synapse.yml, bcftools.yml, globus.yml).

---

### 21. No Environment Specs for Analysis — RESOLVED

**Severity**: Medium
**Status**: Resolved

The Analysis directory had no conda YAML specifications for downstream analysis
(DEG, SCENIC, COMPASS, GSEA).

**Fix**: Created `Analysis/envs/` with YAML specs for all four analysis environments.

---

### 22. Engaging Data Layout Undocumented — RESOLVED

**Severity**: Medium
**Status**: Resolved

After the Openmind-to-Engaging migration, the data layout on Engaging was not
documented anywhere in the repo.

**Fix**: Added comprehensive "Data on MIT Engaging (March 2026 Migration)" section
to `Data_Access/README.md` with directory trees, file counts, and Globus audit trail.

---

### 23. paths.local.sh Missing — Pipeline Cannot Find Data

**Severity**: Critical
**Status**: Resolved

`config/paths.local.sh` was never created after the Openmind-to-Engaging migration.
All path variables defaulted to nonexistent locations (`__UNCONFIGURED__` sentinels
or `${WORKSPACE_ROOT}/Tsai_Data/...` which does not exist on Engaging). The in-repo
`Data/Transcriptomics/` directories are empty placeholders. The 4.8 TB of actual
data at `/orcd/data/lhtsai/001/mabdel03/ROSMAP_Data/Single_Nucleus/` was unreachable
by the pipeline.

**Fix**: Created `config/paths.local.sh` mapping all path variables to the correct
Engaging locations. Verified with `check_paths` (PASS) and `preflight.sh tsai-stage1`
(all green).

---

### 24. copy_data.sbatch Hardcoded Openmind Paths

**Severity**: Medium
**Status**: Resolved

`Data/Transcriptomics/copy_data.sbatch` had hardcoded `REPO_ROOT` pointing to
`/om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq` and all 6 source paths
pointed to decommissioned Openmind scratch directories.

**Fix**: Rewrote to auto-detect `REPO_ROOT` from `BASH_SOURCE`, source
`config/paths.sh`, and use config variables for all source paths.

---

### 25. Conda Environment Name Mismatches on Engaging

**Severity**: Medium
**Status**: Resolved

Two conda environment names in `config/paths.sh` defaults do not match the actual
environment names on Engaging:
- `SINGLECELL_ENV`: default `single_cell_BP` vs actual `single_cell_BP4`
- `BATCHCORR_ENV`: default `BatchCorrection_SingleCell` vs actual `batchCorrectionEnv`

**Fix**: Overridden in `config/paths.local.sh`. Additionally, `batchCorrectionEnv` is
missing packages needed for Stage 3 (scanpy, harmonypy, decoupler, rpy2) and
`single_cell_BP4` is missing `zellkonverter` for Stage 2.

---

### 26. Analysis/SocIsl Scripts Have Hardcoded Openmind Paths

**Severity**: Medium
**Status**: Resolved

All ~60 production scripts in `Analysis/SocIsl/` (DEG, SCENIC, COMPASS, GSEA,
TF) have been migrated to use `config/paths.sh`:

- Shell wrappers: source `config/paths.sh`, use `init_conda`, `activate_env`,
  and config variables (`$SOCISL_OUTPUT_ROOT`, `$NEBULA_ENV`, `$CPLEX_DIR`, etc.)
- R scripts: use `Sys.getenv("SOCISL_OUTPUT_ROOT")` for working directories
- Python scripts: use `os.environ["SOCISL_OUTPUT_ROOT"]` for paths
- SBATCH email headers: replaced with `__SET_YOUR_EMAIL__` placeholder

The `_data_prep/` scripts are marked as DEPRECATED with header comments and a
README. They retain original Openmind paths as they are historical reference
only (the Processing pipeline supersedes them).

New config variables added to `config/paths.sh`:
- `SOCISL_OUTPUT_ROOT`, `ACE_OUTPUT_ROOT`, `RESILIENT_OUTPUT_ROOT`
- `CPLEX_DIR`, `SCENIC_RANKING_DIR`, `SCENIC_TF_LIST`

---

### 27. Conda Environments Missing Packages for Stages 2-3

**Severity**: Medium
**Status**: Resolved

Preflight checks revealed missing packages in existing legacy conda environments
(`batchCorrectionEnv`, `single_cell_BP4`). These are old environments from the
Openmind era that do not match the current pipeline specs.

**Resolution**: The canonical environment specifications in
`Processing/Tsai/Pipeline/envs/` (and the identical DeJager copies) already
contain all required packages. The fix is to rebuild environments from the
official specs rather than patching the legacy envs:

```bash
source config/paths.sh
bash setup/install_envs.sh --processing --method=conda
```

This creates fresh environments at `${CONDA_ENV_BASE}/` with the correct names
and all required packages. After rebuilding, update `config/paths.local.sh` to
point `SINGLECELL_ENV` and `BATCHCORR_ENV` at the new environments, then
validate with `bash config/preflight.sh tsai-stage2 && bash config/preflight.sh tsai-stage3`.

---

## Resolution Tracking

| Issue | Status | Date Fixed |
|-------|--------|------------|
| 1. Hardcoded paths (pipeline) | Resolved | 2026-03-08 |
| 1. Hardcoded paths (preprocessing) | Resolved | 2026-03-08 |
| 2. Demuxlet cleanup | Resolved | 2026-03-11 |
| 3. Scratch dependencies | Resolved | 2026-03-18 |
| 4. Pipeline naming | Open | - |
| 5. Conda paths | Resolved | 2026-03-08 |
| 6. SocIsl scripts | Resolved | 2026-04-08 |
| 10. .gitignore blocking files | Resolved | 2026-03-08 |
| 11. Hemoglobin regex bug | Resolved | 2026-03-08 |
| 12. No QC summary tracking | Resolved | 2026-03-08 |
| 13. Stale parent READMEs | Resolved | 2026-03-08 |
| 14. README hardcoded paths | Resolved | 2026-03-10 |
| 15. Harmony batch variable | Resolved | 2026-03-11 |
| 16. Untracked files | Resolved | 2026-03-14 |
| 17. Hardcoded mabdel03 paths | Resolved | 2026-03-14 |
| 18. Archive scripts non-portable | By design | - |
| 19. DeJager processing not run | Documented | 2026-04-08 |
| 20. No preprocessing env specs | Resolved | 2026-03-16 |
| 21. No analysis env specs | Resolved | 2026-03-16 |
| 22. Engaging data undocumented | Resolved | 2026-03-16 |
| 23. paths.local.sh missing | Resolved | 2026-03-18 |
| 24. copy_data.sbatch Openmind paths | Resolved | 2026-03-18 |
| 25. Conda env name mismatches | Resolved | 2026-03-18 |
| 26. Analysis/SocIsl hardcoded paths | Resolved | 2026-04-08 |
| 27. Conda envs missing packages | Resolved | 2026-04-08 |

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

