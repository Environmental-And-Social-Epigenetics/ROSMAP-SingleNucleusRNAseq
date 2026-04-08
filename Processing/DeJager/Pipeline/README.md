# DeJager Pipeline

Self-contained three-stage processing for the DeJager ROSMAP snRNA-seq cohort,
from CellBender outputs to an integrated, annotated AnnData object.

## Canonical Defaults

- patient identities come from `DEJAGER_PATIENT_MAP` plus `patient_id_overrides.json`
- the canonical Stage 3 Harmony key is **`library_id`**
- the shared marker reference is taken from `DEJAGER_MARKERS_RDS`
- ORA annotation runs on the raw gene space (`use_raw=True`)

## Inputs

- `${DEJAGER_PREPROCESSED}/{library}/processed_feature_bc_matrix_filtered.h5`
- `${DEJAGER_PATIENT_MAP}`
- `${DEJAGER_PATIENT_ID_OVERRIDES}`
- `${DEJAGER_MARKERS_RDS}`

## Outputs

```text
${DEJAGER_PROCESSING_OUTPUTS}/
  01_QC_Filtered/
  02_Doublet_Removed/
  03_Integrated/                 canonical library_id integration
  03_Integrated_patient_id/      sensitivity analysis
  03_Integrated_pool_batch/      sensitivity analysis
  03_Integrated_derived_batch/   sensitivity analysis
  03_Evaluation/                 comparison outputs
  Logs/
```

## Stage Summary

### Stage 1

- script: `01_qc_filter.py`
- wrapper: `01_qc_filter.sh`
- env: `envs/stage1_qc/`
- output: `{library}_qc.h5ad`

### Stage 2

- script: `02_doublet_removal.Rscript`
- wrapper: `02_doublet_removal.sh`
- env: `envs/stage2_doublets/`
- output: `{library}_singlets.h5ad`

### Stage 3

- canonical wrapper: `03_integration_annotation.sh`
- env: `envs/stage3_integration/`
- canonical output: `${DEJAGER_INTEGRATED}`

Sensitivity wrappers:

- `03_integration_annotation_patient_id.sh`
- `03_integration_annotation_pool_batch.sh`
- `03_integration_annotation_derived_batch.sh`

## Quick Start

```bash
source config/paths.sh
bash setup/install_envs.sh
bash config/preflight.sh dejager-stage1

cd Processing/DeJager/Pipeline
./submit_pipeline.sh all
```

Submit selected stages:

```bash
./submit_pipeline.sh 1
./submit_pipeline.sh 2 3
./submit_pipeline.sh 3 3b 3c 3d 3e
```

Stage meanings:

- `3`: canonical `library_id`
- `3b`: `patient_id`
- `3c`: `pool_batch`
- `3d`: `derived_batch`
- `3e`: comparison report

## Direct Test Runs

```bash
python 01_qc_filter.py --sample-ids LIB5001_R1234567-alone
Rscript 02_doublet_removal.Rscript --sample-ids LIB5001_R1234567-alone
python 03_integration_annotation.py --sample-ids LIB5001_R1234567-alone
```

The direct Stage 3 path resolves marker inputs through `DEJAGER_MARKERS_RDS`, so
it no longer depends on an untracked local marker symlink.

## Running the Full Pipeline

**Status (April 2026):** The DeJager processing pipeline has not yet been
executed end-to-end. CellBender outputs (131 library directories) and CellRanger
outputs (47 libraries) are on Engaging and ready to process.

### Prerequisites

1. **Patient assignment CSV** (`DEJAGER_PATIENT_MAP`):
   - File: `cell_to_patient_assignmentsFinal1.csv` (~155 MB)
   - Location: `/home/nkhera/orcd/pool/WGS/` on Engaging
   - This file maps cell barcodes to patient IDs from Demuxlet/Freemuxlet
   - Set in `config/paths.local.sh`:
     ```bash
     export DEJAGER_PATIENT_MAP="/home/nkhera/orcd/pool/WGS/cell_to_patient_assignmentsFinal1.csv"
     ```

2. **CellBender outputs** (`DEJAGER_PREPROCESSED`):
   - 131 directories containing `processed_feature_bc_matrix_filtered.h5`
   - Set in `config/paths.local.sh` to point to actual data location

3. **Conda environments**:
   ```bash
   source config/paths.sh
   bash setup/install_envs.sh --processing
   ```

4. **Preflight check**:
   ```bash
   bash config/preflight.sh dejager-stage1
   bash config/preflight.sh dejager-stage2
   bash config/preflight.sh dejager-stage3
   ```

### Resource Requirements

| Stage | Memory | Cores | Time (est.) |
|-------|--------|-------|-------------|
| Stage 1 (QC) | 50 GB | 1 per sample | 1-2 hours |
| Stage 2 (Doublet) | 100 GB | 4 per sample | 2-4 hours |
| Stage 3 (Integration) | **500 GB** | 30 | 6-12 hours |

### Expected Outputs

After a successful run, verify:
- `${DEJAGER_QC_FILTERED}/` contains `{library}_qc.h5ad` files + `qc_summary.csv`
- `${DEJAGER_DOUBLET_REMOVED}/` contains `{library}_singlets.h5ad` files + `doublet_summary.csv`
- `${DEJAGER_INTEGRATED}/` contains `dejager_integrated.h5ad` and `dejager_annotated.h5ad`
- `dejager_annotated.h5ad` should have `cell_type` in `.obs` from ORA annotation

### Execution

```bash
cd Processing/DeJager/Pipeline
./submit_pipeline.sh all
```

Or stage by stage:
```bash
./submit_pipeline.sh 1    # QC filtering
# Wait for Stage 1 to complete, then:
./submit_pipeline.sh 2    # Doublet removal
# Wait for Stage 2 to complete, then:
./submit_pipeline.sh 3    # Integration (canonical library_id)
```

## Patient Map

`DEJAGER_PATIENT_MAP` is too large for git and must be provided locally. See
[setup/README.md](../../../setup/README.md) and `config/paths.local.sh.template`
for the expected override pattern.
