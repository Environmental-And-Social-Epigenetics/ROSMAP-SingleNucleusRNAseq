# Configuration

This directory contains configuration files for the ROSMAP snRNA-seq pipeline.

## Files

### `paths.sh`

Central configuration for all paths used by the pipeline.  **This is the single
source of truth** — every shell wrapper should `source` it rather than
hard-coding paths.

**Usage:**
```bash
# Source the configuration
source config/paths.sh

# Use the defined variables
echo $CELLRANGER_REF
echo $DEJAGER_FASTQS

# Initialize conda and activate an environment
init_conda
conda activate $QC_ENV

# Verify all paths
check_paths
```

### `paths.local.sh.template`

Template for creating your local path overrides. Copy and edit:

```bash
cp config/paths.local.sh.template config/paths.local.sh
# Edit config/paths.local.sh with your cluster-specific paths
```

### `paths.local.sh` (gitignored, created by you)

Your local overrides.  Automatically sourced by `paths.sh` if it exists.
Variables set here take precedence over the defaults in `paths.sh`.

## Setup Instructions

1. **Create local overrides** — copy the template and fill in your paths:
   ```bash
   cp config/paths.local.sh.template config/paths.local.sh
   ```
   At minimum, set:
   - `CONDA_INIT_SCRIPT` — path to your conda/mamba init script
   - `CONDA_ENV_BASE` — base directory for conda environments
   - `CELLRANGER_PATH` — directory containing the `cellranger` executable
   - `CELLRANGER_REF` — Cell Ranger reference transcriptome directory
   Optional storage overrides include `TRANSCRIPTOMICS_DATA_ROOT`, `DATA_ROOT`, and `SCRATCH_ROOT`.

2. **Create conda environments** — see `setup/install_envs.sh` or manually:
   ```bash
   conda env create -f envs/processing/stage1_qc/environment.yml -p $CONDA_ENV_BASE/qcEnv
   ```

3. **Verify configuration:**
   ```bash
   source config/paths.sh
   check_paths
   ```

## Key Variables

| Variable | Description |
|----------|-------------|
| `CONDA_INIT_SCRIPT` | Path to conda initialization script |
| `CONDA_ENV_BASE` | Base directory for conda environments |
| `TRANSCRIPTOMICS_DATA_ROOT` | Canonical transcriptomics data namespace |
| `DATA_ROOT` | Permanent storage root, defaulting to the repo root |
| `SCRATCH_ROOT` | Temporary scratch storage, defaulting to `Data/.scratch` |
| `CELLRANGER_REF` | Cell Ranger reference transcriptome |
| `CELLRANGER_PATH` | Cell Ranger installation directory |
| `SLURM_MAIL_USER` | Email for SLURM job notifications |
| `SLURM_PARTITION` | Default SLURM partition (empty = cluster default) |

## Conda Environments

After sourcing `paths.sh`, these environment paths are available:

| Variable | Environment | YAML spec |
|----------|-------------|-----------|
| `QC_ENV` | Stage 1 QC filtering | `envs/processing/stage1_qc/environment.yml` |
| `SINGLECELL_ENV` | Stage 2 doublet removal | `envs/processing/stage2_doublets/environment.yml` |
| `BATCHCORR_ENV` | Stage 3 integration | `envs/processing/stage3_integration/environment.yml` |
| `CELLBENDER_ENV` | CellBender ambient RNA removal | `envs/preprocessing/cellbender/environment.yml` |
| `SYNAPSE_ENV` | Synapse client (DeJager download) | `envs/preprocessing/synapse/environment.yml` |
| `PYTHON_ENV` | General Python | use workflow-specific spec |
