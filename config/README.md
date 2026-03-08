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

### `paths.local.sh` (optional, gitignored)

Create this file to override any variable from `paths.sh` without modifying the
tracked file.  It is automatically sourced at the end of `paths.sh` if it
exists.

```bash
# Example paths.local.sh
export CONDA_INIT_SCRIPT="$HOME/miniforge3/etc/profile.d/conda.sh"
export SCRATCH_ROOT="/scratch/$USER"
export QC_ENV="$HOME/conda_envs/qcEnv"
```

## Setup Instructions

1. **Update paths** — edit `paths.sh` (or create `paths.local.sh`) and set at
   minimum:
   - `CONDA_INIT_SCRIPT` — path to your conda/mamba init script
   - `SCRATCH_ROOT` — your cluster's scratch filesystem
   - `DATA_ROOT` — permanent storage location
   - `CELLRANGER_PATH` — directory containing the `cellranger` executable
   - `CELLRANGER_REF` — Cell Ranger reference transcriptome directory

2. **Create conda environments** — see `setup/install_envs.sh` or manually:
   ```bash
   conda env create -f Processing/Tsai/Pipeline/envs/stage1_qc.yml -p $CONDA_ENV_BASE/qcEnv
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
| `DATA_ROOT` | Permanent storage for processed data |
| `SCRATCH_ROOT` | Temporary scratch storage |
| `CELLRANGER_REF` | Cell Ranger reference transcriptome |
| `CELLRANGER_PATH` | Cell Ranger installation directory |
| `SLURM_MAIL_USER` | Email for SLURM job notifications |
| `SLURM_PARTITION` | Default SLURM partition (empty = cluster default) |

## Conda Environments

After sourcing `paths.sh`, these environment paths are available:

| Variable | Environment | YAML spec |
|----------|-------------|-----------|
| `QC_ENV` | Stage 1 QC filtering | `Processing/Tsai/Pipeline/envs/stage1_qc.yml` |
| `SINGLECELL_ENV` | Stage 2 doublet removal | `Processing/Tsai/Pipeline/envs/stage2_doublets.yml` |
| `BATCHCORR_ENV` | Stage 3 integration | `Processing/Tsai/Pipeline/envs/stage3_integration.yml` |
| `CELLBENDER_ENV` | CellBender ambient RNA removal | — |
| `SYNAPSE_ENV` | Synapse client (DeJager download) | — |
| `PYTHON_ENV` | General Python | — |
