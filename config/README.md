# Configuration

This directory contains configuration files for the ROSMAP snRNA-seq pipeline.

## Files

### `paths.sh`

Central configuration for all paths used by the pipeline.

**Usage:**
```bash
# Source the configuration
source config/paths.sh

# Use the defined variables
echo $CELLRANGER_REF
echo $DEJAGER_FASTQS

# Initialize conda
init_conda
conda activate $CELLBENDER_ENV

# Verify all paths
check_paths
```

**Key Variables:**

| Variable | Description |
|----------|-------------|
| `CONDA_INIT_SCRIPT` | Path to conda initialization script |
| `CONDA_ENV_BASE` | Base directory for conda environments |
| `DATA_ROOT` | Permanent storage for processed data |
| `SCRATCH_ROOT` | Temporary scratch storage (UPDATE THIS!) |
| `CELLRANGER_REF` | Cell Ranger reference transcriptome |

## Setup Instructions

1. **Update scratch path:**
   
   Edit `paths.sh` and set `SCRATCH_ROOT` to your cluster's scratch filesystem:
   ```bash
   export SCRATCH_ROOT="/your/scratch/path"
   ```

2. **Verify configuration:**
   ```bash
   source config/paths.sh
   check_paths
   ```

3. **Update individual scripts (optional):**
   
   To use this configuration in scripts, add at the top:
   ```bash
   SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
   REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
   source "$REPO_ROOT/config/paths.sh"
   ```

## Environment Variables

After sourcing `paths.sh`, these conda environments are available:

| Variable | Environment |
|----------|-------------|
| `CELLBENDER_ENV` | CellBender for ambient RNA removal |
| `SYNAPSE_ENV` | Synapse client for data download |
| `PYTHON_ENV` | General Python environment |
| `SINGLECELL_ENV` | Single-cell analysis (scanpy) |
| `BATCHCORR_ENV` | Batch correction (Harmony) |
| `QC_ENV` | Quality control |

