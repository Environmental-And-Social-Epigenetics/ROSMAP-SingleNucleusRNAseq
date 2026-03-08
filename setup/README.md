# First-Time Setup Guide

This guide walks you through setting up the ROSMAP snRNA-seq pipeline on a new
SLURM cluster, from clone to ready-to-run.

## Prerequisites

- **SLURM** cluster access with GPU nodes (for CellBender)
- **conda** or **mamba** installed (miniconda, miniforge, or anaconda)
- **Git** for cloning the repository

## Step 1: Clone the Repository

```bash
cd /your/workspace
git clone <repo-url> ROSMAP-SingleNucleusRNAseq
```

The pipeline expects large data directories (e.g. `Tsai_Data/`) to live
alongside the repo in the same parent directory (`/your/workspace/`).

## Step 2: Configure Paths

Edit `config/paths.sh` — or better, create `config/paths.local.sh` (gitignored)
to override just the variables that differ on your cluster:

```bash
# config/paths.local.sh — example
export CONDA_INIT_SCRIPT="$HOME/miniforge3/etc/profile.d/conda.sh"
export CONDA_ENV_BASE="$HOME/conda_envs"
export SCRATCH_ROOT="/scratch/$USER"
export DATA_ROOT="/data/project/mylab"
export CELLRANGER_PATH="$HOME/apps/cellranger-8.0.0"
export CELLRANGER_REF="$HOME/references/refdata-gex-GRCh38-2020-A"
export SLURM_MAIL_USER="you@university.edu"
export SLURM_PARTITION="your_partition"
```

Then verify:

```bash
source config/paths.sh
check_paths
```

## Step 3: Install Cell Ranger

Cell Ranger v8.0.0 is required for alignment. It is not available through
conda and must be downloaded from 10x Genomics:

1. Register at https://www.10xgenomics.com/support/software/cell-ranger
2. Download Cell Ranger v8.0.0 and extract it
3. Set `CELLRANGER_PATH` in your config to the extracted directory

## Step 4: Download Reference Genome

The pipeline uses the 10x Genomics human reference `refdata-gex-GRCh38-2020-A`:

```bash
# Download from 10x Genomics (see their website for the current URL)
curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
tar -xzf refdata-gex-GRCh38-2020-A.tar.gz
```

Set `CELLRANGER_REF` in your config to the extracted directory.

## Step 5: Create Conda Environments

Run the automated installer:

```bash
source config/paths.sh
bash setup/install_envs.sh
```

This creates three environments for the Processing pipeline stages:

| Environment | YAML spec | paths.sh variable | Purpose |
|-------------|-----------|-------------------|---------|
| `qcEnv` | `envs/stage1_qc.yml` | `QC_ENV` | Stage 1: QC filtering |
| `single_cell_BP` | `envs/stage2_doublets.yml` | `SINGLECELL_ENV` | Stage 2: Doublet removal |
| `BatchCorrection_SingleCell` | `envs/stage3_integration.yml` | `BATCHCORR_ENV` | Stage 3: Integration |

For CellBender (Preprocessing), create the environment separately since it
requires GPU/CUDA support:

```bash
# Install CellBender in a dedicated environment
conda create -n cellbender_env python=3.10
conda activate cellbender_env
pip install cellbender
```

Set `CELLBENDER_ENV` in your config to the environment path.

## Step 6: Synapse Credentials (DeJager only)

The DeJager dataset is downloaded from Synapse.  If you need to run the DeJager
Preprocessing pipeline:

1. Create an account at https://www.synapse.org
2. Install the Synapse client: `pip install synapseclient`
3. Configure credentials: `synapse login -u <username> -p <password> --rememberMe`

## Step 7: Verify Setup

```bash
source config/paths.sh
check_paths

# Test processing pipeline on a small subset
cd Processing/Tsai/Pipeline
python 01_qc_filter.py --list-samples | head -2
```

## Directory Layout

After setup, your workspace should look like:

```
/your/workspace/
├── ROSMAP-SingleNucleusRNAseq/     # This repository
│   ├── config/paths.sh
│   ├── Preprocessing/
│   ├── Processing/
│   └── Analysis/
└── Tsai_Data/                      # Large data (outside repo)
    ├── FASTQs/
    ├── Cellbender_Outputs/
    └── Processing_Outputs/
        ├── 01_QC_Filtered/
        ├── 02_Doublet_Removed/
        ├── 03_Integrated/
        └── Logs/
```
