# First-Time Setup Guide

This guide walks you through setting up the ROSMAP snRNA-seq pipeline on a new
SLURM cluster, from clone to ready-to-run.

## Prerequisites

- **SLURM** cluster access with GPU nodes (for CellBender)
- **conda** or **mamba** installed (miniconda, miniforge, or anaconda)
- **Git** for cloning the repository
- **Singularity** (for DeJager Demuxlet step only): `module load openmind/singularity/3.10.4`
- **Globus CLI** (for Tsai FASTQ transfer only): install in a dedicated conda environment

## Step 1: Clone the Repository

```bash
cd /your/workspace
git clone <repo-url> ROSMAP-SingleNucleusRNAseq
```

The pipeline expects large data directories (e.g. `Tsai_Data/`) to live
alongside the repo in the same parent directory (`/your/workspace/`).

## Step 2: Configure Paths

Copy the template and fill in your cluster-specific paths:

```bash
cp config/paths.local.sh.template config/paths.local.sh
# Edit config/paths.local.sh — set at minimum:
#   CONDA_INIT_SCRIPT, CONDA_ENV_BASE, DATA_ROOT, SCRATCH_ROOT,
#   CELLRANGER_PATH, CELLRANGER_REF, SLURM_MAIL_USER
```

Then verify:

```bash
source config/paths.sh
check_paths
```

See `config/paths.local.sh.template` for MIT Engaging example paths.

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
bash setup/install_envs.sh              # Processing environments only (required)
bash setup/install_envs.sh --analysis   # Also install Analysis environments
bash setup/install_envs.sh --preprocessing  # Also install Preprocessing environments
bash setup/install_envs.sh --all        # Install everything
```

This creates environments from YAML specs in the repo:

**Processing (always installed):**

| Environment | YAML spec | paths.sh variable | Purpose |
|-------------|-----------|-------------------|---------|
| `qcEnv` | `Processing/.../envs/stage1_qc.yml` | `QC_ENV` | Stage 1: QC filtering |
| `single_cell_BP` | `Processing/.../envs/stage2_doublets.yml` | `SINGLECELL_ENV` | Stage 2: Doublet removal |
| `BatchCorrection_SingleCell` | `Processing/.../envs/stage3_integration.yml` | `BATCHCORR_ENV` | Stage 3: Integration |

**Preprocessing (`--preprocessing`):**

| Environment | YAML spec | paths.sh variable | Purpose |
|-------------|-----------|-------------------|---------|
| `Cellbender_env` | `Preprocessing/envs/cellbender.yml` | `CELLBENDER_ENV` | Ambient RNA removal (GPU) |
| `synapse_env` | `Preprocessing/envs/synapse.yml` | `SYNAPSE_ENV` | Synapse downloads |
| `bcftools_env` | `Preprocessing/envs/bcftools.yml` | `BCFTOOLS_ENV` | BAM/VCF filtering |
| `globus_env` | `Preprocessing/envs/globus.yml` | `GLOBUS_ENV` | Data transfers |

**Analysis (`--analysis`):**

| Environment | YAML spec | Purpose |
|-------------|-----------|---------|
| `deg_analysis` | `Analysis/envs/deg.yml` | DEG (DESeq2, limma, edgeR) |
| `scenic_analysis` | `Analysis/envs/scenic.yml` | pySCENIC regulatory networks |
| `compass_analysis` | `Analysis/envs/compass.yml` | COMPASS metabolic analysis |
| `gsea_analysis` | `Analysis/envs/gsea.yml` | GSEA/pathway analysis |

## Step 6: DeJager Patient Map (DeJager only)

The barcode-to-patient mapping file (`cell_to_patient_assignmentsFinal0.csv`,
~155 MB) is required for Stage 1 of the DeJager Processing pipeline but is too
large for git. See
[Processing/DeJager/Pipeline/README.md](../Processing/DeJager/Pipeline/README.md#obtaining-the-patient-map)
for instructions on obtaining it.

Set `DEJAGER_PATIENT_MAP` in `config/paths.local.sh` to point to the file.

## Step 7: Synapse Credentials (DeJager only)

The DeJager dataset is downloaded from Synapse.  If you need to run the DeJager
Preprocessing pipeline:

1. Create an account at https://www.synapse.org
2. Install the Synapse client: `pip install synapseclient`
3. Configure credentials: `synapse login -u <username> -p <password> --rememberMe`

## Step 8: Verify Setup

```bash
source config/paths.sh
check_paths

# Test processing pipeline on a small subset
cd Processing/Tsai/Pipeline
python 01_qc_filter.py --list-samples | head -2
```

## SLURM Partitions

Configure your cluster's partitions in `config/paths.local.sh`:

```bash
export SLURM_PARTITION="your_default_partition"
export SLURM_PARTITION_GPU="your_gpu_partition"
```

If no partition is set, SLURM will use the cluster default.

### MIT Engaging / Openmind Reference

| Pipeline Use | Partition | Notes |
|-------------|-----------|-------|
| Cell Ranger | `mit_preemptable` | Lower-priority, shorter queue; jobs may be preempted |
| CellBender | `mit_normal_gpu` | GPU partition (A100) for ambient RNA removal |
| Stage 3 Integration | `lhtsai` | Lab-specific high-memory partition (500GB) |

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
