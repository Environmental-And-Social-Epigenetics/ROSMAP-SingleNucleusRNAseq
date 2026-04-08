# First-Time Setup Guide

This guide covers the shared setup contract for the ROSMAP transcriptomics repo.

## 1. Clone And Configure Paths

```bash
cd /your/workspace
git clone <repo-url> Transcriptomics
cd Transcriptomics
cp config/paths.local.sh.template config/paths.local.sh
```

Set at minimum:

- `CONDA_INIT_SCRIPT`
- `CONDA_ENV_BASE`
- `DATA_ROOT`
- `SCRATCH_ROOT`
- `CELLRANGER_PATH`
- `CELLRANGER_REF`

Optional but commonly needed:

- `ANALYSIS_OUTPUT_ROOT`
- `DEJAGER_PATIENT_MAP`
- `DEJAGER_SYNAPSE_FASTQ_CSV`
- `NEBULA_ENV`
- `SCCOMP_ENV`

Then validate:

```bash
source config/paths.sh
check_paths
```

## 2. Install Environments

Official environments are organized as per-env directories. Each one contains:

- `environment.yml`
- `requirements.txt`
- `README.md`

Use the installer:

```bash
source config/paths.sh
bash setup/install_envs.sh                  # processing envs
bash setup/install_envs.sh --preprocessing  # preprocessing envs too
bash setup/install_envs.sh --analysis       # analysis envs too
bash setup/install_envs.sh --all            # everything
bash setup/install_envs.sh --analysis --recreate
```

Choose an install method explicitly when needed:

```bash
bash setup/install_envs.sh --all --method=conda
bash setup/install_envs.sh --all --method=requirements
```

### Install Methods

- `conda`: create directly from `environment.yml`
- `requirements`: create the bootstrap env described in each env `README.md`, then install the Python layer from `requirements.txt`

Every official env has a companion `requirements.txt`. Some are
`requirements-complete`, while others are `hybrid` and still need conda or
system packages first.

The installer now validates each env after creation, and it also validates
already-existing envs instead of assuming they are healthy. Use `--recreate`
when you need to rebuild a stale env from the tracked specs.

### Official Environment Inventory

| Area | Env ID | Variable | Type |
|------|--------|----------|------|
| Processing | `stage1_qc` | `QC_ENV` | requirements-complete |
| Processing | `stage2_doublets` | `SINGLECELL_ENV` | hybrid |
| Processing | `stage3_integration` | `BATCHCORR_ENV` | hybrid |
| Preprocessing | `cellbender` | `CELLBENDER_ENV` | hybrid |
| Preprocessing | `synapse` | `SYNAPSE_ENV` | requirements-complete |
| Preprocessing | `bcftools` | `BCFTOOLS_ENV` | hybrid |
| Preprocessing | `globus` | `GLOBUS_ENV` | requirements-complete |
| Analysis | `deg` | `DEG_ANALYSIS_ENV` | hybrid |
| Analysis | `nebula` | `NEBULA_ENV` | hybrid |
| Analysis | `sccomp` | `SCCOMP_ENV` | hybrid |
| Analysis | `scenic` | `SCENIC_ANALYSIS_ENV` | requirements-complete |
| Analysis | `compass` | `COMPASS_ANALYSIS_ENV` | requirements-complete |
| Analysis | `gsea` | `GSEA_ANALYSIS_ENV` | hybrid |

See the env READMEs in:

- [Preprocessing/envs/README.md](../Preprocessing/envs/README.md)
- [Processing/Tsai/Pipeline/envs/README.md](../Processing/Tsai/Pipeline/envs/README.md)
- [Analysis/envs/README.md](../Analysis/envs/README.md)

## 3. External Tools

Not everything is installed through the repo envs.

### Cell Ranger

Cell Ranger v8.0.0 must be installed separately from 10x Genomics and exposed
through `CELLRANGER_PATH`.

### Reference Transcriptome

Set `CELLRANGER_REF` to the extracted `refdata-gex-GRCh38-2020-A` directory.

### Singularity

Needed for DeJager Demuxlet/Freemuxlet helper workflows. Configure
`SINGULARITY_MODULE` if your cluster uses a different module name.

## 4. External Data Contracts

The repo tracks phenotype tables and lightweight metadata, but several large or
protected inputs stay external:

- DeJager Synapse FASTQs
- DeJager barcode-to-patient map
- DeJager WGS / Demuxlet inputs
- cohort-specific raw FASTQs and Cell Ranger outputs

Important DeJager-specific variables:

- `DEJAGER_PATIENT_MAP`
- `DEJAGER_SYNAPSE_FASTQ_CSV`
- `DEJAGER_WGS_DIR`

See [Processing/DeJager/Pipeline/README.md](../Processing/DeJager/Pipeline/README.md)
for the patient-map details.

## 5. Validate The Setup

```bash
source config/paths.sh
check_paths
bash config/preflight.sh env-specs
bash config/preflight.sh ace-tsai-smoke
```

Useful spot checks:

```bash
cd Processing/Tsai/Pipeline
python 01_qc_filter.py --list-samples | head

cd ../../Analysis/ACE/DEG/Tsai
bash smoke_test.sh

cd ../../CellTypeProportion/Tsai
bash smoke_test.sh
```

## 6. Output Locations

Generated outputs should not be written into the repo tree.

- processing outputs: the configured `*_PROCESSING_OUTPUTS` directories
- downstream analysis outputs: `ANALYSIS_OUTPUT_ROOT`

If you leave `ANALYSIS_OUTPUT_ROOT` unset, it defaults to a sibling
`Analysis_Outputs/` directory next to the repo.
