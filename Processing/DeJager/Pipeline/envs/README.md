# Conda Environment Specifications

Each YAML file specifies the minimum dependencies needed for its pipeline stage.

## Automated Setup

The recommended way to create all environments is via the setup script:

```bash
source config/paths.sh
bash setup/install_envs.sh
```

## Manual Setup

To create an environment manually, use the `-p` flag to install it at the
location expected by `config/paths.sh`:

```bash
source config/paths.sh
conda env create -f stage1_qc.yml      -p "${QC_ENV}"
conda env create -f stage2_doublets.yml -p "${SINGLECELL_ENV}"
conda env create -f stage3_integration.yml -p "${BATCHCORR_ENV}"
```

After creating, verify the paths match what `config/paths.sh` defines (or
update `config/paths.local.sh`).

## Environment Mapping

| YAML file | Env name created | `config/paths.sh` variable |
|-----------|------------------|----------------------------|
| `stage1_qc.yml` | `qcEnv` | `QC_ENV` |
| `stage2_doublets.yml` | `single_cell_BP` | `SINGLECELL_ENV` |
| `stage3_integration.yml` | `BatchCorrection_SingleCell` | `BATCHCORR_ENV` |
