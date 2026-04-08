# Preprocessing Environment Specifications

Each supported preprocessing environment now lives in its own directory with:

- `environment.yml` — canonical conda spec
- `requirements.txt` — companion Python package list
- `README.md` — exact setup instructions

Use the automated installer:

```bash
source config/paths.sh
bash setup/install_envs.sh --preprocessing
```

Manual setup is also supported. Each env README documents two paths:

- `conda`: `conda env create -f environment.yml -p <env-path>`
- `requirements`: create a fresh env, install any required conda/system packages, then `pip install -r requirements.txt`

## Environments

| Directory | Env name | Config variable | Pattern | Purpose |
|-----------|----------|-----------------|---------|---------|
| `cellbender/` | `Cellbender_env` | `CELLBENDER_ENV` | hybrid | Ambient RNA removal with CUDA-enabled PyTorch |
| `synapse/` | `synapse_env` | `SYNAPSE_ENV` | requirements-complete | DeJager FASTQ downloads from Synapse |
| `bcftools/` | `bcftools_env` | `BCFTOOLS_ENV` | hybrid | BAM/VCF filtering for Demuxlet |
| `globus/` | `globus_env` | `GLOBUS_ENV` | requirements-complete | Globus data transfers between clusters |

## Notes

- CellBender still requires a GPU and a matching CUDA runtime.
- Synapse requires a personal access token from [synapse.org](https://www.synapse.org/).
- Globus requires interactive `globus login`.
- Cell Ranger is distributed separately from conda; keep using `CELLRANGER_PATH`.
