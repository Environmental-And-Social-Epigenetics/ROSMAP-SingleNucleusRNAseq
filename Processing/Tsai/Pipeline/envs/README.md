# Processing Environment Specifications

The processing pipeline uses three shared environments. Each environment now
has its own directory containing:

- `environment.yml` — canonical conda spec
- `requirements.txt` — companion Python package list
- `README.md` — exact install instructions

Use the automated installer:

```bash
source config/paths.sh
bash setup/install_envs.sh
```

## Environments

| Directory | Env name | Config variable | Pattern | Purpose |
|-----------|----------|-----------------|---------|---------|
| `stage1_qc/` | `qcEnv` | `QC_ENV` | requirements-complete | Stage 1 percentile-based QC filtering |
| `stage2_doublets/` | `single_cell_BP` | `SINGLECELL_ENV` | hybrid | Stage 2 scDblFinder doublet removal |
| `stage3_integration/` | `BatchCorrection_SingleCell` | `BATCHCORR_ENV` | hybrid | Stage 3 integration, Harmony, and ORA annotation |

Manual setup is documented in each per-env README.
