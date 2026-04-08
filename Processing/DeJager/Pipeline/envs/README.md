# DeJager Processing Environment Specifications

The DeJager processing pipeline uses the same three stage environments as the
Tsai pipeline. Each stage now has a dedicated directory containing:

- `environment.yml` — canonical conda spec
- `requirements.txt` — companion Python package list
- `README.md` — exact install instructions

## Environments

| Directory | Env name | Config variable | Pattern | Purpose |
|-----------|----------|-----------------|---------|---------|
| `stage1_qc/` | `qcEnv` | `QC_ENV` | requirements-complete | Stage 1 percentile-based QC filtering |
| `stage2_doublets/` | `single_cell_BP` | `SINGLECELL_ENV` | hybrid | Stage 2 scDblFinder doublet removal |
| `stage3_integration/` | `BatchCorrection_SingleCell` | `BATCHCORR_ENV` | hybrid | Stage 3 integration, Harmony, and ORA annotation |

Use `bash setup/install_envs.sh` for automated creation, or follow the per-env
READMEs in this directory for manual `conda` or `requirements` setup.
