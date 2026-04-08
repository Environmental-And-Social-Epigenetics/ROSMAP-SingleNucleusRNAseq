# Analysis

Downstream phenotype analyses built on the processed, annotated snRNA-seq
objects.

## Shared Analysis Contract

The analysis layer expects processed AnnData inputs from `Processing/`.

Required fields depend on cohort:

| Cohort | Required obs fields |
|--------|---------------------|
| Tsai | `cell_type`, `sample_id` or `projid` |
| DeJager | `cell_type`, `patient_id` |

Phenotype tables come from `Data/Phenotypes/`, with ACE analyses using
`ACE_SCORES_CSV` by default.

Generated outputs belong under `ANALYSIS_OUTPUT_ROOT`, not inside the repo tree.

## Current Workflow Status

| Phenotype / analysis | Status | Notes |
|----------------------|--------|-------|
| ACE DEG | `official` | Tsai + DeJager |
| ACE cell-type proportion | `official` | Tsai + DeJager, sccomp |
| SocIsl | `legacy` | older scripts retained for reference |
| Resilient | `scaffold-only` | structure present, not yet formalized |
| TF / SCENIC outside legacy areas | `scaffold-only` | directories exist, but setup and contracts are incomplete |

## Environments

Install the analysis envs with:

```bash
source config/paths.sh
bash setup/install_envs.sh --analysis
```

Official analysis env variables:

- `DEG_ANALYSIS_ENV`
- `NEBULA_ENV`
- `SCCOMP_ENV`
- `SCENIC_ANALYSIS_ENV`
- `COMPASS_ANALYSIS_ENV`
- `GSEA_ANALYSIS_ENV`

Each env now has its own directory with `environment.yml`, `requirements.txt`,
and `README.md`. See [Analysis/envs/README.md](envs/README.md).

## ACE Entry Points

See [Analysis/ACE/README.md](ACE/README.md) for the shared phenotype contract.

Main launchers:

- `Analysis/ACE/DEG/Tsai/aceDegT.sh`
- `Analysis/ACE/DEG/DeJager/aceDegDJ.sh`
- `Analysis/ACE/CellTypeProportion/Tsai/acePropT.sh`
- `Analysis/ACE/CellTypeProportion/DeJager/acePropDJ.sh`

Smoke tests:

- `Analysis/ACE/DEG/Tsai/smoke_test.sh`
- `Analysis/ACE/DEG/DeJager/smoke_test.sh`
- `Analysis/ACE/CellTypeProportion/Tsai/smoke_test.sh`
- `Analysis/ACE/CellTypeProportion/DeJager/smoke_test.sh`

## Output Policy

Analysis code and documentation stay in git. Large generated artifacts do not.

Examples of outputs that should live under `ANALYSIS_OUTPUT_ROOT`:

- split `.h5ad` files
- DEG result tables
- sccomp result objects
- figures
- SLURM logs
