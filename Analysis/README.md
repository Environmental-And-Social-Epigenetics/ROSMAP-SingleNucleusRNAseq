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

Status vocabulary (consistent across all READMEs in this repo):

- `production` — fully implemented and run for the applicable cohort(s); trustworthy
- `implemented` — code complete and runnable but not yet validated/run end-to-end
- `scaffold` — directory/structure exists but the analysis is not implemented
- `migrated` — ported from legacy (e.g. SocIsl) with paths updated, not re-validated

| Phenotype / analysis | Status | Notes |
|----------------------|--------|-------|
| ACE DEG | `production` | Tsai + DeJager |
| ACE cell-type proportion | `production` | Tsai + DeJager, sccomp |
| SocIsl | `migrated` | older scripts ported from legacy, not re-validated |
| Resilient | `scaffold` | structure present, not yet formalized |
| TF / SCENIC outside migrated areas | `scaffold` | directories exist, but setup and contracts are incomplete |

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
and `README.md`. See [envs/README.md](../envs/README.md) (analysis specs live
under `envs/analysis/`).

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
