# Analysis Environment Specifications

Each supported analysis environment now has its own directory containing:

- `environment.yml` — canonical conda spec
- `requirements.txt` — companion Python package list
- `README.md` — exact install instructions

Use the automated installer:

```bash
source config/paths.sh
bash setup/install_envs.sh --analysis
```

## Environments

| Directory | Env name | Config variable | Pattern | Purpose |
|-----------|----------|-----------------|---------|---------|
| `deg/` | `deg_analysis` | `DEG_ANALYSIS_ENV` | hybrid | General differential expression analysis |
| `scenic/` | `scenic_analysis` | `SCENIC_ANALYSIS_ENV` | requirements-complete | pySCENIC regulatory network inference |
| `compass/` | `compass_analysis` | `COMPASS_ANALYSIS_ENV` | requirements-complete | COMPASS metabolic modeling |
| `gsea/` | `gsea_analysis` | `GSEA_ANALYSIS_ENV` | hybrid | Gene-set enrichment in R |
| `nebula/` | `nebulaAnalysis7` | `NEBULA_ENV` | hybrid | ACE pseudobulk DEG and H5AD preprocessing |
| `sccomp/` | `sccompAnalysis` | `SCCOMP_ENV` | hybrid | ACE cell-type proportion modeling |

## Notes

- SCENIC still requires external motif ranking reference files.
- COMPASS still requires IBM CPLEX outside the conda environment.
- The ACE envs are now first-class supported environments in setup and docs.
