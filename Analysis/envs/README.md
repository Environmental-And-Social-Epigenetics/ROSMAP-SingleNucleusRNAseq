# Analysis Conda Environment Specifications

Environment YAML files for downstream analysis pipelines. Create environments with:

```bash
conda env create -f <spec.yml> -p ${CONDA_ENV_BASE}/<env_name>
```

Or use the automated installer with the `--analysis` flag:

```bash
bash setup/install_envs.sh --analysis
```

## Environments

| File | Name | Purpose | Key Packages |
|------|------|---------|-------------|
| `deg.yml` | deg_analysis | Differential expression analysis | DESeq2, edgeR, limma, scanpy |
| `scenic.yml` | scenic_analysis | Regulatory network inference | pySCENIC, loompy |
| `compass.yml` | compass_analysis | Metabolic flux analysis | COMPASS (requires CPLEX) |
| `gsea.yml` | gsea_analysis | Gene set enrichment | WebGestaltR, clusterProfiler |

## External Dependencies

- **SCENIC**: Requires ~3.5 GB of motif ranking reference files from [aertslab](https://resources.aertslab.org/cistarget/)
- **COMPASS**: Requires IBM CPLEX solver ([academic license](https://www.ibm.com/academic/))
