# Preprocessing Conda Environment Specifications

Environment YAML files for the preprocessing pipeline. Create environments with:

```bash
conda env create -f <spec.yml> -p ${CONDA_ENV_BASE}/<env_name>
```

## Environments

| File | Name | Config Variable | Purpose |
|------|------|----------------|---------|
| `cellbender.yml` | Cellbender_env | `$CELLBENDER_ENV` | Ambient RNA removal (GPU required) |
| `synapse.yml` | synapse_env | `$SYNAPSE_ENV` | DeJager FASTQ download from Synapse |
| `bcftools.yml` | bcftools_env | `$BCFTOOLS_ENV` | BAM/VCF filtering for Demuxlet |
| `globus.yml` | globus_env | `$GLOBUS_ENV` | Globus data transfers between clusters |

## Notes

- **CellBender** requires a GPU. Install PyTorch with CUDA support matching your cluster.
- **Synapse** requires a personal access token from [synapse.org](https://www.synapse.org/).
- **Globus** requires MIT Kerberos authentication via `globus login`.
- **Cell Ranger** (v8.0.0) is distributed as a standalone binary, not via conda. Download from [10x Genomics](https://www.10xgenomics.com/support/software/cell-ranger).
