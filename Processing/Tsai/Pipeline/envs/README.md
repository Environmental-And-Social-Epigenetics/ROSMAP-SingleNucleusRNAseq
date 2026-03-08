# Conda Environment Specifications

Each YAML file specifies the minimum dependencies needed for its pipeline stage.
To recreate an environment:

```bash
conda env create -f stage1_qc.yml
conda env create -f stage2_doublets.yml
conda env create -f stage3_integration.yml
```

After creating, update the corresponding paths in `config/paths.sh` to point to
your new environments.
