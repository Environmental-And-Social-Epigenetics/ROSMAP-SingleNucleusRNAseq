# SocIsl Cell-Type Proportion Analysis

Bayesian compositional analysis of cell type proportions across social isolation
phenotype levels, using sccomp.

## Method

1. Count cells per sample x cell type from integrated annotated H5AD
2. Merge with phenotype metadata (social_isolation_avg)
3. Run sccomp Bayesian GLM: `composition ~ social_isolation_avg + sex`
4. Report cell types with significant proportion shifts

## Running

```bash
source config/paths.sh

# Tsai
cd Analysis/SocIsl/CellTypeProportion/Tsai
sbatch run_prep.sh derived_batch
sbatch run_sccomp.sh derived_batch all fine

# DeJager
cd Analysis/SocIsl/CellTypeProportion/DeJager
sbatch run_prep.sh library_id
sbatch run_sccomp.sh library_id all fine
```

## Outputs

Under `${SOCISL_OUTPUT_ROOT}/CellTypeProportion/`:
- `data/` — count tables and metadata CSVs
- `results_{integration}/` — sccomp model RDS files
