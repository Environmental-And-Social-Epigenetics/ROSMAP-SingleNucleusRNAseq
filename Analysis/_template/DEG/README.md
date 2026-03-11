# DEG Analysis — {Phenotype}

Pseudobulk differential expression analysis comparing {phenotype groups} within each cell type.

## Method

- Aggregate counts per patient per cell type (pseudobulk)
- DESeq2 or edgeR for differential testing
- Covariates: age, sex, PMI, batch

## Subdirectories

- `DeJager/` — DEG analysis on DeJager dataset
- `Tsai/` — DEG analysis on Tsai dataset
