# ACE DEG - DeJager

Canonical entry point:

```bash
bash aceDegDJ.sh
```

This launcher:

1. builds per-cell-type split `.h5ad` inputs from the canonical DeJager integrated object
2. runs pseudobulk DEG for the official ACE phenotype trio
3. writes outputs to `${ANALYSIS_OUTPUT_ROOT}/ACE/DEG/DeJager`

Useful helpers:

- `prep_celltype_splits.py`
- `run_deg.sh`
- `smoke_test.sh`

Canonical integration: `library_id`

Optional sensitivity integrations: `patient_id`, `pool_batch`, `derived_batch`
