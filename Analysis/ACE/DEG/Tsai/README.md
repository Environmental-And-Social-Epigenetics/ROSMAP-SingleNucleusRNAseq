# ACE DEG - Tsai

Canonical entry point:

```bash
bash aceDegT.sh
```

This launcher:

1. builds per-cell-type split `.h5ad` inputs from the canonical Tsai integrated object
2. runs pseudobulk DEG for the official ACE phenotype trio
3. writes outputs to `${ANALYSIS_OUTPUT_ROOT}/ACE/DEG/Tsai`

Useful helpers:

- `prep_celltype_splits.py`
- `run_deg.sh`
- `smoke_test.sh`

Canonical integration: `derived_batch`

Optional sensitivity integration: `projid`

Legacy exploratory helpers live under `legacy/` and are not part of the
supported workflow or smoke-test surface.
