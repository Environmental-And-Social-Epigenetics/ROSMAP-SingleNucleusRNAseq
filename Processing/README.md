# Processing

The processing layer converts CellBender-corrected matrices into integrated,
annotated AnnData objects for downstream analysis.

## Shared Three-Stage Structure

1. **Stage 1**: percentile-based QC filtering
2. **Stage 2**: `scDblFinder` doublet removal
3. **Stage 3**: normalization, HVG selection, PCA, Harmony, clustering, UMAP, ORA annotation

## Cohort-Specific Defaults

| Parameter | Tsai | DeJager |
|-----------|------|---------|
| Canonical Harmony key | `derived_batch` | `library_id` |
| Patient identity source | `patient_metadata.csv` | barcode map + overrides |
| Stage 3 output prefix | `tsai_` | `dejager_` |

Alternative DeJager Stage 3 batch keys remain available as sensitivity analyses:

- `patient_id`
- `pool_batch`
- `derived_batch`

## Shared Stage Parameters

### Stage 1

| Metric | Threshold |
|--------|-----------|
| `log1p_total_counts` | below 4.5th or above 96th percentile |
| `log1p_n_genes_by_counts` | below 5th percentile |
| `pct_counts_mt` | above 10% |

### Stage 2

| Parameter | Value |
|-----------|-------|
| Method | `scDblFinder` |
| Seed | `123` |
| Filter | retain singlets only |

### Stage 3

| Step | Parameters |
|------|------------|
| HVG selection | `flavor="seurat_v3"`, `n_top_genes=3000`, `layer="counts"` |
| PCA | `n_comps=30`, `svd_solver="arpack"` |
| Neighbors | `n_neighbors=30`, `n_pcs=30`, `metric="cosine"` |
| Leiden | `0.2`, `0.5`, `1.0` |
| UMAP | `min_dist=0.15`, `random_state=0` |
| ORA annotation | Mohammadi 2020 PFC markers, `use_raw=True` |

## Quick Start

```bash
cd Processing/Tsai/Pipeline && ./submit_pipeline.sh all
cd Processing/DeJager/Pipeline && ./submit_pipeline.sh all
```

See:

- [Tsai/Pipeline/README.md](Tsai/Pipeline/README.md)
- [DeJager/Pipeline/README.md](DeJager/Pipeline/README.md)
