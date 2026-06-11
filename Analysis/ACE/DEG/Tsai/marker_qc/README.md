# Marker QC for Tsai snRNA cell-type annotation

This directory holds a canonical-marker visualization of the Tsai integrated
single-nucleus object, used to sanity-check the cell-type annotation that feeds
the ACE DEG pipeline.

## Why this exists

Some ACE DEG hits look cell-type-implausible (e.g. NPY/SST among microglia,
glial-flavored genes among PV-basket). This QC asks two questions per cell
type:

1. **Do its own canonical markers express where expected?** (`MISSING_OWN` flag
   in [tables/marker_flags.csv](tables/marker_flags.csv))
2. **Do markers from *other* lineages express at suspiciously high rates?**
   (`UNEXPECTED_OTHER` flag in the same file)

## How cell-type labels were assigned

`adata.obs["cell_type"]` is set by **over-representation analysis (ORA) using
the Mohammadi 2020 PFC marker set**, not by label transfer / scANVI / celltypist:

- Markers loaded from
  [Resources/Brain_Human_PFC_Markers_Mohammadi2020.rds](../../../../../Processing/Tsai/Pipeline/Resources/Brain_Human_PFC_Markers_Mohammadi2020.rds)
  via rpy2 in
  [03_integration_annotation.py:157-174](../../../../../Processing/Tsai/Pipeline/03_integration_annotation.py#L157-L174).
- `decoupler.run_ora(adata, ...)` against `adata.raw`
  ([line 229](../../../../../Processing/Tsai/Pipeline/03_integration_annotation.py#L229)).
- Sources ranked per `leiden_res0_5` cluster via
  `dc.rank_sources_groups(..., method="t-test_overestim_var")`
  ([lines 245-250](../../../../../Processing/Tsai/Pipeline/03_integration_annotation.py#L245-L250)).
- Top-1 source per cluster becomes the cell-type label for every cell in that
  cluster
  ([lines 256-260](../../../../../Processing/Tsai/Pipeline/03_integration_annotation.py#L256-L260)).

**Important consequence:** a label like `In-PV_Basket` means "this cluster's
ORA score was highest for the Mohammadi PV-basket signature," not "these cells
are marker-validated PV-basket neurons." Cross-validating with the canonical
panel in this QC is therefore essential.

Per-cluster ranking tables (`cluster_annotation_rankings.csv`,
`cluster_annotation_top3.csv`) are written to the integration output dir at
`${TSAI_INTEGRATED}/`; check those for clusters where the top vs second source
are close.

## What the script produces

Run via `sbatch make_marker_dotplots.sh` (or `bash make_marker_dotplots.sh`
interactively):

### Figures (PDFs in `figures/`)

- `dotplot_all_categories.{pdf,png}` — primary dotplot. Rows = cell types,
  columns = canonical markers grouped/bracketed by the 5 lineage categories.
- `dotplot_all_categories_swapped.pdf` — same data, axes swapped (portrait).
- `dotplot_<category>.pdf` — per-category dotplot at higher resolution
  (categories: `PV_inhibitory`, `Microglia`, `Oligo_myelin`, `Astrocyte`,
  `Excitatory`).
- `stacked_violin_<category>.pdf` and `stacked_violin_all_categories.pdf` —
  matched stacked-violin views.

All plots use `use_raw=True` and `standard_scale="var"`. The dotplot dot size
encodes fraction of cells expressing the gene; color encodes the column-scaled
mean expression.

### Tables (CSVs in `tables/`)

- `markers_present.csv` — which canonical markers are present in
  `adata.raw.var_names` (markers absent from the assay are flagged here).
- `marker_fraction_expressing.csv` — per cell type × marker, fraction of cells
  with raw expression > 0. Uses the FULL adata (not the subsample) for accurate
  percentages.
- `marker_mean_expression.csv` — per cell type × marker, mean
  log-normalized expression.
- `marker_flags.csv` — flagged (cell_type, gene, category, fraction, flag_type)
  rows where `flag_type ∈ {MISSING_OWN, UNEXPECTED_OTHER}`. Default threshold
  30% — override with `THRESHOLD=0.25 sbatch make_marker_dotplots.sh`.
- `run_metadata.json` — input path, cell counts, seed, scanpy/anndata
  versions, timestamp.

## Inputs & assumptions

- Default input: `${TSAI_INTEGRATED}/tsai_annotated.h5ad` (resolved from
  [config/paths.sh:195](../../../../../config/paths.sh#L195)). Override with
  `--input-h5ad` or `INPUT_H5AD=...`.
- `adata.raw` must be populated (it is, from
  [03_integration_annotation.py:552](../../../../../Processing/Tsai/Pipeline/03_integration_annotation.py#L552)).
  Plots use `adata.raw` because several markers (KCNS3, MOBP, ERBB4, ATP1A2,
  …) are not in the 3000-HVG matrix at `adata.X`.
- Subsampling: `--per-celltype-cap 20000` (default) caps each cell type
  independently. This keeps plotting tractable without underweighting rare
  types (Endo, In-PV_Chandelier) compared to a global `sc.pp.subsample`.
  Fractions in CSVs use the FULL adata regardless of cap.

## Knobs

| Variable / flag      | Default | Meaning                                |
|----------------------|---------|----------------------------------------|
| `--per-celltype-cap` | 20000   | Max cells per cell type for plotting   |
| `--threshold`        | 0.30    | Fraction cutoff for marker_flags       |
| `--seed`             | 0       | Subsample RNG seed                     |
| `INPUT_H5AD` env     | `${TSAI_INTEGRATED}/tsai_annotated.h5ad` | Override input |

## Follow-up (lower priority)

- **Donor-level marker QC.** Group by (`projid`, `cell_type`) instead of
  `cell_type` alone — surfaces whether the unexpected-marker signal in e.g.
  Mic is concentrated in a small subset of donors. A sibling script
  `donor_marker_qc.py` (not yet written) would add this.
- **Residual-doublet check.** `${TSAI_DOUBLET_REMOVED}/doublet_summary.csv`
  already records per-sample scDblFinder scores. Adding a per-`cell_type`
  aggregation (median residual doublet score among kept singlets) would
  highlight cell types disproportionately enriched in barely-sub-threshold
  cells.
