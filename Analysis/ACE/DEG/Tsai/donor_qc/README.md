# Donor-level + annotation QC for Tsai snRNA cell-type labels

Diagnoses **why** the Mohammadi-ORA cell-type labels in
`tsai_annotated.h5ad` look off for some cell types (Mic, OPC, deep-layer Ex),
and **whether the bad signal is concentrated in a small donor subset** or
cohort-wide.

## How to run

```bash
# Smoke (50k cells, Section A only — ~2 min):
SMOKE_FLAG=--smoke sbatch donor_qc.sh

# Single problematic cell type drill-down (~10 min):
CELLTYPE=Mic sbatch donor_qc.sh

# Full pass (~30–90 min, all 6 sections):
sbatch donor_qc.sh

# Build the report after the full pass:
python build_report.py
```

## What the six diagnostic sections answer

| § | Check | Output | Question it answers |
|---|-------|--------|---------------------|
| A | Per-donor QC table + cohort distribution | `tables/donor_qc.csv`, `tables/donor_outliers.csv`, `figures/cohort_qc_violins.{pdf,png}` | Are there donors with systematically bad QC metrics? |
| B | Cluster × donor crosstab + inverse-Simpson | `tables/cluster_donor_{counts,composition,skew}.csv`, `figures/cluster_donor_heatmap.{pdf,png}` | Is any cluster dominated by 1–2 donors? |
| C | Per-cluster top1 vs top2 ORA gap + per-cell ORA gap | `tables/cluster_annotation_gap.csv`, `tables/per_cell_ora_gap_quantiles.csv`, `figures/per_cell_ora_gap_by_celltype.{pdf,png}` | How confident was each cluster's annotation call? Which cells sit between identities? |
| D | Mohammadi-panel (pipeline's own panel) fraction-expressing + canonical vs Mohammadi comparison | `tables/marker_fraction_expressing_mohammadi_all.csv`, `tables/mohammadi_own_source_fractions.csv`, `tables/panel_comparison.csv` | Are the canonical-panel flags real mis-annotation or just panel-mismatch artifacts? |
| E | Per-donor failing-marker fractions + drop-one-donor sensitivity | `tables/donor_problem_markers_<ct>.csv`, `tables/drop_one_donor_sensitivity.csv`, `figures/drop_one_donor_sensitivity_<ct>.{pdf,png}`, `figures/donor_qc_vs_problem_marker_<ct>.pdf` | **Is the problem concentrated in a few donors?** |
| F | Per-cell-type doublet score distributions | `tables/doublet_score_by_celltype.csv`, `figures/doublet_score_by_celltype.{pdf,png}` | Are problematic cell types doublet-contaminated? |

## What each check teaches (brief)

- **Section A (per-donor QC):** Donors with high MT% have stressed/dying
  nuclei → noisy expression and unstable cluster placement. Low n_genes →
  poor library or ambient-RNA-heavy nuclei. High doublet rate →
  under-dissociated samples. Donors flagged on ≥2 metrics are prime
  suspects.

- **Section B (cluster × donor):** A healthy cluster draws cells from many
  donors at roughly the per-donor base rate. A cluster dominated by 1–2
  donors is suspect — the "cell type" label is likely reflecting those
  donors' technical state, not a real shared biology.

- **Section C (ORA confidence):** ORA scores how well each marker signature
  matches a cell's expressed genes. A confident assignment has a large gap
  between the top-ranked signature and the next-best. Cells with small gaps
  sit between cell-type identities — most mis-annotation lives there.

- **Section D (panel comparison):** The first-pass marker QC used canonical
  markers the mentor specified, *not* the Mohammadi markers the pipeline
  used. So a "missing own marker" flag can mean two different things:
  - **Mohammadi passes, canonical fails:** the label is internally
    consistent with the pipeline; the canonical-panel flag is a
    panel-mismatch artifact (different reference panel = different
    expected markers).
  - **Both fail:** the cell type is likely genuinely mis-annotated.

- **Section E (drop-one-donor):** Iteratively removes the worst donor (by
  mean fraction expressing the failing markers, weighted by donor's cell
  count) and recomputes the cohort-wide fraction. A sharp jump in the first
  few drops means a small donor subset is driving the problem (fix: exclude
  them). A slow climb means the problem is cohort-wide (fix: redo the
  annotation pipeline).

- **Section F (doublets):** scDblFinder.score is a per-sample-relative
  call. Cell types concentrated in cells with high doublet scores are
  likely partly capturing unresolved doublets. OPC and Mic are classic
  offenders because both are small populations adjacent to large ones
  (Oli, neurons).

## Inputs (all already on disk)

- `${TSAI_INTEGRATED}/tsai_annotated.h5ad` — has per-cell QC metrics
  (`pct_counts_mt`, `n_genes_by_counts`, `total_counts`, `pct_counts_ribo`),
  doublet scores (`scDblFinder.score`, `scDblFinder.class`), the labels
  (`cell_type`, `leiden_res0_5`, `projid`), and the ORA estimate matrix
  (`obsm['ora_estimate']`).
- `${TSAI_QC_FILTERED}/qc_summary.csv` — per-sample post-QC counts.
- `${TSAI_DOUBLET_REMOVED}/doublet_summary.csv` — per-sample doublet rates.
- `${TSAI_INTEGRATED}/cluster_annotation_rankings.csv` — Mohammadi-signature
  rankings per Leiden cluster (for the top1-top2 gap in Section C).
- `Brain_Human_PFC_Markers_Mohammadi2020.rds` — the marker source the
  pipeline used (loaded via rpy2 for Section D).

## Outputs

After a full run + `build_report.py`:

- `donor_qc_report.md` — single shareable markdown for the mentor.
- `tables/` — every CSV referenced in the report.
- `figures/` — every PDF/PNG referenced in the report.
