# Cell-type "confusion" QC for ACE DEG hits

Detects whether a **canonical marker gene for cell type A turns up as a strong DEG
hit in a different cell type B** — a signal of cell-type annotation
contamination/ambiguity in the integrated object's `cell_type` labels.

This operates on the **DEG result CSVs** (which genes ACE significantly moves per
cell type). It is **distinct from `../marker_qc/`**, which checks marker *expression
fractions* in the annotated h5ad. Here the question is: among the genes ACE
significantly moves in cell type B, are any strong hits canonical markers of a
*different* lineage?

## Run

```bash
bash celltype_confusion.sh                 # MaleContAD, padj<0.05 & |log2FC|>=1, both sources
ARM=MaleNoADadj bash celltype_confusion.sh # a different model arm
PADJ=0.05 LFC=1.0 SOURCE=panel bash celltype_confusion.sh
```

Lightweight (reads small CSVs only) — runs on a login node. Needs `$NEBULA_ENV`
(pandas + scipy + matplotlib).

## Two marker sources

- **`panel`** — the collaborator canonical panel (PV/inhibitory, microglia,
  oligo/myelin, astrocyte, excitatory). Small and lineage-specific → **primary
  evidence**, because it is independent of how the cells were annotated.
- **`mohammadi`** — the Mohammadi 2020 marker list, which was the actual basis for
  the ORA annotation. Cross-firing here is **partly circular** (a marker driving an
  annotation will naturally be "expressed" in that type) → treat as **corroborating
  only**. Restricted to *specific* markers (default: genes marking ≤2 types,
  `--spec-min 0.5`) to avoid the "every gene is somebody's marker" noise.

## Scoring

A gene `g` is a **confusion event** in DEG cell type `B` when it is a strong hit
(`padj < 0.05 AND |log2FC| ≥ 1`, configurable) and is canonical to a cell
type/category that does **not** include `B`. Events are split by:

- **cross_lineage = True** (headline): owner lineage ≠ B's lineage — e.g. an
  inhibitory marker (SST) strong in Microglia, or an excitatory marker in astrocytes.
- **cross_lineage = False** (noted, down-weighted): within-lineage neighbor — e.g.
  an Ex-subtype marker strong in another Ex-subtype, which is largely expected.

`severity = specificity × (2.0 if cross_lineage else 0.5)`.

Per cell type, a **hypergeometric enrichment** test asks whether cross-lineage
specific markers are over-represented among the strong hits (vs. their share of the
tested gene universe). `fold_enrichment > 1` with `padj < 0.05` flags a genuinely
"confused" cell type. Cell types with `< 5` strong hits get `NA` fold/p (noted).

## Outputs

`tables/`
- `confusion_events_<source>_<arm>.csv` — one row per (cell_type, gene) event:
  `deg_cell_type, deg_lineage, gene, owner, owner_lineage, marker_source, rank,
   padj, log2FC, baseMean, cross_lineage, specificity, severity, hit_reason`.
- `confusion_summary_<source>_<arm>.csv` — per cell type: hit counts, cross-lineage
  fraction, hypergeometric `M/K/n/observed_k/expected_k/fold_enrichment/pvalue/padj`,
  `top_cross_genes`, `note`.
- `confusion_matrix_<source>_<arm>.csv` — cell type (rows) × marker owner (cols):
  count of strong hits that are markers of each owner. Off-diagonal = confusion.

`figures/`
- `confusion_heatmap_<source>_<arm>.{png,pdf}` — heatmap of the matrix; bright
  off-diagonal cells = contamination.

`run_metadata.json` — parameters, paths, marker-universe sizes.

## Interpreting

- **Diagonal should dominate.** A clean cell type's strong hits are mostly its own
  category's markers, with few/zero cross-lineage flags.
- **Known positive control:** in `MaleContAD`, `SST` (Mohammadi In-SST–specific) is
  the #2 ACE DEG in `Mic`. It appears in `confusion_events_mohammadi_MaleContAD.csv`
  with `deg_cell_type=Mic, cross_lineage=True`. (It is **not** in the collaborator
  panel, so the panel source won't flag it — an expected source difference.)
- Cross-reference flagged cell types with the `ace_vs_ad_beta_concordance` outliers
  (Mic, OPC, Oli) from `../report/`.

## Caveats

- The Mohammadi panel folds inhibitory subtypes; the collaborator panel maps all
  `In-*` to a single `PV_inhibitory` category (low resolution for inhibitory
  subtypes). Mohammadi gives subtype resolution there.
- `Inh` has no single Mohammadi key → its "own" markers are the union of all `In-*`
  sources; only excitatory/glial markers count as cross-lineage there.
