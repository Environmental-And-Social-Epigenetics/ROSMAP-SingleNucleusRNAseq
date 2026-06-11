# ACE Mouse LNB Overlap - Tsai

Direction-aware overlap analysis between the Tsai human ACE pseudobulk DEG
results and mouse LNB-vs-control DEG results.

Canonical entry point:

```bash
bash aceMouseOverlapT.sh
```

Functional concordance entry point:

```bash
bash aceMouseFunctionalConcordanceT.sh
```

Focused male total-adverse functional concordance entry point:

```bash
bash aceMouseFunctionalTotAdverseMaleT.sh
```

By default this writes outputs to:

```text
${ANALYSIS_OUTPUT_ROOT}/ACE/MouseOverlap/Tsai
```

## Inputs

- Mouse DEG CSVs from `Analysis/ACE/Mouse_Results`.
- Human Tsai ACE DEG `.rda` files from the local `Analysis/ACE/DEG/Tsai`
  result folders, with `${ACE_OUTPUT_ROOT}/DEG/Tsai` as a fallback.
- Mouse-human orthologs from the MGI `HOM_MouseHumanSequence.rpt` report. The
  downloaded raw report, prepared one-to-one table, and metadata manifest are
  cached under the output directory.

## Sign Conventions

The mouse files use the supervisor-provided convention that negative
`log2FoldChange` is up in LNB. This workflow therefore reports:

```text
mouse_lnb_log2fc = -log2FoldChange
```

Mouse `LNB_up` is compared with human positive ACE `log2FoldChange` as the
main concordant direction. Human positive effects mean higher expression with
higher ACE phenotype values in the DESeq2 pseudobulk model.

## Comparisons

The workflow computes the full Cartesian set among:

- Mouse datasets: `NeuN_Dark`, `NeuN_Light`, `PV_Dark`, `PV_Light`
- Human targets: `Exc`, `Inh`, `In-PV_Basket`, `In-PV_Chandelier`
- Human phenotypes: `tot_adverse_exp`, `ace_aggregate`, `early_hh_ses`
- Human integrations: `derived_batch`, `projid`
- Human strata: `Fem`, `Male`

Comparison roles are labeled in the output:

- `primary`: `NeuN` to `Exc`, and `PV` to broad `Inh`
- `pv_sensitivity`: `PV` to `In-PV_Basket` or `In-PV_Chandelier`
- `exploratory_off_target`: all other cell-type pairings

The primary result to inspect first is `derived_batch`, `tot_adverse_exp`,
`primary`.

## Outputs

```text
human_deg_flat.csv
mouse_deg_normalized.csv
ortholog_mapping_used.csv
ortholog_mapping_metadata.json
overlap_summary.csv
overlap_genes.csv
inhibitory_proportion_context.csv
pathology_context.md
figures/
```

`overlap_summary.csv` includes directional Fisher and hypergeometric tests,
BH-adjusted p-values, Jaccard index, Spearman effect concordance, and signed
product rank-test summaries. `overlap_genes.csv` contains the gene-level rows
for every directional overlap.

## Functional Concordance

The exact ortholog DEG overlap is intentionally stringent. The functional
concordance workflow asks a broader question: whether mouse LNB and human ACE
shift the same GO Biological Process or Reactome pathways, even when different
individual genes pass the DEG threshold.

This workflow reuses:

```text
mouse_deg_normalized.csv
human_deg_flat.csv
ortholog_mapping_used.csv
```

Mouse ranks are computed as `mouse_lnb_stat = -mouse_stat`, so positive values
mean higher expression in LNB. Human ranks use the DESeq2 `stat`, so positive
values mean higher expression with the ACE phenotype. If a statistic is
missing, the fallback rank is `sign(log2FC) * -log10(pvalue)`.

Primary interpretation remains male human `derived_batch/tot_adverse_exp`:

- `NeuN_Dark` and `NeuN_Light` vs male `Exc`
- `PV_Dark` and `PV_Light` vs male broad `Inh`

Sensitivity rows include male PV subtype comparisons, female human strata, the
`projid` integration, and the `ace_aggregate` and `early_hh_ses` phenotypes.

Functional outputs are written under:

```text
${ACE_OUTPUT_ROOT}/MouseOverlap/Tsai/functional
```

Key files:

```text
functional_ranked_genes.csv
functional_gsea_results.csv
functional_concordance_summary.csv
functional_concordance_pathways.csv
functional_leading_edge_overlap.csv
figures/
```

`functional_concordance_summary.csv` reports NES correlations, pathway sign
concordance, concordant/discordant significant pathway counts, and directional
Fisher tests. `functional_concordance_pathways.csv` contains the joined
pathway-level mouse and human NES/FDR values.

## Male Total-Adverse Functional Concordance

`aceMouseFunctionalTotAdverseMaleT.sh` is the focused v2 workflow for comparing
the male LNB mice against the completed male human
`derived_batch/tot_adverse_exp` WebGestaltR GSEA outputs. It reads the human
RDS files directly from:

```text
${ACE_OUTPUT_ROOT}/GSEA/Tsai/results_derived_batch/tot_adverse_exp/Male
```

and therefore does not rely on `gsea_summary.csv`, which may be stale or empty.
Mouse genes are one-to-one ortholog mapped to human symbols and ranked as:

```text
mouse_lnb_stat = -mouse_stat
```

so positive mouse NES means enrichment among genes higher in LNB. The workflow
runs the same eight WebGestaltR libraries used in the human GSEA:

- GO Biological Process, Cellular Component, and Molecular Function
- KEGG, Panther, Reactome, WikiPathway
- Transcription Factor target

Mouse dark and light conditions remain separate. Primary comparisons are
`NeuN_Dark` and `NeuN_Light` against an aggregate of male human excitatory
subtype GSEA results, and `PV_Dark` and `PV_Light` against male human broad
`Inh`. PV basket/chandelier and individual excitatory subtype comparisons are
included as sensitivity rows.

Outputs are written under:

```text
${ACE_OUTPUT_ROOT}/MouseOverlap/Tsai/functional_tot_adverse_male
```

Key files:

```text
mouse_ranked_genes.csv
mouse_gsea_results.csv
human_gsea_reference.csv
human_ex_aggregate_reference.csv
functional_concordance_summary.csv
functional_concordance_pathways.csv
functional_concordant_terms.csv
functional_discordant_terms.csv
functional_leading_edge_overlap.csv
figures/
```

## Pathology Context

The human DEG model already adjusts for `niareagansc`, so DEG overlap is not a
naive pathology-unadjusted comparison. Inhibitory neuron proportion results are
summarized from the existing Tsai sccomp outputs when present. A non-AD-only
cell-proportion analysis is explicitly treated as follow-up work, not part of
this v1 overlap module.

## Smoke Test

```bash
bash smoke_test.sh
```

The smoke test builds tiny synthetic mouse, human, and MGI-like ortholog tables
and verifies mouse sign flipping, one-to-one ortholog filtering, universe
construction, and Fisher overlap counts.

Functional smoke test:

```bash
bash functional_smoke_test.sh
```

This test uses a tiny custom GMT to verify the functional sign convention,
one-to-one ortholog filtering, use of all finite ranked genes, and concordant
versus discordant pathway classification.

Focused male total-adverse smoke test:

```bash
bash functional_tot_adverse_male_smoke_test.sh
```

This test additionally verifies separate dark/light mouse datasets and
excitatory subtype aggregation.
