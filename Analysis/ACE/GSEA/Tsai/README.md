# ACE GSEA - Tsai Cohort

Gene Set Enrichment Analysis (GSEA) for Adverse Childhood Experience (ACE)
differential expression results from the Tsai cohort prefrontal cortex
snRNA-seq data.

## What This Does

GSEA identifies coordinated changes in biologically meaningful gene sets,
going beyond individual gene significance to detect pathway-level
dysregulation associated with ACE phenotypes.  Unlike over-representation
analysis (ORA), which requires an arbitrary significance cutoff, GSEA uses
the full ranked gene list and is therefore better powered for detecting
subtle, distributed effects across gene programs.

## Method

1. **Input**: DESeq2 pseudobulk DEG results (`.rda` files from
   `Analysis/ACE/DEG/Tsai/`)
2. **Gene ranking**: Each gene is scored as
   `sign(log2FoldChange) * -log10(pvalue)`.  This preserves effect direction
   and weights by statistical significance.
3. **Enrichment**: WebGestaltR runs ranked GSEA (Kolmogorov-Smirnov-like
   statistic with 1,000 permutations) against 8 curated databases.
4. **Output**: Per-cell-type enrichment results, a combined summary CSV, and
   publication-quality figures.

## Databases

| Database | Description |
|----------|-------------|
| `geneontology_Biological_Process_noRedundant` | GO-BP terms with redundancy removed via affinity propagation |
| `geneontology_Cellular_Component_noRedundant` | GO-CC terms (non-redundant) |
| `geneontology_Molecular_Function_noRedundant` | GO-MF terms (non-redundant) |
| `pathway_KEGG` | KEGG metabolic and signaling pathways |
| `pathway_Panther` | PANTHER classification pathways |
| `pathway_Reactome` | Reactome curated human pathways |
| `pathway_Wikipathway` | WikiPathways community-curated pathways |
| `network_Transcription_Factor_target` | Transcription factor target gene sets |

## Prerequisites

- Completed DEG analysis: `.rda` files must exist under
  `${ACE_OUTPUT_ROOT}/DEG/Tsai/results_${INTEGRATION}/${PHENOTYPE}/`
- Conda environment: `${GSEA_ANALYSIS_ENV}` with R, WebGestaltR, dplyr,
  ggplot2, tidyr, RColorBrewer

## Running

### Full pipeline (all phenotypes and sexes)

```bash
bash aceGseaT.sh                    # default: derived_batch integration
bash aceGseaT.sh derived_batch      # explicit integration
```

This submits 6 SLURM jobs (3 phenotypes x 2 sexes).

### Single phenotype/sex (manual)

```bash
sbatch run_enrichment.sh derived_batch tot_adverse_exp Female \
  "${ACE_OUTPUT_ROOT}/DEG/Tsai/results_derived_batch" \
  "${ACE_OUTPUT_ROOT}/GSEA/Tsai"
```

### Smoke test

```bash
bash smoke_test.sh
```

Creates a mock DESeq2 result and runs 2 databases in smoke mode.

## Output Structure

```
${ACE_OUTPUT_ROOT}/GSEA/Tsai/
  logs/                                  # SLURM job logs
  results_derived_batch/
    tot_adverse_exp/
      Female/
        gsea_summary.csv                 # Combined results (all cell types + databases)
        Exc_ranked_genes.csv             # Ranked gene list for Exc
        Exc_Biological_Process_noRedundant.rds   # Per-database RDS
        Exc_KEGG.rds
        ...
      Male/
        ...
    early_hh_ses/
      ...
    ace_aggregate/
      ...
```

### `gsea_summary.csv` Columns

| Column | Description |
|--------|-------------|
| `cell_type` | Cell type (e.g. Exc, Ast, Inh, Mic, Oli, Opc) |
| `sex` | Female or Male |
| `database` | Short database name |
| `geneSet` | Gene set identifier |
| `description` | Human-readable gene set name |
| `NES` | Normalized Enrichment Score |
| `pValue` | Nominal p-value |
| `FDR` | Benjamini-Hochberg adjusted p-value |
| `size` | Number of genes in the gene set |
| `leadingEdgeNum` | Number of leading edge genes driving enrichment |

## Interpretation

- **NES > 0**: Gene set is enriched among genes *upregulated* with higher
  ACE phenotype values (positive log2FC).
- **NES < 0**: Gene set is enriched among genes *downregulated* with higher
  ACE phenotype values (negative log2FC).
- **FDR < 0.05**: Strong evidence of enrichment.
- **FDR < 0.2**: Suggestive enrichment (the default WebGestaltR reporting
  threshold).
- **Leading edge genes**: The subset of genes that drive the enrichment
  signal.  These are the core genes contributing most to the NES.

## Scripts

| File | Purpose |
|------|---------|
| `aceGseaT.sh` | Main launcher -- submits SLURM jobs for all phenotypes/sexes |
| `run_enrichment.sh` | SLURM job wrapper -- activates env and calls R script |
| `gsea_analysis.R` | Core analysis -- gene ranking, WebGestaltR GSEA, summary CSV |
| `gsea_visualize.R` | Visualization -- bubble heatmap, barplots, convergence table |
| `smoke_test.sh` | End-to-end smoke test with mock data |
