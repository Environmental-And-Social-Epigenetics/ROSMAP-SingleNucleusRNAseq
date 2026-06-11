# ACE GSEA -- Gene Set Enrichment Analysis

WebGestaltR ranked GSEA for ACE differential expression results.

## Method

Gene-level DESeq2 results are ranked by `sign(log2FC) * -log10(pvalue)` and
tested against 8 pathway/network databases (GO-BP, GO-CC, GO-MF, KEGG,
Panther, Reactome, Wikipathway, TF targets) using the WebGestaltR GSEA
implementation.

## Cohorts

- **[Tsai/](Tsai/)** -- Tsai cohort prefrontal cortex snRNA-seq
- **DeJager/** -- DeJager cohort (to be implemented)

## Prerequisites

- Completed DEG analysis (`Analysis/ACE/DEG/`)
- `${GSEA_ANALYSIS_ENV}` conda environment with R and WebGestaltR

## Quick Start

```bash
# Tsai cohort (submits 6 SLURM jobs)
cd Analysis/ACE/GSEA/Tsai
bash aceGseaT.sh
```

## Outputs

Results are written to `${ACE_OUTPUT_ROOT}/GSEA/{Tsai,DeJager}/` with a
`gsea_summary.csv` per phenotype/sex combination containing NES, FDR, and
leading edge information for each cell type and database.
