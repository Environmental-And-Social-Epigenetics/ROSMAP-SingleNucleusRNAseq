# ACE GSEA — Gene Set Enrichment Analysis

WebGestaltR pathway enrichment analysis for ACE DEG results.

## Method

1. Take DEG results from `Analysis/ACE/DEG/` (DESeq2 outputs)
2. Rank genes by signed significance: `sign(log2FC) * -log10(pvalue)`
3. Run WebGestaltR GSEA against multiple pathway databases:
   - Gene Ontology (BP, CC, MF — redundancy-reduced)
   - KEGG, Reactome, Panther, Wikipathway
   - Transcription factor targets
4. Generate bubble plots showing enrichment landscape across cell types

## Prerequisites

- `${GSEA_ANALYSIS_ENV}` conda environment
- Completed DEG analysis (requires `.rda` result files)

## Running

```bash
source config/paths.sh
cd Analysis/ACE/GSEA/Tsai
sbatch run_enrichment.sh
```

## Outputs

Under `${ACE_OUTPUT_ROOT}/GSEA/Tsai/`:
- `overallWebGestaltRResult{Sex}{CellType}.rds` — enrichment results per cell type
- `{dataset}{Sex}Heatmap{pathway}.png` — bubble plot heatmaps

## Subdirectories

- `DeJager/` — GSEA on DeJager DEG results
- `Tsai/` — GSEA on Tsai DEG results
