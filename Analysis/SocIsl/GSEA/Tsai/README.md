# SocIsl GSEA — Tsai

WebGestaltR gene set enrichment analysis for Tsai SocIsl DEG results.

## Method

1. Load DESeq2 results (`.rda` files from `DEG/Tsai/`)
2. Rank genes: `sign(log2FC) * -log10(pvalue)`
3. Run WebGestaltR GSEA against: GO (BP, CC, MF), KEGG, Panther, Reactome,
   Wikipathway, TF targets
4. Aggregate results into bubble plots across cell types

## Running

```bash
source config/paths.sh
cd Analysis/SocIsl/GSEA/Tsai
sbatch gsea.sh
```

## Scripts

- `gseaResults.Rscript` — Main enrichment + visualization (DeJager data)
- `tsaiGseaResults.Rscript` — Tsai-specific enrichment
- `tsaiGseaResults_v2.Rscript` — Updated Tsai visualization

## Outputs

- `overallWebGestaltRResult{Sex}{CellType}.rds` — per cell type enrichment
- `plotGsea{Sex}{Pathway}.csv` — aggregated results for plotting
