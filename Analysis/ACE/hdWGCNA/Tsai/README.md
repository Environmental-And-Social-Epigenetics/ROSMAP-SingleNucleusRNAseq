# ACE hdWGCNA Co-Expression Networks - Tsai

Canonical entry point:

```bash
bash aceWgcnaT.sh
```

This launcher submits one SLURM job per cell type, processing priority cell
types with the most DEGs first.

## Method

1. Load cell-type-specific h5ad and convert to Seurat object
2. Construct metacells (k=25 cells per metacell) via `construct_metacells()`
3. Normalize, find variable features, scale data on metacells
4. WGCNA network construction:
   - Select soft threshold power via `pickSoftThreshold()`
   - Construct TOM and cluster into modules via `blockwiseModules()`
5. Calculate module eigengenes (ME) per metacell
6. Aggregate MEs to patient level (mean across metacells per patient)
7. Correlation: ME vs ACE phenotypes (Pearson r, FDR-corrected)
8. Module-DEG overlap: fraction of module genes that are DEGs
9. Module GO enrichment via clusterProfiler

## Priority Cell Types

| Sex | Cell Types (by DEG count) |
|-----|---------------------------|
| Male | Inh (353), Mic (327), Ast (196), In-PV_Basket (154) |
| Female | Oli (163), Ex-L2_3 (16), Ast (16), Exc (15) |

## Prerequisites

- Cell-type split h5ads: `${ACE_OUTPUT_ROOT}/DEG/Tsai/celltype_splits_${INTEGRATION}/`
- ACE phenotype CSV: `${ACE_SCORES_CSV}`
- `${WGCNA_ENV}` conda environment

## Setup

```bash
conda create -n wgcna_analysis r-base=4.3 r-seurat=5.0 r-devtools \
  r-ggplot2 r-wgcna r-biocmanager r-igraph
# In R:
# BiocManager::install(c("zellkonverter", "clusterProfiler", "org.Hs.eg.db"))
# devtools::install_github("smorabit/hdWGCNA", ref = "dev")
```

## Output

Results written to `${ACE_OUTPUT_ROOT}/hdWGCNA/Tsai/results_${INTEGRATION}/`.

```
results_derived_batch/
├── Male_Inh/
│   ├── module_assignments.csv       # Gene-to-module mapping
│   ├── module_eigengenes.csv        # ME values per patient
│   ├── module_trait_correlations.csv # ME vs ACE phenotype (r, p, padj)
│   ├── module_go_enrichment.csv     # GO-BP enrichment per module
│   ├── module_deg_overlap.csv       # DEG enrichment per module
│   ├── soft_threshold_plot.pdf
│   ├── dendrogram.pdf
│   └── wgcna_object.rds
├── Male_Mic/
└── ...
```

Canonical integration: `derived_batch`
