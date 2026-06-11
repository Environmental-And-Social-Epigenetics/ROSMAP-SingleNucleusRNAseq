# ACE Cell-Cell Communication - Tsai

Canonical entry point:

```bash
bash aceCellChatT.sh
```

This launcher:

1. Splits individuals into ACE-high and ACE-low groups (median split)
2. Runs CellChat on each group to infer ligand-receptor interactions
3. Compares communication probability between groups
4. Generates differential interaction results and visualizations

## Method

CellChat v2 is used to infer cell-cell communication from single-nucleus
RNA-seq data. The analysis:

1. Loads the full annotated h5ad (all cell types) for a given sex
2. Splits individuals by median of `tot_adverse_exp` into high/low groups
3. For each group: creates a CellChat object, identifies overexpressed
   ligands/receptors, computes communication probability (triMean)
4. Compares interaction strengths between groups
5. Focuses on biologically motivated axes: microglia-neuron,
   astrocyte-neuron, oligodendrocyte-neuron signaling

## Prerequisites

- Annotated integrated h5ad: `${TSAI_INTEGRATED}/tsai_annotated.h5ad`
- ACE phenotype CSV: `${ACE_SCORES_CSV}`
- `${CELLCHAT_ENV}` conda environment with CellChat, Seurat v5, zellkonverter

## Setup

The CellChat environment requires manual installation:

```bash
conda create -n cellchat_analysis r-base=4.3 r-seurat=5.0 r-devtools \
  r-ggplot2 r-patchwork r-ggalluvial r-nmf r-biocmanager
# In R:
# BiocManager::install("zellkonverter")
# devtools::install_github("jinworks/CellChat")
```

## Output

Results written to `${ACE_OUTPUT_ROOT}/CellChat/Tsai/results_${INTEGRATION}/`.

```
results_derived_batch/
├── Male/
│   ├── cellchat_high.rds          # CellChat object for ACE-high group
│   ├── cellchat_low.rds           # CellChat object for ACE-low group
│   ├── differential_interactions.csv
│   ├── pathway_changes.csv
│   └── focus_axes_results.csv     # Results for targeted interaction axes
└── Female/
    └── ...
```

## Interpretation

- Positive interaction change = stronger signaling in ACE-high group
- Negative interaction change = weaker signaling in ACE-high group
- Focus on complement (C3-C3AR1) for synaptic pruning hypothesis
- Focus on CX3CL1-CX3CR1 for microglia-neuron homeostatic signaling

Canonical integration: `derived_batch`
