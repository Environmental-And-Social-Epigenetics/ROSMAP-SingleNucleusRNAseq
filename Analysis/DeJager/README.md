# DeJager Analysis

Downstream analysis for the DeJager ROSMAP snRNA-seq dataset.

## Overview

This directory will contain analysis scripts for the DeJager dataset, including differential expression, SCENIC regulatory network inference, and transcription factor analysis.

## Directory Structure

```
DeJager/
├── DEG_Analysis/           # Differential expression analysis
├── SCENIC/                 # Single-cell regulatory network inference
└── Transcription_Factors/  # Transcription factor activity analysis
```

## Planned Analyses

### DEG_Analysis

**Planned comparisons:**
- AD vs Control (by cell type)
- Resilient vs Non-resilient
- APOE genotype effects

**Methods:**
- Pseudobulk DESeq2/edgeR
- Cell type-specific analysis
- Pathway enrichment (GSEA, GO)

### SCENIC

**Goals:**
- Identify active regulons per cell type
- Compare TF activity across conditions
- Discover disease-specific regulatory networks

### Transcription_Factors

**Goals:**
- TF activity inference (DoRothEA)
- TF-target network analysis
- Integration with GWAS hits

## Data Dependencies

Requires processed data from `Processing/DeJager/`:
- Annotated h5ad files with cell type labels
- Patient metadata
- Clinical phenotypes

## Status

> **Note**: This directory contains placeholder structure. Scripts will be added as analyses are developed.

## Related Data

Individual patient analyses from the Processing phase are available at:
```
Processing/DeJager/Individual_Analyses/
```

