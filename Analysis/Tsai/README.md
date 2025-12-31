# Tsai Analysis

Downstream analysis for the Tsai ROSMAP snRNA-seq dataset.

## Overview

This directory will contain analysis scripts for the Tsai dataset (ACE, Resilient, SocIsl cohorts), including differential expression, SCENIC regulatory network inference, and transcription factor analysis.

## Directory Structure

```
Tsai/
├── DEG_Analysis/           # Differential expression analysis
├── SCENIC/                 # Single-cell regulatory network inference
└── Transcription_Factors/  # Transcription factor activity analysis
```

## Cohorts

| Cohort | Focus Area |
|--------|------------|
| ACE | Adverse Childhood Experiences |
| Resilient | Cognitive Resilience |
| SocIsl | Social Isolation |

## Planned Analyses

### DEG_Analysis

**Planned comparisons:**
- AD vs Control (by cell type, by cohort)
- Resilient vs Non-resilient
- High ACE vs Low ACE
- Socially isolated vs Connected

**Methods:**
- Pseudobulk DESeq2/edgeR
- Cell type-specific analysis
- Cross-cohort meta-analysis
- Pathway enrichment (GSEA, GO, KEGG)

### SCENIC

**Goals:**
- Identify active regulons per cell type
- Compare TF activity across conditions
- Discover cohort-specific regulatory networks
- Cross-cohort regulon conservation

### Transcription_Factors

**Goals:**
- TF activity inference (DoRothEA)
- TF-target network analysis
- Integration with epigenetic data
- Cross-cohort TF comparison

## Data Dependencies

Requires processed data from `Processing/Tsai/`:
- Annotated h5ad files with cell type labels
- Patient metadata
- Clinical phenotypes

## Gene Markers

Reference marker gene sets are available locally:
```
/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Gene_Markers/
├── Brain_Human_PFC_Markers_Mohammadi2020.rds
├── human_mouse_ortho_HMD_HumanPhenotype.rds
└── mouse_pfc_markers_liftover_from_human.rds
```

## Status

> **Note**: This directory contains placeholder structure. Scripts will be added as analyses are developed.

## Related Work

Existing analysis notebooks in the original location:
```
/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Processing/
```

