# Analysis

This directory contains scripts for downstream analysis of processed snRNA-seq data.

## Overview

The Analysis phase takes fully processed, annotated snRNA-seq data and performs biological analyses including:

1. **Differential Expression (DEG)**: Identify differentially expressed genes between conditions
2. **SCENIC**: Single-cell regulatory network inference
3. **Transcription Factor Analysis**: TF activity and regulatory analysis

## Directory Structure

```
Analysis/
├── DeJager/
│   ├── DEG_Analysis/           # Differential expression
│   ├── SCENIC/                 # Regulatory network inference
│   └── Transcription_Factors/  # TF analysis
└── Tsai/
    ├── DEG_Analysis/           # Differential expression
    ├── SCENIC/                 # Regulatory network inference
    └── Transcription_Factors/  # TF analysis
```

## Analysis Types

### Differential Expression (DEG_Analysis)

Identify genes differentially expressed between:
- Disease status (AD vs Control)
- Cell types
- Other phenotypic groups

**Common Methods:**
- Pseudobulk analysis (recommended for snRNA-seq)
- Wilcoxon rank-sum test
- MAST (for zero-inflated data)
- edgeR/DESeq2 (pseudobulk)

**Typical Workflow:**
1. Aggregate counts to pseudobulk
2. Filter lowly expressed genes
3. Normalize counts
4. Run differential expression
5. Multiple testing correction
6. Pathway enrichment

### SCENIC

Single-Cell rEgulatory Network Inference and Clustering.

**Purpose:**
- Infer gene regulatory networks
- Identify active transcription factors
- Discover regulons (TF + target genes)

**Components:**
1. **GRNBoost2**: Co-expression network inference
2. **RcisTarget**: TF-motif enrichment
3. **AUCell**: Regulon activity scoring

**Prerequisites:**
- Processed count matrix
- Species-specific TF-motif databases

### Transcription Factor Analysis

Analyze transcription factor activity and regulation.

**Methods:**
- DoRothEA TF activity inference
- ChEA3 TF enrichment
- Motif analysis

## Data Requirements

### Input

Processed AnnData objects from the Processing phase:
- Cell type annotations (`obs['cell_type']`)
- Patient/sample IDs (`obs['patient_id']`)
- Phenotype information (`obs['condition']`)
- Normalized counts (`X` or `layers['normalized']`)
- Raw counts (`layers['counts']`)

### Gene Markers

Reference marker gene sets are available in:
```
/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Analysis/Tsai/Gene_Markers/
├── Brain_Human_PFC_Markers_Mohammadi2020.rds
├── human_mouse_ortho_HMD_HumanPhenotype.rds
└── mouse_pfc_markers_liftover_from_human.rds
```

## Resource Requirements

| Analysis | Cores | Memory | Time |
|----------|-------|--------|------|
| DEG (pseudobulk) | 8 | 64GB | 1-2h |
| SCENIC | 32+ | 256GB+ | 24-48h |
| TF Analysis | 8 | 64GB | 2-4h |

## Best Practices

### Differential Expression

1. Use pseudobulk for proper statistical modeling
2. Account for batch effects in the model
3. Use appropriate multiple testing correction (FDR < 0.05)
4. Validate top hits with literature

### SCENIC

1. Run on high-quality cells only
2. Use appropriate species-specific databases
3. Consider running separately per cell type
4. Validate regulons with known biology

## Status

> **Note**: This directory contains placeholder structure. Analysis scripts will be added as they are developed.

## Next Steps

1. Add DEG analysis scripts
2. Implement SCENIC pipeline
3. Add TF analysis workflows
4. Create visualization scripts

