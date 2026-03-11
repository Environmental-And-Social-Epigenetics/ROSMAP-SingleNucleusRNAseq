# Analysis

Downstream biological analyses of processed snRNA-seq data, organized by phenotype.

## Directory Structure

Each phenotype gets its own directory containing analysis type subdirectories, each split by dataset:

```
Analysis/
├── _template/          # Copy this to start a new phenotype analysis
├── ACE/                # Adverse Childhood Experiences
│   ├── DEG/
│   │   ├── DeJager/
│   │   └── Tsai/
│   ├── TF/
│   │   ├── DeJager/
│   │   └── Tsai/
│   └── SCENIC/
│       ├── DeJager/
│       └── Tsai/
├── Resilient/          # Cognitive Resilience
│   └── (same structure)
└── SocIsl/             # Social Isolation
    └── (same structure)
```

## Adding a New Phenotype Analysis

1. Copy the template: `cp -r _template/ NewPhenotype/`
2. Edit `NewPhenotype/README.md` with the phenotype definition and patient selection criteria
3. Add analysis scripts to the appropriate `DEG/`, `TF/`, or `SCENIC/` subdirectories, split by dataset

## Analysis Types

| Type | Method | Description |
|------|--------|-------------|
| **DEG** | Pseudobulk (DESeq2/edgeR) | Differential expression comparing conditions within cell types |
| **TF** | DoRothEA | Transcription factor activity inference and TF-target networks |
| **SCENIC** | pySCENIC | Single-cell regulatory network inference — active TFs and regulons per cell type |

## Data Requirements

- Annotated AnnData objects from the Processing phase (`obs['cell_type']`, `obs['projid']`)
- Clinical phenotype data from `Data/Phenotypes/` (see `${PHENOTYPE_DIR}` in `config/paths.sh`)
- Gene markers: `Processing/Tsai/Pipeline/Resources/Brain_Human_PFC_Markers_Mohammadi2020.rds`

## Current Phenotype Analyses

| Phenotype | Description | Status |
|-----------|-------------|--------|
| ACE | Adverse Childhood Experiences | Placeholder |
| Resilient | Cognitive Resilience despite AD pathology | Placeholder |
| SocIsl | Social Isolation | Placeholder |

## Resource Requirements

| Analysis | Cores | Memory | Time |
|----------|-------|--------|------|
| DEG (pseudobulk) | 8 | 64GB | 1-2h |
| SCENIC | 32+ | 256GB+ | 24-48h |
| TF Analysis | 8 | 64GB | 2-4h |
