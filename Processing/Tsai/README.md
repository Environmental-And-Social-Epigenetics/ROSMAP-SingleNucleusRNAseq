# Tsai Processing

Processing scripts for the Tsai ROSMAP snRNA-seq dataset.

## Primary Pipeline

The current production pipeline lives in **`Pipeline/`** and is documented in
[Pipeline/README.md](Pipeline/README.md).  It runs in three stages:

1. **QC Filtering** — percentile-based outlier removal + mitochondrial filter
2. **Doublet Removal** — scDblFinder
3. **Integration & Annotation** — Harmony batch correction, Leiden clustering,
   ORA-based cell type annotation (Mohammadi 2020 markers)

```bash
cd Pipeline
sbatch 01_qc_filter.sh
sbatch 02_doublet_removal.sh
sbatch 03_integration_annotation.sh
```

## Directory Structure

```
Tsai/
├── Pipeline/                   # Primary 3-stage pipeline (current)
│   ├── 01_qc_filter.*
│   ├── 02_doublet_removal.*
│   ├── 03_integration_annotation.*
│   ├── Resources/              # Marker gene references
│   └── envs/                   # Conda environment specs
├── QC/                         # Legacy per-sample scripts (reference only)
│   ├── Doublets/
│   └── Outliers/
├── Batch_Correction/           # Legacy batch correction scripts (reference only)
├── Norm_FeatSel_DimRed_Cluster.ipynb   # Prototyping notebook
└── Tsai_Sample_SingleCell_Pipeline.ipynb
```

## Resource Requirements

| Stage | Cores | Memory | Time |
|-------|-------|--------|------|
| 1 — QC Filtering | 4 | 32GB | ~1h per sample (array job) |
| 2 — Doublet Removal | 4 | 32GB | ~1h per sample (array job) |
| 3 — Integration | 32 | 500GB | 4-8h (single job) |

## Legacy Scripts

The `QC/` and `Batch_Correction/` directories contain earlier development
scripts that preceded the unified Pipeline.  They are kept for reference but
are **not** part of the current workflow.

