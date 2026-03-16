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

## Quick Start: From Processing Output to Analysis

### What the Processing Pipeline Produces

The Processing pipeline (`Processing/Tsai/Pipeline/` or `Processing/DeJager/Pipeline/`) outputs
an annotated AnnData object with these key fields in `obs`:

| Field | Description | Source |
|-------|-------------|--------|
| `cell_type` | ORA-annotated cell type (Ast, Exc, Inh, Mic, Oli, OPC) | Stage 3 annotation |
| `projid` | Patient identifier | Metadata CSV |
| `leiden_res0_2`, `leiden_res0_5`, `leiden_res1` | Leiden clusters at multiple resolutions | Stage 3 clustering |
| Clinical columns (e.g., `msex`, `age_death`, `pmi`) | From `patient_metadata.csv` | Attached in Stage 3 |

The final files are:
- **Tsai**: `${TSAI_PROCESSING_OUTPUTS}/03_Integrated/tsai_annotated.h5ad` (~83 GB)
- **DeJager**: `${DEJAGER_PROCESSING_OUTPUTS}/03_Integrated/dejager_annotated.h5ad` (not yet generated)

### Loading Data for Analysis

```python
import scanpy as sc
import pandas as pd

# Load the annotated object from Processing
adata = sc.read_h5ad("Tsai_Data/Processing_Outputs/03_Integrated/tsai_annotated.h5ad")

# Attach phenotype data (e.g., social isolation scores)
pheno = pd.read_csv("Data/Phenotypes/dataset_652_basic_12-23-2021.csv")
adata.obs = adata.obs.merge(pheno[['projid', 'social_isolation_avg']], on='projid', how='left')
```

### Running Your First Analysis (DEG Example)

```bash
source config/paths.sh
init_conda
conda activate "${CONDA_ENV_BASE}/deg_analysis"

# Run Social Isolation DEG on Tsai dataset
cd Analysis/SocIsl/DEG/Tsai
sbatch socIslDegT.sh
```

### Conda Environments for Analysis

Install with `bash setup/install_envs.sh --analysis`. See [Analysis/envs/README.md](envs/README.md) for details.

| Environment | Variable | Used For |
|-------------|----------|----------|
| `deg_analysis` | `${CONDA_ENV_BASE}/deg_analysis` | DEG (limma, DESeq2, edgeR) |
| `scenic_analysis` | `${CONDA_ENV_BASE}/scenic_analysis` | pySCENIC regulatory networks |
| `compass_analysis` | `${CONDA_ENV_BASE}/compass_analysis` | COMPASS metabolic flux |
| `gsea_analysis` | `${CONDA_ENV_BASE}/gsea_analysis` | WebGestaltR gene set enrichment |

---

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

| Phenotype | Description | Status | Notes |
|-----------|-------------|--------|-------|
| ACE | Adverse Childhood Experiences | Placeholder | Directory structure only |
| Resilient | Cognitive Resilience despite AD pathology | Placeholder | Directory structure only |
| SocIsl | Social Isolation | Legacy scripts migrated | DEG, SCENIC, COMPASS, GSEA scripts + result CSVs from Openmind |

## Resource Requirements

| Analysis | Cores | Memory | Time |
|----------|-------|--------|------|
| DEG (pseudobulk) | 8 | 64GB | 1-2h |
| SCENIC | 32+ | 256GB+ | 24-48h |
| TF Analysis | 8 | 64GB | 2-4h |
