# Resilient — Cognitive Resilience

Analysis of cognitive resilience in Alzheimer's disease — patients who maintain
cognitive function despite significant AD neuropathology.

## Phenotype Contract

**Source CSV:** `${RESILIENT_PHENOTYPE_CSV}` (default: `Data/Phenotypes/ROSMAP_clinical.csv`)

**Required columns:** `projid`, `cogdx`, `braaksc`, `ceradsc`, `msex`,
`age_death`, `pmi`, `niareagansc`

**Resilience group derivation:**

| Group | cogdx | braaksc | ceradsc | Description |
|-------|-------|---------|---------|-------------|
| Resilient | 1-2 (NCI/MCI) | >= 4 | <= 2 | High pathology + preserved cognition |
| Non-resilient | 4-5 (dementia) | >= 4 | <= 2 | High pathology + dementia |
| Control | 1 (NCI) | <= 2 | >= 3 | Low pathology + no impairment |

**DESeq2 design formula:** `~ age_death + pmi + resilience_group + niareagansc`

**sccomp formula:** `composition ~ resilience_group + sex`

## Key Comparisons

- Resilient vs non-resilient (matched for pathology), stratified by cell type
- Resilient vs control (different pathology, similar cognition)
- Cell type-specific resilience signatures

## Supported Workflows

| Analysis | Tsai | DeJager | Status |
|----------|------|---------|--------|
| DEG (pseudobulk DESeq2) | `DEG/Tsai/` | `DEG/DeJager/` | Implemented |
| Cell-type proportion (sccomp) | `CellTypeProportion/Tsai/` | `CellTypeProportion/DeJager/` | Implemented |

## Quick Start

```bash
source config/paths.sh

# Tsai DEG pipeline
cd Analysis/Resilient/DEG/Tsai
bash pipeline.sh --integration derived_batch --stage all

# Tsai cell-type proportion
cd Analysis/Resilient/CellTypeProportion/Tsai
sbatch run_prep.sh
sbatch run_sccomp.sh
```

## Outputs

All outputs go under `${RESILIENT_OUTPUT_ROOT}` (default:
`${ANALYSIS_OUTPUT_ROOT}/Resilient/`):

```
Resilient/
├── DEG/
│   ├── Tsai/
│   │   ├── celltype_splits_{integration}/
│   │   └── results_{integration}/resilience_group/
│   └── DeJager/
│       └── (same structure)
└── CellTypeProportion/
    ├── Tsai/
    │   ├── data/
    │   └── results_{integration}/
    └── DeJager/
```

## Directory Structure

```
Resilient/
├── README.md
├── DEG/
│   ├── Tsai/
│   │   ├── prep_celltype_splits.py
│   │   ├── resilDegT.Rscript
│   │   ├── resilDegT.sh
│   │   ├── run_deg.sh
│   │   ├── pipeline.conf
│   │   └── pipeline.sh
│   └── DeJager/
│       ├── prep_celltype_splits.py
│       ├── resilDegDJ.Rscript
│       └── resilDegDJ.sh
└── CellTypeProportion/
    ├── Tsai/
    │   ├── prep_counts.py
    │   ├── sccomp_analysis.R
    │   ├── run_prep.sh
    │   └── run_sccomp.sh
    └── DeJager/
        ├── prep_counts.py
        ├── sccomp_analysis.R
        ├── run_prep.sh
        └── run_sccomp.sh
```
