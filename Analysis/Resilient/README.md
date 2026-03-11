# Resilient — Cognitive Resilience

Analysis of cognitive resilience in Alzheimer's disease — patients who maintain cognitive function despite significant AD neuropathology.

## Phenotype Definition

Resilience is defined using the combination of:
- **Cognitive diagnosis** (`cogdx`): clinical cognitive status at last visit
- **Neuropathology** (`braaksc`, `ceradsc`): postmortem AD pathology burden

Resilient individuals have high pathology but preserved cognition.

## Patient Selection

Phenotype data: `Data/Phenotypes/ROSMAP_clinical.csv`

| Group | Criteria | Description |
|-------|----------|-------------|
| Resilient | High pathology + no/mild cognitive impairment | Maintained function despite AD pathology |
| Non-resilient | High pathology + dementia | Expected cognitive decline given pathology |
| Control | Low pathology + no cognitive impairment | No significant AD pathology |

## Key Comparisons

- Resilient vs non-resilient (matched for pathology), stratified by cell type
- Resilient vs control (different pathology, similar cognition)
- Cell type-specific resilience signatures

## Directory Structure

```
Resilient/
├── DEG/
│   ├── DeJager/
│   └── Tsai/
├── TF/
│   ├── DeJager/
│   └── Tsai/
└── SCENIC/
    ├── DeJager/
    └── Tsai/
```
