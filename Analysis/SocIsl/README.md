# SocIsl — Social Isolation

Analysis of the impact of social isolation on gene expression in Alzheimer's disease brain tissue.

## Phenotype Definition

Social isolation is measured using ROSMAP variables:
- **Social network size** (`soc_net_bl`, `soc_net_lv`)
- **Social isolation score** (`social_isolation_avg`, `social_isolation_lv`)

## Patient Selection

Phenotype data: `Data/Phenotypes/dataset_652_basic_12-23-2021.csv`

| Group | Criteria | Description |
|-------|----------|-------------|
| Socially isolated | High `social_isolation_avg` | Patients with elevated social isolation scores |
| Non-isolated | Low `social_isolation_avg` | Patients with typical social engagement |

## Key Comparisons

- Isolated vs non-isolated, stratified by cell type
- Continuous association of social isolation score with gene expression
- Interaction with AD pathology

## Directory Structure

```
SocIsl/
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
