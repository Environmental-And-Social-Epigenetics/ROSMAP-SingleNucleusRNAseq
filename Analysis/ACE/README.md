# ACE — Adverse Childhood Experiences

Analysis of the impact of adverse childhood experiences on gene expression in Alzheimer's disease brain tissue.

## Phenotype Definition

ACE measures early-life adversity using five components from the ROSMAP dataset:
- **Emotional neglect** (`emotional_neglect`)
- **Family problems/separation** (`family_pro_sep`)
- **Financial need** (`financial_need`)
- **Parental intimidation** (`parental_intimidation`)
- **Parental violence** (`parental_violence`)

Total ACE score: `tot_adverse_exp` (sum of components)

## Patient Selection

Phenotype data: `Data/Phenotypes/TSAI_DEJAGER_all_patients_wACEscores.csv` (296 patients with ACE scores)

| Group | Criteria | Description |
|-------|----------|-------------|
| ACE-exposed | `tot_adverse_exp > 0` | Patients reporting any adverse childhood experiences |
| Non-exposed | `tot_adverse_exp == 0` | No reported adverse experiences |

## Key Comparisons

- ACE-exposed vs non-exposed, stratified by cell type
- Dose-response: correlation of total ACE score with gene expression
- ACE component-specific effects (e.g., parental violence alone)

## Directory Structure

```
ACE/
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
