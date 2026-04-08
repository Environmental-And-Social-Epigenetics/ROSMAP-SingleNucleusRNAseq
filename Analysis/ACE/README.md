# ACE

ACE analysis studies adverse childhood experience phenotypes in the ROSMAP
single-nucleus datasets.

## Official ACE Phenotype Models

These are the formalized primary models for both cohorts:

| Model | Variable | Definition |
|-------|----------|------------|
| Total adversity | `tot_adverse_exp` | sum of the five ACE component scores |
| Early household SES | `early_hh_ses` | tracked continuous SES measure from the phenotype table |
| Aggregate score | `ace_aggregate` | `zscore(tot_adverse_exp) + zscore(early_hh_ses)` |

ACE component columns are also available for exploratory composition models:

- `emotional_neglect`
- `family_pro_sep`
- `financial_need`
- `parental_intimidation`
- `parental_violence`

## Shared Input Contract

### Processed AnnData

- `cell_type`
- Tsai: `sample_id` or `projid`
- DeJager: `patient_id`

### Phenotype CSV

Default: `${ACE_SCORES_CSV}`

Required columns:

- `projid`
- `tot_adverse_exp`
- `early_hh_ses`
- `msex`
- `age_death`
- `pmi`
- `niareagansc`
- the five ACE component columns above

## Supported ACE Workflows

| Workflow | Tsai | DeJager |
|----------|------|---------|
| DEG | `DEG/Tsai/aceDegT.sh` | `DEG/DeJager/aceDegDJ.sh` |
| Cell-type proportion | `CellTypeProportion/Tsai/acePropT.sh` | `CellTypeProportion/DeJager/acePropDJ.sh` |

Canonical integration choices:

- **Tsai**: `derived_batch`
- **DeJager**: `library_id`

Sensitivity integrations can still be passed explicitly to the launchers.

## Output Layout

All ACE outputs are written beneath:

```text
${ANALYSIS_OUTPUT_ROOT}/ACE/
  DEG/
    Tsai/
    DeJager/
  CellTypeProportion/
    Tsai/
    DeJager/
```

The repo directories store code only.

## Validation Layers

- Smoke: fixture-based checks that exercise the ACE prep and analysis entrypoints
  without requiring the full processed cohort objects
- Full: production-readiness checks against the canonical annotated cohort inputs

Preflight entrypoints:

- `bash config/preflight.sh ace-tsai-smoke`
- `bash config/preflight.sh ace-dejager-smoke`
- `bash config/preflight.sh ace-tsai-full`
- `bash config/preflight.sh ace-dejager-full`
