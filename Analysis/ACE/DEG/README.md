# ACE DEG

Pseudobulk differential expression for the ACE phenotype set.

## Primary Models

- `tot_adverse_exp`
- `early_hh_ses`
- `ace_aggregate`

The DEG workflow uses the shared ACE phenotype table and writes split inputs,
result tables, and logs under `ANALYSIS_OUTPUT_ROOT/ACE/DEG/`.

Fixture-based smoke tests are available in each cohort directory as
`smoke_test.sh`.

## Cohort Entry Points

- [Tsai/README.md](Tsai/README.md)
- [DeJager/README.md](DeJager/README.md)
