# ACE Transcription Factor & Pathway Activity

Decoupler-based TF activity (DoRothEA, CollecTRI) and pathway activity
(PROGENy) analysis for the ROSMAP ACE phenotype study.

## Cohorts

| Cohort | Directory | Status |
|--------|-----------|--------|
| Tsai | `Tsai/` | Active |
| DeJager | `DeJager/` | Planned |

## Quick Start

```bash
# Tsai cohort — submit all phenotype jobs
bash Tsai/aceTfActT.sh

# Smoke test (no SLURM, uses fixtures)
bash Tsai/smoke_test.sh
```

## Approach

Uses [decoupler](https://decoupler-py.readthedocs.io/) to infer
transcription factor and pathway activities from pseudobulked snRNA-seq
data.  Activities are then tested for association with ACE phenotypes via
OLS regression with covariates (age_death, pmi, niareagansc), stratified
by sex.

Three complementary databases provide different coverage/confidence
trade-offs:

- **DoRothEA** (levels A-C): high-confidence curated TF regulons
- **CollecTRI**: broader community-collected TF regulons
- **PROGENy** (top 500): pathway-responsive gene signatures for 14
  canonical signaling pathways

## Relation to SCENIC

This analysis complements the SCENIC analysis (`../SCENIC/`).  SCENIC
uses GRNBoost2 + cisTarget to discover TF regulons de novo from the data,
then scores them with AUCell.  The decoupler approach here uses curated
prior-knowledge regulons and a different scoring method (MLM).  The
convergence script (`tf_convergence_analysis.py`) can cross-reference
results from both methods.

## Output Layout

All outputs go to `${ACE_OUTPUT_ROOT}/TFActivity/<Cohort>/`.
See `Tsai/README.md` for detailed output structure.
