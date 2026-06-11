# ACE Microglial State Analysis

Tests whether microglia in ACE-exposed individuals are shifted along the
homeostatic-to-disease-associated (DAM) continuum, even though total
microglial numbers are unchanged (null sccomp result).

## Method

Microglia are scored for established state signatures (homeostatic, DAM stage 1,
DAM stage 2, interferon-response, inflammatory) using scanpy.tl.score_genes().
State scores are aggregated to the patient level and tested against ACE
phenotypes via OLS regression with covariates. Diffusion maps capture
continuous state variation.

## Cohort Directories

| Cohort | Entry Point |
|--------|-------------|
| Tsai | `Tsai/aceMicStateT.sh` |
| DeJager | `DeJager/aceMicStateDJ.sh` |

## Biological Motivation

The DEG analysis revealed 327 downregulated genes in male microglia. Cell-type
proportions showed no significant changes. Together, this suggests microglia
are functionally impaired without being lost. This analysis tests whether they
are shifted toward a disease-associated activation state.
