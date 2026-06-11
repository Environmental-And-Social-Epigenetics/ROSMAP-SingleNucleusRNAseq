# ACE Microglial State Analysis - Tsai

Canonical entry point:

```bash
bash aceMicStateT.sh
```

## Method

1. Load microglia h5ad from cell type splits
2. Score cells for established state gene signatures:
   - **Homeostatic**: P2RY12, CX3CR1, TMEM119, CSF1R, HEXB, SELPLG
   - **DAM stage 1**: TYROBP, APOE, B2M, FTH1 (Trem2-independent)
   - **DAM stage 2**: TREM2, AXL, CST7, LPL, SPP1, ITGAX (Trem2-dependent)
   - **Interferon-response**: ISG15, IFI44L, IFIT1, MX1, OAS1
   - **Inflammatory**: IL1B, TNF, CCL2, CCL3, CXCL10
3. Compute diffusion map for continuous state visualization
4. Aggregate state scores to patient level
5. OLS regression: `state_score ~ tot_adverse_exp + age_death + pmi + niareagansc`
6. K-means clustering on state scores (k=3-5) for state proportions
7. Test ACE association with state proportions

## Prerequisites

- Microglia h5ad: `${ACE_OUTPUT_ROOT}/DEG/Tsai/celltype_splits_${INTEGRATION}/Mic.h5ad`
- ACE phenotype CSV: `${ACE_SCORES_CSV}`
- `${DECOUPLER_ENV}` or `${NEBULA_ENV}` with scanpy

## Output

Results written to `${ACE_OUTPUT_ROOT}/MicState/Tsai/results_${INTEGRATION}/`.

```
results_derived_batch/
├── Male/
│   ├── state_scores_per_cell.csv
│   ├── state_scores_per_patient.csv
│   ├── state_regression_results.csv
│   ├── state_cluster_proportions.csv
│   ├── state_proportion_test.csv
│   └── diffmap_coordinates.csv
└── Female/
    └── ...
```

## Interpretation

- Higher homeostatic score = more surveilling, protective state
- Higher DAM score = disease-associated, potentially dysfunctional
- ACE shifting microglia toward DAM without changing cell numbers supports
  the "functional impairment without cell loss" hypothesis

Canonical integration: `derived_batch`
