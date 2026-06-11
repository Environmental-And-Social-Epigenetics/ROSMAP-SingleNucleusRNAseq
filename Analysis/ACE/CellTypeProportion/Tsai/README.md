# ACE Cell-Type Proportion - Tsai

Canonical entry point:

```bash
bash acePropT.sh
```

By default this runs the Tsai sccomp workflow on the canonical
`derived_batch`-integrated object and writes outputs to:

```text
${ANALYSIS_OUTPUT_ROOT}/ACE/CellTypeProportion/Tsai
```

Optional sensitivity integrations can be passed as positional arguments:

```bash
bash acePropT.sh derived_batch projid
```

## Shared sccomp implementation (Tsai == DeJager)

`Tsai/sccomp_analysis.R` is the **single source of truth** for the ACE
cell-type proportion analysis. Both cohorts call the identical statistical
core: the modern, maintained sccomp pipeline

```text
sccomp_estimate() |> sccomp_remove_outliers() |> sccomp_test()
```

plus `sccomp_proportional_fold_change()`. The deprecated monolithic
`sccomp_glm()` is no longer used.

The DeJager `sccomp_analysis.R` is a thin wrapper that `source()`s this file
(mirroring the DEG parity pattern where `aceDegDJ.Rscript` source()s
`aceDegT.Rscript`). The script is cohort-aware via:

- `--cohort {tsai,dejager}` (or `ACE_PROP_COHORT`) — selects the default
  integration label (`tsai` → `derived_batch`, `dejager` → `library_id`);
- `--output-root` (or `ACE_PROP_OUTPUT_ROOT`) — selects the cohort output tree,
  standardized to `${ACE_OUTPUT_ROOT}/CellTypeProportion/<cohort>`.

The model formula, covariates, priors (`bimodal_mean_variability_association`,
`mcmc_seed = 42`), smoke-test mode, and sex stratification are identical across
cohorts. The only cohort-specific input is the cell-grouping id column, which is
normalized to `sample`/`projid` upstream by each cohort's `prep_counts.py`
(Tsai keys on `projid`/`sample_id`, DeJager on `patient_id`), so the R core
reads identical columns for both.

Key scripts:

- `prep_counts.py`
- `run_prep.sh`
- `run_sccomp.sh`
- `run_visualize.sh`
- `smoke_test.sh`
- `sccomp_analysis.R` (canonical, cohort-aware)
- `sccomp_visualize.R`
