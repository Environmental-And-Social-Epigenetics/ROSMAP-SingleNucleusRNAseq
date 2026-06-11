# ACE Cell-Type Proportion - DeJager

Canonical entry point:

```bash
bash acePropDJ.sh
```

By default this runs the DeJager sccomp workflow on the canonical
`library_id`-integrated object and writes outputs to:

```text
${ANALYSIS_OUTPUT_ROOT}/ACE/CellTypeProportion/DeJager
```

Optional sensitivity integrations can be passed as positional arguments:

```bash
bash acePropDJ.sh library_id patient_id pool_batch derived_batch
```

## Shared sccomp implementation (DeJager == Tsai)

DeJager and Tsai run the **identical** sccomp statistical core. The DeJager
`sccomp_analysis.R` is a thin wrapper that `source()`s the single source of
truth at `../Tsai/sccomp_analysis.R` (mirroring the DEG parity pattern where
`aceDegDJ.Rscript` source()s `aceDegT.Rscript`). It only sets the default
cohort to `dejager` (→ integration `library_id`) before sourcing.

Both cohorts therefore call the modern, maintained sccomp pipeline

```text
sccomp_estimate() |> sccomp_remove_outliers() |> sccomp_test()
```

plus `sccomp_proportional_fold_change()`. The deprecated monolithic
`sccomp_glm()` is no longer used.

Cohort-specific bits stay parameterized: the integration label (via `--cohort`
/ `ACE_PROP_COHORT`), the output tree (via `--output-root` /
`ACE_PROP_OUTPUT_ROOT`, standardized to
`${ACE_OUTPUT_ROOT}/CellTypeProportion/DeJager`), and the cell-grouping id
column, which `prep_counts.py` normalizes to `sample`/`projid` upstream
(DeJager keys on `patient_id`). The model formula, covariates, priors, smoke
mode, and sex stratification are identical across cohorts.

Key scripts:

- `prep_counts.py`
- `run_prep.sh`
- `run_sccomp.sh`
- `run_visualize.sh`
- `smoke_test.sh`
- `sccomp_analysis.R` (thin wrapper around `../Tsai/sccomp_analysis.R`)
- `sccomp_visualize.R`
