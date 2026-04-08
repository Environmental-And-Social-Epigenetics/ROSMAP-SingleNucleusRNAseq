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

Key scripts:

- `prep_counts.py`
- `run_prep.sh`
- `run_sccomp.sh`
- `run_visualize.sh`
- `smoke_test.sh`
- `sccomp_analysis.R`
- `sccomp_visualize.R`
