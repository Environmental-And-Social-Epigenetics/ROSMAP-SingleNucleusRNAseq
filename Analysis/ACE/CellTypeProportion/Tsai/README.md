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

Key scripts:

- `prep_counts.py`
- `run_prep.sh`
- `run_sccomp.sh`
- `run_visualize.sh`
- `smoke_test.sh`
- `sccomp_analysis.R`
- `sccomp_visualize.R`
