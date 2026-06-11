# Tsai ACE SCENIC Legacy

DEPRECATED: these scripts are superseded by the modular `scenic_analysis.py`
pipeline and are retained for reference only.

| File | Was | Now |
|---|---|---|
| `aceScenic.py` | Monolithic Tsai SCENIC script (hardcoded cell types, paths, plotting inline) | Replaced by `../scenic_analysis.py` (cohort-aware, one cell-type x sex per call, `--pool-size` default 50) |
| `aceScenic.sh` | Thin SLURM wrapper that ran `aceScenic.py` | Replaced by `../aceScenicT.sh` -> `../run_scenic.sh` |

Do not run these. The supported entry point is:

```bash
source config/paths.sh
bash Analysis/ACE/SCENIC/Tsai/aceScenicT.sh [INTEGRATION] [PHENOTYPE]
```

These files are not covered by the smoke test and reflect older path/output
assumptions. After relocation into `legacy/`, their internal `REPO_ROOT` and
`aceScenic.py` path derivations are no longer correct; both have an early
`exit 1` / deprecation banner to prevent accidental execution.
