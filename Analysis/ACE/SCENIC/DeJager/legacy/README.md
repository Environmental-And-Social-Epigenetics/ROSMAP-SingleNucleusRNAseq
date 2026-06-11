# DeJager ACE SCENIC Legacy

DEPRECATED: these scripts are superseded by the modular `scenic_analysis.py`
pipeline (shared with Tsai) and are retained for reference only.

| File | Was | Now |
|---|---|---|
| `aceScenic.py` | Monolithic DeJager SCENIC script with an INTERNALLY INCONSISTENT micropool `pool_size` (30 males / 39 females) | Replaced by the shared `../../Tsai/scenic_analysis.py` run with `--cohort dejager` and a SINGLE `--pool-size` (default 50, identical to Tsai) |
| `aceScenic.sh` | Thin SLURM wrapper that ran `aceScenic.py` | Replaced by `../aceScenicDJ.sh` -> `../run_scenic.sh` (delegates to the Tsai modular pipeline) |

Do not run these. The supported entry point is:

```bash
source config/paths.sh
bash Analysis/ACE/SCENIC/DeJager/aceScenicDJ.sh [INTEGRATION] [PHENOTYPE]
```

The pool_size inconsistency (30/39) is the single most important reason this
script is deprecated: the unified pipeline uses ONE pool_size for both cohorts
and both sexes, so the DeJager and Tsai cohorts now run identical methods.
After relocation into `legacy/`, the internal `REPO_ROOT`/`aceScenic.py` path
derivations are no longer correct; both files have an early `exit 1` /
deprecation banner to prevent accidental execution.
