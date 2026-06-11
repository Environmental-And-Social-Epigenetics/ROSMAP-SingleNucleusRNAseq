# HOWTO — Start a New Phenotype Analysis

This directory (`Analysis/_template/`) is a copy-paste starting point for
bootstrapping a **new phenotype analysis**. Follow the steps below; replace
every `{Phenotype}` / `TODO` placeholder with your real values.

A complete, working example to model your code on is the ACE analysis:
`Analysis/ACE/` (in particular the DEG launcher `Analysis/ACE/DEG/Tsai/run_deg.sh`
and its smoke test `Analysis/ACE/DEG/Tsai/smoke_test.sh`).

---

## 1. Copy the template

From the repo root:

```bash
cp -r Analysis/_template Analysis/<NewPhenotype>
# e.g. cp -r Analysis/_template Analysis/Resilient
```

Then fill in `Analysis/<NewPhenotype>/README.md` (phenotype definition, patient
selection, key comparisons, phenotype data source).

## 2. Standard layout

Each analysis follows `<AnalysisType>/<Cohort>/`:

```
Analysis/<NewPhenotype>/
├── README.md                  Phenotype definition, patient selection, comparisons
├── DEG/                       Differential expression
│   ├── Tsai/                  Tsai cohort
│   └── DeJager/               DeJager cohort
├── SCENIC/                    Regulatory network inference (pySCENIC)
│   ├── Tsai/
│   └── DeJager/
└── TF/                        Transcription factor / activity analysis
    ├── Tsai/
    └── DeJager/
```

`<AnalysisType>` is one of `DEG`, `SCENIC`, `TF`, ... (add more as needed).
`<Cohort>` is `Tsai` or `DeJager`. Keep launchers and analysis scripts inside the
leaf cohort directory (e.g. `Analysis/<NewPhenotype>/DEG/Tsai/`).

## 3. The launcher pattern

Every launcher follows the same canonical structure (see
`DEG/Tsai/run_analysis.sh.template` in this template, modeled on
`Analysis/ACE/DEG/Tsai/run_deg.sh`):

1. **Resolve the repo root robustly.** Under SLURM, `BASH_SOURCE` points to a
   temp copy in `/var/spool/slurmd/`, so prefer `SLURM_SUBMIT_DIR` when set.
2. **`source config/paths.sh`** — this exports `ANALYSIS_OUTPUT_ROOT`, the
   per-analysis output roots (`<PHENO>_OUTPUT_ROOT`), the env-path variables
   (`NEBULA_ENV`, `DEG_ANALYSIS_ENV`, ...), and the `activate_env` /
   `init_conda` helper functions.
3. **`activate_env "${SOME_ENV}"`** — activates the conda env. `activate_env`
   safely wraps conda activation in `set +u` / restore, so you don't need your
   own nounset dance.
4. **`sbatch`** the launcher (or run it directly for a quick test).

Minimal sketch:

```bash
source config/paths.sh
activate_env "${DEG_ANALYSIS_ENV}"   # or NEBULA_ENV, SCENIC_ANALYSIS_ENV, ...
# ... run your Rscript / python here ...
```

To submit on the cluster:

```bash
cd Analysis/<NewPhenotype>/DEG/Tsai
sbatch run_analysis.sh <args...>
```

## 4. Where outputs go

**Never write generated outputs into the repo.** All results, intermediate
files, and figures go under `ANALYSIS_OUTPUT_ROOT` (defined in
`config/paths.sh`, default `${WORKSPACE_ROOT}/Analysis_Outputs`):

```
${ANALYSIS_OUTPUT_ROOT}/<Phenotype>/<AnalysisType>/<Cohort>/...
# e.g. ${ANALYSIS_OUTPUT_ROOT}/Resilient/DEG/Tsai/results_derived_batch/...
```

If you want a dedicated per-phenotype root variable (like `ACE_OUTPUT_ROOT` or
`SOCISL_OUTPUT_ROOT`), add one to `config/paths.sh`:

```sh
export <PHENO>_OUTPUT_ROOT="${<PHENO>_OUTPUT_ROOT:-${ANALYSIS_OUTPUT_ROOT}/<Phenotype>}"
```

Only small, hand-curated result CSVs and source code belong in the repo. Figures
and large matrices (`.h5ad`, `.rds`, `.feather`, ...) stay under the output root
and are `.gitignore`d.

## 5. Add a smoke test

Each launcher should have a companion `smoke_test.sh` that runs the analysis
end-to-end on a tiny fixture and asserts an output file exists. Model it on
`Analysis/ACE/DEG/Tsai/smoke_test.sh`. The pattern:

1. `source config/paths.sh`.
2. Redirect `ANALYSIS_OUTPUT_ROOT` to a throwaway `.../Smoke/...` directory so
   the test never pollutes real outputs.
3. Generate (or point at) a small fixture input.
4. `activate_env` and run the launcher in `--smoke` mode (small, fast).
5. `test -f <expected_output>` and `echo "<...> smoke test passed"`.

Run it with `bash Analysis/<NewPhenotype>/DEG/Tsai/smoke_test.sh` before
submitting full jobs.

## 6. Sanity checks before committing

```bash
bash -n Analysis/<NewPhenotype>/DEG/Tsai/run_analysis.sh   # shell syntax check
git status                                                  # confirm no outputs staged
```
