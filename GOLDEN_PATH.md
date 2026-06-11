# GOLDEN_PATH.md — The single linear tutorial

This is the one document to read top-to-bottom if you are new to the ROSMAP
snRNA-seq repository. It stitches the scattered per-area READMEs into a single
ordered path so you do not have to assemble the run yourself. Every command
below is copy-pasteable and exists in the repo.

If you only have 30 seconds: most newcomers should do
[Configure](#2-configure-5-minutes) →
[Install environments](#3-install-environments) →
[Smoke path](#4-smoke-path-no-lab-data-needed). That verifies your whole setup
end-to-end on generated fixtures without any protected data.

All commands assume you start from the repo root:

```bash
cd /orcd/data/lhtsai/001/mabdel03/ROSMAP_Code/Transcriptomics
```

---

## 1. Before you start — prerequisites and an honest scope note

**Honest scope note.** The FULL pipeline (preprocessing → processing → analysis
on real ROSMAP data) is a large, gated, multi-day undertaking. It requires:

- **Protected ROSMAP snRNA-seq data access** — the raw FASTQs (Tsai on the
  cluster, DeJager from Synapse) and downstream count matrices are gated and are
  **not** in this repo.
- **Cell Ranger 8.0.0 + the GRCh38-2020-A reference** — installed separately
  from 10x Genomics (see step 1 prerequisites below).
- **A SLURM cluster with GPUs** — preprocessing alone runs Cell Ranger for ~2
  days per batch and GPU CellBender; processing and analysis are submitted as
  SLURM array jobs.
- **Multi-day compute and large RAM** — e.g. integration and DEG steps request
  64–500 GB and many hours.

**If you do NOT have lab data**, do not try the full run. Instead do
[Configure](#2-configure-5-minutes), [Install environments](#3-install-environments),
and the [Smoke path](#4-smoke-path-no-lab-data-needed). The smoke tests run on
generated fixtures and verify your environment end-to-end without any protected
input.

### External prerequisites

| Prerequisite | Needed for | How to get it |
|--------------|-----------|---------------|
| **conda / miniforge** | everything (all envs) | Install miniforge; point `CONDA_INIT_SCRIPT` at its `conda.sh`. |
| **Cell Ranger 8.0.0** | preprocessing (alignment) | <https://www.10xgenomics.com/support/software/cell-ranger>; set `CELLRANGER_PATH`. |
| **GRCh38-2020-A reference** | preprocessing (alignment) | <https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz>; set `CELLRANGER_REF`. |
| **Protected ROSMAP data access** | the full run on real data | Gated; lab-internal. Not required for the smoke path. |
| **SCENIC ranking databases** (~3.5 GB) | optional — ACE/SocIsl SCENIC | <https://resources.aertslab.org/cistarget/>; set `SCENIC_RANKING_DIR`. |
| **IBM CPLEX (academic license)** | optional — COMPASS metabolic flux | <https://www.ibm.com/academic/>; set `CPLEX_DIR`. |
| **WGS / Demuxlet inputs + Singularity** | optional — DeJager demultiplexing only | Lab-internal; set `DEJAGER_WGS_DIR`. |

---

## 2. Configure (5 minutes)

Clone (skip if you already have the repo), then create your local path overrides
from the template and validate. `config/paths.sh` is the single source of truth
for every path; you customize it only through the gitignored
`config/paths.local.sh`.

```bash
git clone <repo-url> Transcriptomics
cd Transcriptomics

cp config/paths.local.sh.template config/paths.local.sh
$EDITOR config/paths.local.sh
```

In `config/paths.local.sh`, set at minimum:

- `CONDA_INIT_SCRIPT` — your conda/miniforge `etc/profile.d/conda.sh`
- `CONDA_ENV_BASE` — where conda envs live (or will be created)
- `CELLRANGER_PATH` — Cell Ranger 8.0.0 install dir (starts as `__UNCONFIGURED__`)
- `CELLRANGER_REF` — GRCh38-2020-A reference dir (starts as `__UNCONFIGURED__`)

The template includes a ready-to-copy **MIT Engaging** block near the bottom with
correct `/orcd/...` paths for `CONDA_INIT_SCRIPT`, `CONDA_ENV_BASE`, `DATA_ROOT`,
`ANALYSIS_OUTPUT_ROOT`, `CELLRANGER_PATH`, `CELLRANGER_REF`, and partitions —
start from that if you are on Engaging.

Then source the config and validate:

```bash
source config/paths.sh
check_paths
```

**Expected behavior:** `check_paths` will **ERROR** on `CELLRANGER_PATH` and
`CELLRANGER_REF` while they still hold the `__UNCONFIGURED__` sentinel value.
That is by design — the sentinel is an invalid path so an unconfigured clone
fails loudly. Fill in those two variables and re-run until you see
`All paths verified.` (or `PASSED with N warning(s).`). The optional sentinels
(`DEJAGER_WGS_DIR`, `CPLEX_DIR`, `SCENIC_RANKING_DIR`) print a `NOTE:` instead of
an error — they are only needed for specific optional workflows.

Validate the repo contract (env specs, variant definitions, tracked files):

```bash
PYTHONPATH=src python -m rosmap_tx.validate envs variants tracked-files
```

---

## 3. Install environments

The installer creates and **validates** each conda environment from the tracked
specs under `envs/`. Run it once, after `config/paths.sh` is sourced.

```bash
source config/paths.sh
bash setup/install_envs.sh --all
```

`--all` installs processing + preprocessing + analysis envs. Install only the
group you need with the per-group flags:

```bash
bash setup/install_envs.sh                  # processing envs only (default)
bash setup/install_envs.sh --preprocessing  # processing + preprocessing
bash setup/install_envs.sh --analysis       # processing + analysis
bash setup/install_envs.sh --all            # everything
bash setup/install_envs.sh --analysis --recreate   # rebuild stale envs
```

For the [smoke path](#4-smoke-path-no-lab-data-needed) you only need the analysis
group (`--analysis`), which builds `NEBULA_ENV` and `SCCOMP_ENV`. Choose an
install method explicitly if needed: `--method=conda` (default, from
`environment.yml`) or `--method=requirements` (bootstrap then
`pip install -r requirements.txt`). The installer prints the `export` lines to
copy into `config/paths.local.sh` if your env names differ from the defaults.

---

## 4. Smoke path (no lab data needed)

This is the recommended path for anyone without protected lab data. The ACE
smoke tests generate fixtures and run the real prep + analysis entrypoints
end-to-end, so a successful run proves your config and environments work.

Run the Tsai ACE DEG and cell-type-proportion smoke tests:

```bash
bash Analysis/ACE/DEG/Tsai/smoke_test.sh
bash Analysis/ACE/CellTypeProportion/Tsai/smoke_test.sh
```

DeJager equivalents:

```bash
bash Analysis/ACE/DEG/DeJager/smoke_test.sh
bash Analysis/ACE/CellTypeProportion/DeJager/smoke_test.sh
```

Each script calls `Analysis/ACE/_smoke/generate_fixtures.py` to synthesize a tiny
annotated cohort + phenotype CSV under `${ANALYSIS_OUTPUT_ROOT}/ACE/Smoke/...`,
then runs the actual launchers in `--smoke` mode and asserts the summary CSV was
written. On success you will see e.g. `Tsai ACE DEG smoke test passed`. These run
in minutes on the analysis envs and never touch real data.

You can also run the preflight smoke checks:

```bash
bash config/preflight.sh env-specs
bash config/preflight.sh ace-tsai-smoke
```

---

## 5. Full run (with data)

Only attempt this with protected data access, Cell Ranger, and SLURM. The flow
is linear: **Preprocessing → Processing → Analysis**. Run preprocessing once per
cohort you need.

### 5a. Preprocessing — FASTQ → CellBender counts

Follow the per-cohort README; these are multi-day SLURM jobs.

- **Tsai** — [Preprocessing/Tsai/README.md](Preprocessing/Tsai/README.md).
  Locate FASTQs → Cell Ranger `count` → GPU CellBender for all 480 patients (16
  batches of 30). Entry point:
  ```bash
  cd Preprocessing/Tsai/02_Cellranger_Counts
  python Scripts/generate_batch_scripts.py
  sbatch Scripts/pipeline_slurm_wrapper.sh
  ```
  **What you get:** CellBender-filtered `.h5` per `projid` under
  `${TSAI_PREPROCESSED}`. **Runtime/RAM:** Cell Ranger ~2 days, 64 GB, 16 cores;
  CellBender ~4 h, 64 GB, 1 GPU (per batch).

- **DeJager** — [Preprocessing/DeJager/README.md](Preprocessing/DeJager/README.md).
  Synapse download → Cell Ranger → CellBender → Demuxlet/Freemuxlet
  demultiplexing (requires `DEJAGER_WGS_DIR` + Singularity).
  **What you get:** CellBender counts plus per-cell patient assignments at
  `${DEJAGER_PATIENT_MAP}`. **Runtime/RAM:** Cell Ranger ~47 h, 128 GB;
  CellBender ~47 h, 128 GB, A100; Demuxlet pileup ~36 h, 500 GB.

### 5b. Processing — three-stage QC / doublets / integration

Both cohorts share the same three-stage pipeline (Stage 1 percentile QC →
Stage 2 `scDblFinder` doublet removal → Stage 3 normalization + Harmony +
clustering + ORA annotation). The submit wrapper chains the stages with SLURM
dependencies. See [Processing/Tsai/Pipeline/README.md](Processing/Tsai/Pipeline/README.md).

```bash
cd Processing/Tsai/Pipeline
./submit_pipeline.sh all
```

Run a named pipeline variant instead of the primary (outputs land in
variant-namespaced dirs):

```bash
./submit_pipeline.sh --variant <id> all       # see config/variants.yaml / docs/VARIANTS.md
```

DeJager uses the same wrapper under `Processing/DeJager/Pipeline/`.

**What you get:** per-stage outputs under `${TSAI_PROCESSING_OUTPUTS}`, ending in
`03_Integrated/<variant>/tsai_annotated.h5ad` (the annotated object every
analysis consumes). **Runtime/RAM:** Stage 1/2 run as a `--array=1-478%32` SLURM
array; Stage 3 is a single large integration job. Each successful launcher run
writes a `run_manifest.json` beside its output directory (see step 7).

### 5c. Analysis — phenotype DEG / proportions / etc.

Downstream analyses build on the annotated object. The ACE DEG launcher takes a
**variant** (use `primary` to resolve to the declared primary) and a
**phenotype** (`tot_adverse_exp`, `early_hh_ses`, or `ace_aggregate`). See
[Analysis/README.md](Analysis/README.md) and [Analysis/ACE/README.md](Analysis/ACE/README.md).

```bash
# pseudobulk DEG, primary integration, total-adversity phenotype:
cd Analysis/ACE/DEG/Tsai
sbatch run_deg.sh primary tot_adverse_exp
```

The canonical one-shot launcher `bash Analysis/ACE/DEG/Tsai/aceDegT.sh` builds
the per-cell-type splits and runs the official phenotype trio.
**What you get:** DEG result tables + volcano plots under
`${ANALYSIS_OUTPUT_ROOT}/ACE/DEG/Tsai/results_<variant>/<phenotype>/`.
**Runtime/RAM:** `run_deg.sh` requests 200 GB, 4 cores, up to 24 h.

Cell-type proportion (sccomp) follows the same pattern:

```bash
cd Analysis/ACE/CellTypeProportion/Tsai
sbatch run_prep.sh        # build per-cohort counts
sbatch run_sccomp.sh      # sccomp compositional model (8 cores, 64 GB)
```

GSEA, TF Activity, SCENIC, MicState launchers are listed per workflow in
[Analysis/ACE/README.md](Analysis/ACE/README.md).

---

## 6. Choosing / switching the pipeline variant

A **variant** is a named version of the processing pipeline (different QC
thresholds, Harmony key, HVG count, etc.). One variant is declared the
**primary** per dataset, and every default run plus all downstream analysis
resolve to it. To switch the primary, edit the single `primary:` block at the top
of `config/variants.yaml` (e.g. change `tsai: "derived_batch"`), then validate
with `PYTHONPATH=src python -m rosmap_tx.validate variants`. After that,
`./submit_pipeline.sh all` and `run_deg.sh primary <phenotype>` automatically use
the new primary. Full details, the variant schema, and output-namespacing rules
are in [docs/VARIANTS.md](docs/VARIANTS.md).

---

## 7. Where results go and how to read them

Generated outputs never go inside the repo tree. Two roots matter:

- **Processing outputs** → `${TSAI_PROCESSING_OUTPUTS}` /
  `${DEJAGER_PROCESSING_OUTPUTS}` (`01_QC_Filtered/`, `02_Doublet_Removed/`,
  `03_Integrated/<variant>/`, `Logs/`).
- **Analysis outputs** → `${ANALYSIS_OUTPUT_ROOT}` (defaults to a sibling
  `Analysis_Outputs/` next to the repo). ACE results land under
  `${ANALYSIS_OUTPUT_ROOT}/ACE/<Workflow>/<Cohort>/` — DEG tables, sccomp result
  objects, GSEA pathways, figures, and SLURM logs.

**Provenance:** every successful processing launcher run writes a
`run_manifest.json` beside the output directory recording the dataset, stage,
variant, and resolved inputs — read it to know exactly which configuration
produced a given output. Figures land in the per-stage `figures/` subdirectories
(e.g. `03_Integrated/<variant>/figures/` for UMAPs); result tables land beside
them as `.csv`. The data namespace and copy/symlink modes are specified in
[docs/DATA_CONTRACT.md](docs/DATA_CONTRACT.md).

---

## 8. Troubleshooting / next reads

- [KNOWN_ISSUES.md](KNOWN_ISSUES.md) — flagged problems and gotchas to check
  first when something fails.
- [docs/RUNBOOK.md](docs/RUNBOOK.md) — the terse operator command reference
  (this tutorial complements it; the runbook is the cheat sheet).
- [docs/DATA_CONTRACT.md](docs/DATA_CONTRACT.md) — the canonical data namespace
  and how to point it at external storage.
- [docs/VARIANTS.md](docs/VARIANTS.md) — pipeline variants in depth.
- [setup/README.md](setup/README.md) — the full first-time setup contract and
  environment inventory.
