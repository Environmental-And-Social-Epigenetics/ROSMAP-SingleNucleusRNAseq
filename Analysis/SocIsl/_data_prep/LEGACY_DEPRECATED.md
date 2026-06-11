# LEGACY / DEPRECATED — SocIsl `_data_prep/` scripts

**Status: DEPRECATED — historical reference only. Do NOT use for new work.**

The ~25 scripts under `Tsai/` and `DeJager/` in this directory were the original
ad-hoc preprocessing and data-preparation workflow used during the Social
Isolation analysis on the MIT Openmind cluster (2024-2025). They are retained
in the repository purely as a historical record of how the original results were
produced.

## Why they are deprecated

- **The `Processing/` pipeline supersedes them.** All QC filtering, doublet
  removal, integration, and cell-type annotation are now handled by the
  standardized pipeline (`Processing/*/Pipeline/`, e.g.
  `Processing/Tsai/Pipeline/submit_pipeline.sh all`), which produces the
  annotated `tsai_annotated.h5ad` / `dejager_annotated.h5ad` inputs that all
  current Analysis code consumes. If you have those annotated objects, you do
  **not** need anything in this directory.
- **They contain hardcoded Openmind paths.** Every script references
  decommissioned Openmind directories (e.g. `/om2/user/mabdel03/...`,
  `/om/scratch/.../SocialIsolation/`, `/net/vast-storage/...`) and Openmind
  conda activation (`/om2/user/mabdel03/anaconda/etc/profile.d/conda.sh`).
  **They will not run on Engaging without manual path updates.** See the
  "Updating Hardcoded Paths" table in `../README.md` for the old → new mapping.

## What to do instead

For any new phenotype or re-analysis work:

1. Run the `Processing/` pipeline to generate the annotated h5ad inputs.
2. Use the current Analysis launchers (which `source config/paths.sh`, call
   `activate_env`, and write to `${ANALYSIS_OUTPUT_ROOT}`), modeled on
   `Analysis/ACE/DEG/Tsai/run_deg.sh` and the skeleton in `Analysis/_template/`.

Each script in this directory already carries a one-line `DEPRECATED` header
comment. Treat them as read-only historical artifacts; do not extend them.
