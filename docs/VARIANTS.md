# Pipeline Variants

A **variant** is a named version of the processing pipeline. Variants let you
experiment — different QC thresholds, doublet method, integration parameters,
Harmony key — while keeping ONE declared **primary** pipeline that every default
run and all downstream analysis resolve to. Everything lives in one git branch;
you switch the primary by editing a single line.

All variant IDs are defined in [`config/variants.yaml`](../config/variants.yaml)
and are the only supported names.

## Switching the primary pipeline

Edit the one `primary:` block at the top of `config/variants.yaml`:

```yaml
primary:
  tsai: "derived_batch"     # change this to switch the Tsai primary
  dejager: "library_id"     # change this to switch the DeJager primary
```

After changing it:

- `./submit_pipeline.sh all` (no `--variant`) runs the new primary.
- `python -m rosmap_tx.processing --variant primary ...` resolves to the new primary.
- `run_deg.sh primary <phenotype>` (and the other analysis launchers) use the new primary.
- `python -m rosmap_tx.config --primary tsai` prints the current primary id.

Validate the change with `PYTHONPATH=src python -m rosmap_tx.validate variants`.

## Vocabulary

| Term | Meaning |
|------|---------|
| `primary` | The single chosen variant per dataset. The owner-facing word. |
| `canonical` / `main` / `official` | Input aliases that resolve to `primary`. |
| `role:` | A label on each variant (`primary` / `sensitivity` / `experiment`); informational only. |

`status.official` in `config/analysis_models.yaml` is a different, unrelated
concept (which analysis *models* are production-grade) — it is not a variant.

## Variant schema

Each entry under `datasets.<ds>.variants.<id>` may declare:

| Field | Purpose |
|-------|---------|
| `role` | `primary` / `sensitivity` / `experiment` (label only). |
| `harmony_batch_key` | obs column used for Harmony (Stage 3). |
| `skip_harmony` | `true` to use the PCA embedding directly (no correction). |
| `output_env` | env var (from `config/paths.sh`) holding the Stage 3 dir. |
| `overrides` | optional subtree mirroring `config/pipeline.yaml`, deep-merged over the defaults so a variant can change ANY stage parameter. |
| `touches` | list of stages this variant changes vs the primary (e.g. `[stage1, stage3]`). Defaults to `[stage3]`. |
| `analysis` | optional per-workflow model gating, e.g. `analysis: { ace: { models: [sex_stratified_ace] } }`. |

## Output namespacing (how outputs never collide)

Each stage writes to a per-variant **leaf** directory:

- Stage 3 is always per-variant: `03_Integrated/<variant_id>/`.
- Stages 1 and 2 use a `shared` leaf **unless** the variant changes that stage
  (or an upstream stage). Divergence is **sticky**: once a variant changes an
  earlier stage, every later stage is also namespaced, because its input
  changed. This means cheap variants (e.g. only a different HVG count) reuse the
  primary's Stage 1/2 outputs instead of recomputing them.

Worked example — `strict_qc` (`touches: [stage1, stage3]`):

| Stage | Reads | Writes |
|-------|-------|--------|
| 1 | `Cellbender_Output` | `01_QC_Filtered/strict_qc` |
| 2 | `01_QC_Filtered/strict_qc` | `02_Doublet_Removed/strict_qc` (sticky) |
| 3 | `02_Doublet_Removed/strict_qc` | `03_Integrated/strict_qc` |

Worked example — `hvg5000` (`touches: [stage3]`):

| Stage | Reads | Writes |
|-------|-------|--------|
| 1 | `Cellbender_Output` | `01_QC_Filtered/shared` |
| 2 | `01_QC_Filtered/shared` | `02_Doublet_Removed/shared` |
| 3 | `02_Doublet_Removed/shared` | `03_Integrated/hvg5000` |

## Running a variant

```bash
# the primary, end to end
cd Processing/Tsai/Pipeline && ./submit_pipeline.sh all

# a named variant, end to end (outputs land in variant-namespaced dirs)
./submit_pipeline.sh --variant strict_qc all

# inspect resolved dirs without running anything
PYTHONPATH=src python -m rosmap_tx.processing --dataset tsai --stage 3 \
  --variant strict_qc --print-output-dir
```

## Adding a new variant

1. Add an entry under `datasets.<ds>.variants` in `config/variants.yaml`.
2. If it changes Stage 1 or Stage 2 params, add those stages to `touches` (the
   validator hard-fails an `overrides` that is not listed in `touches`, since
   that would silently collide with the primary's `shared` output).
3. Run `PYTHONPATH=src python -m rosmap_tx.validate variants`.
4. Run it: `./submit_pipeline.sh --variant <id> all`.

## Current Tsai variants

| Variant | Role | Harmony key | Notes |
|---------|------|-------------|-------|
| `derived_batch` | primary | `derived_batch` | Flowcell-derived batch correction. |
| `projid` | sensitivity | `projid` | Earlier patient-level correction. |
| `no_harmony` | sensitivity | none | PCA baseline with Harmony skipped. |
| `strict_qc` | experiment | `derived_batch` | Stricter QC + 2000 HVGs (touches stage1, stage3). |
| `hvg5000` | experiment | `derived_batch` | 5000 HVGs / 50 PCs (touches stage3). |

## Current DeJager variants

| Variant | Role | Harmony key | Notes |
|---------|------|-------------|-------|
| `library_id` | primary | `library_id` | Per-library correction. |
| `patient_id` | sensitivity | `patient_id` | Patient-level correction. |
| `pool_batch` | sensitivity | `pool_batch` | B-number pool correction. |
| `derived_batch` | sensitivity | `derived_batch` | Flowcell-derived correction. |
| `sequencing_date` | sensitivity | `sequencing_date` | YYMMDD library-prefix correction. |
| `no_harmony` | sensitivity | none | PCA baseline with Harmony skipped. |
