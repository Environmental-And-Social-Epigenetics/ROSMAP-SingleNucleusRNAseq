# Repository Audit

Date: 2026-05-26

## Findings

- The active refactor branch is `5-26-2026_Refactor`; the clean pre-refactor branch is `5-26-2026_PRE-Refactor`.
- Source and generated outputs are mixed. The largest tracked issue is 1,026 generated Tsai Cell Ranger/CellBender batch scripts under `Preprocessing/Tsai/02_Cellranger_Counts/Batch_Scripts/`.
- Generated run products also appear as untracked `Project_*` directories, ACE report/result folders, `__pycache__`, logs, PDFs, and `ckpt.tar.gz`.
- Processing parameters are duplicated across Tsai and DeJager scripts instead of being declared once.
- Stage 3 integration variants are encoded through wrapper names and output directory names.
- DeJager Stage 3 wrappers were inconsistent about the validated environment (`BATCHCORR_ENV`) and the activated environment (`CONDA_ENV_BASE/decoupler_compat`).
- Some local or stale absolute paths remain in helper/report/legacy areas and should be converted to `config/paths.sh` variables before those workflows are promoted to production.

## Refactor Decisions

- Biological scripts remain the execution engines; shared Python code only handles config, launch command construction, validation, and manifests.
- Canonical data lives under `Data/Transcriptomics/{Tsai,DeJager}` or symlinks from there to external storage.
- Canonical processing variants live in `config/variants.yaml` and write to `Processing_Outputs/03_Integrated/{variant_id}`.
- Large data and generated products are ignored and removed from git tracking without deleting local files.

