# Runbook

## Setup

```bash
cp config/paths.local.sh.template config/paths.local.sh
$EDITOR config/paths.local.sh
bash setup/install_envs.sh --all
python -m pip install -e .
```

## Validate

```bash
python -m rosmap_tx.validate paths data envs variants tracked-files
bash config/preflight.sh env-specs
```

## Processing

Existing SLURM entrypoints remain supported:

```bash
Processing/Tsai/Pipeline/submit_pipeline.sh all
Processing/DeJager/Pipeline/submit_pipeline.sh all
```

The shared launcher can also run individual stages from the appropriate environment:

```bash
python -m rosmap_tx.processing --dataset tsai --stage 1 --variant canonical
python -m rosmap_tx.processing --dataset dejager --stage 3 --variant library_id
python -m rosmap_tx.processing --dataset dejager --stage 3 --variant no_harmony
```

Each successful launcher run writes `run_manifest.json` beside the output directory.

