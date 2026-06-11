from __future__ import annotations

import copy
import hashlib
import os
import re
import shlex
import subprocess
from pathlib import Path
from typing import Any

import yaml


# Input aliases that all resolve to the dataset's declared primary variant.
PRIMARY_ALIASES = {"primary", "canonical", "main", "official"}


DATASET_DIRS = {
    "tsai": "Tsai",
    "dejager": "DeJager",
}


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def config_dir(root: Path | None = None) -> Path:
    return (root or repo_root()) / "config"


def load_yaml(name: str, root: Path | None = None) -> dict[str, Any]:
    path = config_dir(root) / name
    with path.open("r", encoding="utf-8") as handle:
        loaded = yaml.safe_load(handle) or {}
    if not isinstance(loaded, dict):
        raise TypeError(f"{path} must contain a YAML mapping")
    return loaded


def yaml_sha256(names: list[str] | None = None, root: Path | None = None) -> dict[str, str]:
    root = root or repo_root()
    names = names or [
        "datasets.yaml",
        "pipeline.yaml",
        "variants.yaml",
        "analysis_models.yaml",
    ]
    hashes: dict[str, str] = {}
    for name in names:
        path = config_dir(root) / name
        if path.exists():
            hashes[str(path.relative_to(root))] = hashlib.sha256(path.read_bytes()).hexdigest()
    return hashes


def paths_env(root: Path | None = None) -> dict[str, str]:
    """Return the environment produced by sourcing config/paths.sh."""
    root = root or repo_root()
    paths_sh = root / "config" / "paths.sh"
    command = f"source {shlex.quote(str(paths_sh))} >/dev/null && env -0"
    completed = subprocess.run(
        ["bash", "-lc", command],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=False,
    )
    env: dict[str, str] = {}
    for item in completed.stdout.split(b"\0"):
        if not item or b"=" not in item:
            continue
        key, value = item.split(b"=", 1)
        env[key.decode()] = value.decode()
    return env


_VAR_PATTERN = re.compile(r"\$\{([A-Za-z_][A-Za-z0-9_]*)\}|\$([A-Za-z_][A-Za-z0-9_]*)")


def expand_vars(value: str, env: dict[str, str]) -> str:
    def repl(match: re.Match[str]) -> str:
        key = match.group(1) or match.group(2)
        return env.get(key, match.group(0))

    return _VAR_PATTERN.sub(repl, value)


def dataset_dir(dataset: str) -> str:
    key = dataset.lower()
    if key not in DATASET_DIRS:
        valid = ", ".join(sorted(DATASET_DIRS))
        raise ValueError(f"Unknown dataset '{dataset}'. Valid datasets: {valid}")
    return DATASET_DIRS[key]


def dataset_prefix(dataset: str) -> str:
    return dataset.upper()


def deep_merge(base: dict[str, Any], override: dict[str, Any] | None) -> dict[str, Any]:
    """Recursively merge ``override`` into a deep copy of ``base``.

    Nested dicts are merged; scalars and lists from ``override`` replace the
    corresponding value in ``base``. ``base`` is never mutated.
    """
    out = copy.deepcopy(base)
    for key, value in (override or {}).items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = deep_merge(out[key], value)
        else:
            out[key] = copy.deepcopy(value)
    return out


def primary_variant_id(dataset: str, root: Path | None = None) -> str:
    """Return the dataset's declared primary variant id.

    Reads the top-level ``primary:`` block, falling back to the per-dataset
    ``canonical:`` key for backward compatibility.
    """
    variants = load_yaml("variants.yaml", root)
    dataset_key = dataset.lower()
    top = variants.get("primary", {})
    if dataset_key in top:
        return top[dataset_key]
    block = variants.get("datasets", {}).get(dataset_key)
    if not block:
        raise ValueError(f"No variants configured for dataset '{dataset}'")
    if "canonical" in block:
        return block["canonical"]
    raise ValueError(
        f"Dataset '{dataset}' has no primary defined (add it to the top-level "
        f"primary: block in variants.yaml)"
    )


def resolve_variant(dataset: str, requested: str, root: Path | None = None) -> tuple[str, dict[str, Any]]:
    variants = load_yaml("variants.yaml", root)
    dataset_key = dataset.lower()
    dataset_block = variants.get("datasets", {}).get(dataset_key)
    if not dataset_block:
        raise ValueError(f"No variants configured for dataset '{dataset}'")

    variant_id = requested
    if requested in PRIMARY_ALIASES:
        variant_id = primary_variant_id(dataset, root)

    records = dataset_block.get("variants", {})
    if variant_id not in records:
        valid = ", ".join(sorted(records))
        raise ValueError(f"Unknown {dataset} variant '{requested}'. Valid variants: primary, {valid}")
    return variant_id, records[variant_id]


def resolve_pipeline_variant(
    dataset: str, requested: str, stage: str, root: Path | None = None
) -> tuple[str, dict[str, Any], dict[str, Any]]:
    """Resolve a variant and return its deep-merged parameters for ``stage``.

    Returns ``(variant_id, record, merged_stage_params)`` where
    ``merged_stage_params`` is ``pipeline.yaml[stageN]`` deep-merged with the
    variant's ``overrides.stageN`` (the variant's values win).
    """
    variant_id, record = resolve_variant(dataset, requested, root)
    pipeline = load_yaml("pipeline.yaml", root)
    base = pipeline.get(f"stage{stage}", {})
    overrides = (record.get("overrides") or {}).get(f"stage{stage}", {})
    return variant_id, record, deep_merge(base, overrides)


def stage_leaf(record: dict[str, Any], stage: str, variant_id: str) -> str:
    """Output-dir leaf for a stage under a variant.

    Stage 3 is always per-variant. A stage gets its own leaf (the variant id) if
    the variant touches THIS stage OR any UPSTREAM stage -- divergence is sticky:
    once a variant changes an earlier stage, every later stage's output also
    differs (its input changed) and must be namespaced to avoid colliding with
    the primary's ``shared`` output. A stage upstream of every change shares the
    primary's output via the literal ``shared`` leaf.
    """
    stage_n = int(stage)
    if stage_n == 3:
        return variant_id
    touches = record.get("touches", ["stage3"])
    # touched at or before this stage -> diverged here or upstream -> own leaf
    for s in range(1, stage_n + 1):
        if f"stage{s}" in touches:
            return variant_id
    return "shared"


def stage_output_dir(dataset: str, stage: str, leaf: str, env: dict[str, str]) -> Path:
    """Per-stage output directory for stages 1 and 2 (Stage 3 uses
    ``variant_output_dir``). The configured base dir gets the ``leaf`` appended."""
    prefix = dataset_prefix(dataset)
    bases = {
        "1": env[f"{prefix}_QC_FILTERED"],
        "2": env[f"{prefix}_DOUBLET_REMOVED"],
    }
    if str(stage) not in bases:
        raise ValueError(f"stage_output_dir only supports stages 1 and 2, got {stage!r}")
    return Path(bases[str(stage)]) / leaf


def default_variant_output_dir(dataset: str, variant_id: str, env: dict[str, str]) -> Path:
    prefix = dataset_prefix(dataset)
    base = Path(env[f"{prefix}_PROCESSING_OUTPUTS"])
    return base / "03_Integrated" / variant_id


def variant_output_dir(dataset: str, variant_id: str, record: dict[str, Any], env: dict[str, str]) -> Path:
    if "output_env" in record and record["output_env"] in env:
        return Path(env[record["output_env"]])
    if "output_path" in record:
        return Path(expand_vars(str(record["output_path"]), env))
    return default_variant_output_dir(dataset, variant_id, env)


def unconfigured(value: str | os.PathLike[str]) -> bool:
    return "__UNCONFIGURED__" in str(value)


def analysis_models_for(dataset: str, requested: str, workflow: str, root: Path | None = None) -> list[str] | None:
    """Return the list of analysis models a variant gates for a workflow, or None
    if the variant does not restrict models (run the workflow's full default set)."""
    _variant_id, record = resolve_variant(dataset, requested, root)
    analysis = (record.get("analysis") or {}).get(workflow.lower())
    if isinstance(analysis, dict) and isinstance(analysis.get("models"), list):
        return analysis["models"]
    return None


def _cli(argv: list[str] | None = None) -> int:
    """Tiny query CLI so shell launchers can read the variant config without
    re-implementing YAML parsing. Examples:
        python -m rosmap_tx.config --primary tsai
        python -m rosmap_tx.config --resolve tsai strict_qc
        python -m rosmap_tx.config --analysis-models tsai hvg5000 ace
    """
    import argparse

    parser = argparse.ArgumentParser(description="Query the ROSMAP variant config.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--primary", metavar="DATASET", help="Print the dataset's primary variant id.")
    group.add_argument("--resolve", nargs=2, metavar=("DATASET", "VARIANT"),
                       help="Print the resolved variant id (primary/canonical aliases honored).")
    group.add_argument("--analysis-models", nargs=3, metavar=("DATASET", "VARIANT", "WORKFLOW"),
                       help="Print the gated analysis models (one per line), if any.")
    args = parser.parse_args(argv)

    if args.primary:
        print(primary_variant_id(args.primary))
    elif args.resolve:
        dataset, variant = args.resolve
        variant_id, _ = resolve_variant(dataset, variant)
        print(variant_id)
    elif args.analysis_models:
        dataset, variant, workflow = args.analysis_models
        models = analysis_models_for(dataset, variant, workflow)
        if models:
            print("\n".join(models))
    return 0


if __name__ == "__main__":
    raise SystemExit(_cli())

