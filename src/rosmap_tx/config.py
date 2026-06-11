from __future__ import annotations

import hashlib
import os
import re
import shlex
import subprocess
from pathlib import Path
from typing import Any

import yaml


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


def resolve_variant(dataset: str, requested: str, root: Path | None = None) -> tuple[str, dict[str, Any]]:
    variants = load_yaml("variants.yaml", root)
    dataset_key = dataset.lower()
    dataset_block = variants.get("datasets", {}).get(dataset_key)
    if not dataset_block:
        raise ValueError(f"No variants configured for dataset '{dataset}'")

    variant_id = requested
    if requested == "canonical":
        variant_id = dataset_block["canonical"]

    records = dataset_block.get("variants", {})
    if variant_id not in records:
        valid = ", ".join(sorted(records))
        raise ValueError(f"Unknown {dataset} variant '{requested}'. Valid variants: canonical, {valid}")
    return variant_id, records[variant_id]


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

