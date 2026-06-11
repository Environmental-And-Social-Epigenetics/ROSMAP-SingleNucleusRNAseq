from __future__ import annotations

import json
import os
import shlex
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

from .config import repo_root, yaml_sha256


def _git(root: Path, *args: str) -> str | None:
    try:
        completed = subprocess.run(
            ["git", *args],
            cwd=root,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
        )
    except subprocess.CalledProcessError:
        return None
    return completed.stdout.strip()


def git_state(root: Path | None = None) -> dict[str, Any]:
    root = root or repo_root()
    status = _git(root, "status", "--short")
    return {
        "commit": _git(root, "rev-parse", "HEAD"),
        "branch": _git(root, "branch", "--show-current"),
        "dirty": bool(status),
    }


def command_string(command: list[str]) -> str:
    return " ".join(shlex.quote(part) for part in command)


def count_stage_inputs(input_dir: Path, stage: str) -> int | None:
    if not input_dir.exists():
        return None
    patterns = {
        "1": "processed_feature_bc_matrix_filtered.h5",
        "2": "*_qc.h5ad",
        "3": "*_singlets.h5ad",
    }
    pattern = patterns.get(str(stage))
    if not pattern:
        return None
    if stage == "1":
        return sum(1 for path in input_dir.glob(f"*/{pattern}") if path.is_file())
    return sum(1 for path in input_dir.glob(pattern) if path.is_file())


def write_run_manifest(
    output_dir: Path,
    *,
    dataset: str,
    stage: str,
    variant: str,
    command: list[str],
    parameters: dict[str, Any],
    inputs: dict[str, str],
    outputs: dict[str, str],
    env: dict[str, str],
    root: Path | None = None,
) -> Path:
    root = root or repo_root()
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "run_manifest.json"
    input_dir = Path(inputs["input_dir"]) if "input_dir" in inputs else None
    payload = {
        "created_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "dataset": dataset,
        "stage": str(stage),
        "variant": variant,
        "command": command_string(command),
        "argv": sys.argv,
        "git": git_state(root),
        "config_sha256": yaml_sha256(root=root),
        "environment": {
            "CONDA_PREFIX": os.environ.get("CONDA_PREFIX", ""),
            "QC_ENV": env.get("QC_ENV", ""),
            "SINGLECELL_ENV": env.get("SINGLECELL_ENV", ""),
            "BATCHCORR_ENV": env.get("BATCHCORR_ENV", ""),
        },
        "inputs": inputs,
        "outputs": outputs,
        "parameters": parameters,
        "sample_count": count_stage_inputs(input_dir, str(stage)) if input_dir else None,
    }
    manifest_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest_path

