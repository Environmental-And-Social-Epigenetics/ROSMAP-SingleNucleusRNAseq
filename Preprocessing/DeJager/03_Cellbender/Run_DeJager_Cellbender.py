#!/usr/bin/env python3
"""
Submit CellBender jobs for DeJager libraries from tracked config.

Usage:
    source config/paths.sh
    python Preprocessing/DeJager/03_Cellbender/Run_DeJager_Cellbender.py --dry-run
    python Preprocessing/DeJager/03_Cellbender/Run_DeJager_Cellbender.py --submit
"""

from __future__ import annotations

import argparse
import os
import shlex
import subprocess
from pathlib import Path


def env_path(name: str, default: str) -> Path:
    return Path(os.environ.get(name, default)).expanduser()


def default_paths() -> dict[str, Path]:
    script_path = Path(__file__).resolve()
    repo_root = script_path.parents[3]
    workspace_root = repo_root.parent
    input_dir = env_path("DEJAGER_COUNTS", str(workspace_root / "DeJager_Data" / "Counts"))
    output_dir = env_path(
        "DEJAGER_PREPROCESSED",
        str(workspace_root / "DeJager_Data" / "Cellbender_Output"),
    )
    return {
        "repo_root": repo_root,
        "input_dir": input_dir,
        "output_dir": output_dir,
        "log_dir": output_dir.parent / "Cellbender_Logs",
    }


def parse_args() -> argparse.Namespace:
    paths = default_paths()
    parser = argparse.ArgumentParser(description="Submit DeJager CellBender jobs via sbatch --wrap.")
    parser.add_argument("--input-dir", type=Path, default=paths["input_dir"])
    parser.add_argument("--output-dir", type=Path, default=paths["output_dir"])
    parser.add_argument("--log-dir", type=Path, default=paths["log_dir"])
    parser.add_argument("--library-ids", type=str, default="", help="Comma-separated library IDs to process.")
    parser.add_argument("--limit", type=int, default=0, help="Process only the first N discovered libraries.")
    parser.add_argument("--fpr", type=float, default=0.0)
    parser.add_argument("--time", type=str, default="47:00:00")
    parser.add_argument("--cpus", type=int, default=32)
    parser.add_argument("--mem", type=str, default="128G")
    parser.add_argument("--gpu", type=str, default="a100:1")
    parser.add_argument("--submit", action="store_true", help="Actually submit jobs to SLURM.")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def discover_libraries(input_dir: Path) -> list[str]:
    libraries: list[str] = []
    for library_dir in sorted(input_dir.iterdir()):
        raw_h5 = library_dir / "outs" / "raw_feature_bc_matrix.h5"
        if library_dir.is_dir() and raw_h5.exists():
            libraries.append(library_dir.name)
    return libraries


def selected_libraries(args: argparse.Namespace) -> list[str]:
    libraries = discover_libraries(args.input_dir)
    if args.library_ids.strip():
        requested = [lib.strip() for lib in args.library_ids.split(",") if lib.strip()]
        missing = sorted(set(requested) - set(libraries))
        if missing:
            raise SystemExit(f"Requested libraries not found in {args.input_dir}: {', '.join(missing)}")
        libraries = requested
    if args.limit > 0:
        libraries = libraries[: args.limit]
    return libraries


def build_sbatch_command(args: argparse.Namespace, library_id: str) -> list[str]:
    conda_init = os.environ.get("CONDA_INIT_SCRIPT", "")
    cellbender_env = os.environ.get("CELLBENDER_ENV", "")
    mail_user = os.environ.get("SLURM_MAIL_USER", "")
    partition = os.environ.get("SLURM_PARTITION_GPU", "") or os.environ.get("SLURM_PARTITION", "")

    input_h5 = args.input_dir / library_id / "outs" / "raw_feature_bc_matrix.h5"
    output_dir = args.output_dir / library_id
    output_h5 = output_dir / "processed_feature_bc_matrix.h5"
    log_out = args.log_dir / f"{library_id}_%j.out"
    log_err = args.log_dir / f"{library_id}_%j.err"

    wrap = " && ".join(
        [
            f"source {shlex.quote(conda_init)}",
            f"conda activate {shlex.quote(cellbender_env)}",
            f"mkdir -p {shlex.quote(str(output_dir))}",
            (
                "cellbender remove-background --cuda "
                f"--input {shlex.quote(str(input_h5))} "
                f"--fpr {args.fpr} "
                f"--output {shlex.quote(str(output_h5))}"
            ),
        ]
    )

    cmd = [
        "sbatch",
        "--job-name",
        f"dej_cb_{library_id}",
        "--time",
        args.time,
        "--ntasks",
        str(args.cpus),
        "--mem",
        args.mem,
        "--gres",
        f"gpu:{args.gpu}",
        "--mail-type",
        "FAIL",
        "--output",
        str(log_out),
        "--error",
        str(log_err),
    ]
    if partition:
        cmd.extend(["--partition", partition])
    if mail_user:
        cmd.extend(["--mail-user", mail_user])
    cmd.extend(["--wrap", wrap])
    return cmd


def main() -> None:
    args = parse_args()
    args.input_dir = args.input_dir.resolve()
    args.output_dir = args.output_dir.resolve()
    args.log_dir = args.log_dir.resolve()

    libraries = selected_libraries(args)
    if not libraries:
        raise SystemExit(f"No Cell Ranger raw matrices found under {args.input_dir}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.log_dir.mkdir(parents=True, exist_ok=True)

    print(f"Found {len(libraries)} libraries to process")
    for library_id in libraries:
        output_h5 = args.output_dir / library_id / "processed_feature_bc_matrix.h5"
        if output_h5.exists() and not args.overwrite:
            print(f"[skip] {library_id}: {output_h5} already exists")
            continue

        cmd = build_sbatch_command(args, library_id)
        print("[plan]", " ".join(shlex.quote(part) for part in cmd))
        if args.submit and not args.dry_run:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            print(f"[submit] {library_id}: {result.stdout.strip()}")


if __name__ == "__main__":
    main()
