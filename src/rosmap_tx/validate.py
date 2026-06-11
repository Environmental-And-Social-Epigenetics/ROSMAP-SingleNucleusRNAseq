from __future__ import annotations

import argparse
import fnmatch
import subprocess
import sys
from pathlib import Path

from .config import load_yaml, paths_env, repo_root, resolve_variant, unconfigured, variant_output_dir


FORBIDDEN_TRACKED_PATTERNS = [
    "*.log",
    "*.out",
    "*.err",
    "*.pyc",
    "*.h5",
    "*.h5ad",
    "*.h5seurat",
    "*.bam",
    "*.bai",
    "*.cram",
    "*.crai",
    "*.fastq",
    "*.fastq.gz",
    "*.fq",
    "*.fq.gz",
    "*.vcf",
    "*.vcf.gz",
    "*.mtx",
    "*.mtx.gz",
    "*.zarr",
    "*.pdf",
    "ckpt.tar.gz",
    "Project_*/*",
    "*/Project_*/*",
    "*/__pycache__/*",
    "Preprocessing/Tsai/02_Cellranger_Counts/Batch_Scripts/*/cellranger/*.sh",
    "Preprocessing/Tsai/02_Cellranger_Counts/Batch_Scripts/*/cellbender/*.sh",
    "Preprocessing/Tsai/03_Cellbender/Batch_Scripts/*/cellbender/*.sh",
    "Preprocessing/Tsai/01_FASTQ_Location/03_Globus_Transfer/Batch_Files/*.txt",
    "Preprocessing/Tsai/02_Cellranger_Counts/Scripts/resubmit_fixed.sh",
    "Analysis/ACE/**/tables/*",
    "Analysis/ACE/**/reports/*",
    "Analysis/ACE/**/*_report.md",
]

TRACKED_ALLOWLIST = {
    "Processing/Tsai/Pipeline/Resources/Brain_Human_PFC_Markers_Mohammadi2020.rds",
}


def _print_result(ok: bool, label: str, detail: str = "") -> int:
    status = "PASS" if ok else "FAIL"
    print(f"[{status}] {label}")
    if detail:
        print(detail)
    return 0 if ok else 1


def check_paths() -> int:
    env = paths_env()
    required = [
        "REPO_ROOT",
        "TRANSCRIPTOMICS_DATA_ROOT",
        "TSAI_PROCESSING_OUTPUTS",
        "DEJAGER_PROCESSING_OUTPUTS",
        "TSAI_INTEGRATED",
        "DEJAGER_INTEGRATED",
    ]
    missing = [name for name in required if not env.get(name) or unconfigured(env[name])]
    return _print_result(not missing, "path variables resolve", "\n".join(missing))


def check_data() -> int:
    root = repo_root()
    manifest = root / "Data" / "manifest.tsv"
    data_root = root / "Data" / "Transcriptomics"
    missing = []
    for rel in ["Tsai/FASTQs", "Tsai/Cellbender_Output", "Tsai/Processing_Outputs", "DeJager/FASTQs", "DeJager/Cellbender_Output", "DeJager/Processing_Outputs"]:
        if not (data_root / rel).exists():
            missing.append(str(data_root / rel))
    if not manifest.exists():
        missing.append(str(manifest))
    return _print_result(not missing, "data namespace present", "\n".join(missing))


def check_envs() -> int:
    root = repo_root()
    expected = [
        "envs/preprocessing/cellbender",
        "envs/preprocessing/synapse",
        "envs/preprocessing/bcftools",
        "envs/preprocessing/globus",
        "envs/processing/stage1_qc",
        "envs/processing/stage2_doublets",
        "envs/processing/stage3_integration",
        "envs/analysis/deg",
        "envs/analysis/nebula",
        "envs/analysis/sccomp",
        "envs/analysis/scenic",
        "envs/analysis/compass",
        "envs/analysis/gsea",
    ]
    missing = []
    for rel in expected:
        for name in ["environment.yml", "requirements.txt", "README.md"]:
            path = root / rel / name
            if not path.exists():
                missing.append(str(path.relative_to(root)))
    return _print_result(not missing, "canonical env specs complete", "\n".join(missing))


def _tracked_files() -> list[str]:
    completed = subprocess.run(
        ["git", "ls-files", "-z"],
        cwd=repo_root(),
        check=True,
        stdout=subprocess.PIPE,
    )
    return [item.decode() for item in completed.stdout.split(b"\0") if item]


def check_tracked_files() -> int:
    offenders: list[str] = []
    for path in _tracked_files():
        if path in TRACKED_ALLOWLIST:
            continue
        if any(fnmatch.fnmatch(path, pattern) for pattern in FORBIDDEN_TRACKED_PATTERNS):
            offenders.append(path)
    detail = "\n".join(offenders[:200])
    if len(offenders) > 200:
        detail += f"\n... {len(offenders) - 200} more"
    return _print_result(not offenders, "tracked generated-data boundary", detail)


def check_variants() -> int:
    env = paths_env()
    variants = load_yaml("variants.yaml")
    errors: list[str] = []
    for dataset, block in variants.get("datasets", {}).items():
        try:
            resolve_variant(dataset, "canonical")
        except Exception as exc:  # noqa: BLE001
            errors.append(str(exc))
        for variant_id, record in block.get("variants", {}).items():
            try:
                out = variant_output_dir(dataset, variant_id, record, env)
                if unconfigured(out):
                    errors.append(f"{dataset}/{variant_id} output is unconfigured: {out}")
                if record.get("skip_harmony") and record.get("harmony_batch_key"):
                    errors.append(f"{dataset}/{variant_id} sets both skip_harmony and harmony_batch_key")
                if not record.get("skip_harmony") and not record.get("harmony_batch_key"):
                    errors.append(f"{dataset}/{variant_id} must define harmony_batch_key or skip_harmony")
            except Exception as exc:  # noqa: BLE001
                errors.append(f"{dataset}/{variant_id}: {exc}")
    return _print_result(not errors, "variants resolve", "\n".join(errors))


CHECKS = {
    "paths": check_paths,
    "data": check_data,
    "envs": check_envs,
    "tracked-files": check_tracked_files,
    "variants": check_variants,
}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Validate ROSMAP transcriptomics repository contracts.")
    parser.add_argument("checks", nargs="+", choices=sorted(CHECKS))
    args = parser.parse_args(argv)
    status = 0
    for check in args.checks:
        status |= CHECKS[check]()
    return status


if __name__ == "__main__":
    raise SystemExit(main())
