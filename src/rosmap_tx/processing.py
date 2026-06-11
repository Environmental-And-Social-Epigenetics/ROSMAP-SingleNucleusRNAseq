from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
from typing import Any

from .config import (
    dataset_dir,
    load_yaml,
    paths_env,
    repo_root,
    resolve_variant,
    variant_output_dir,
)
from .manifest import command_string, write_run_manifest


def _add_if_present(command: list[str], flag: str, value: Any) -> None:
    if value is not None:
        command.extend([flag, str(value)])


def _stage_params(stage: str) -> dict[str, Any]:
    pipeline = load_yaml("pipeline.yaml")
    return pipeline.get(f"stage{stage}", {})


def _script(dataset: str, name: str) -> Path:
    return repo_root() / "Processing" / dataset_dir(dataset) / "Pipeline" / name


def build_command(args: argparse.Namespace, passthrough: list[str], env: dict[str, str]) -> tuple[list[str], dict[str, Any], dict[str, str], dict[str, str], Path, str]:
    dataset = args.dataset.lower()
    stage = str(args.stage)
    variant_id = "canonical"
    parameters: dict[str, Any] = {}

    if stage == "1":
        qc = _stage_params("1").get("qc_filter", {})
        input_dir = Path(env["TSAI_PREPROCESSED" if dataset == "tsai" else "DEJAGER_PREPROCESSED"])
        output_dir = Path(env["TSAI_QC_FILTERED" if dataset == "tsai" else "DEJAGER_QC_FILTERED"])
        command = [sys.executable, str(_script(dataset, "01_qc_filter.py")), "--input-dir", str(input_dir)]
        if args.list_samples:
            command.append("--list-samples")
        else:
            command.extend(["--output-dir", str(output_dir)])
        if dataset == "tsai":
            command.extend(["--metadata-csv", env["TSAI_METADATA_CSV"]])
        else:
            command.extend([
                "--patient-map-csv", env["DEJAGER_PATIENT_MAP"],
                "--patient-id-overrides-json", env["DEJAGER_PATIENT_ID_OVERRIDES"],
            ])
        for key, flag in [
            ("counts_low_pct", "--counts-low-pct"),
            ("counts_high_pct", "--counts-high-pct"),
            ("genes_low_pct", "--genes-low-pct"),
            ("mt_pct_threshold", "--mt-pct-threshold"),
        ]:
            _add_if_present(command, flag, qc.get(key))
        parameters = qc.copy()
    elif stage == "2":
        dbl = _stage_params("2").get("doublet_removal", {})
        input_dir = Path(env["TSAI_QC_FILTERED" if dataset == "tsai" else "DEJAGER_QC_FILTERED"])
        output_dir = Path(env["TSAI_DOUBLET_REMOVED" if dataset == "tsai" else "DEJAGER_DOUBLET_REMOVED"])
        command = ["Rscript", str(_script(dataset, "02_doublet_removal.Rscript")), "--input-dir", str(input_dir)]
        if args.list_samples:
            command.append("--list-samples")
        else:
            command.extend(["--output-dir", str(output_dir)])
        if dataset == "tsai":
            command.extend(["--metadata-csv", env["TSAI_METADATA_CSV"]])
        _add_if_present(command, "--threads", args.threads or dbl.get("threads"))
        parameters = dbl.copy()
    elif stage == "3":
        integration = _stage_params("3").get("integration", {})
        variant_id, record = resolve_variant(dataset, args.variant)
        input_dir = Path(env["TSAI_DOUBLET_REMOVED" if dataset == "tsai" else "DEJAGER_DOUBLET_REMOVED"])
        output_dir = variant_output_dir(dataset, variant_id, record, env)
        command = [
            sys.executable,
            str(_script(dataset, "03_integration_annotation.py")),
            "--input-dir", str(input_dir),
            "--output-dir", str(output_dir),
            "--markers-rds", env["TSAI_MARKERS_RDS" if dataset == "tsai" else "DEJAGER_MARKERS_RDS"],
            "--annotation-cluster-key", str(integration.get("annotation_cluster_key", "leiden_res0_5")),
            "--n-pcs", str(integration.get("n_pcs", 30)),
            "--n-neighbors", str(integration.get("n_neighbors", 30)),
            "--n-hvgs", str(integration.get("n_hvgs", 3000)),
            "--neighbor-metric", str(integration.get("neighbor_metric", "cosine")),
            "--umap-min-dist", str(integration.get("umap_min_dist", 0.15)),
            "--harmony-theta", str(integration.get("harmony_theta", 2.0)),
        ]
        if dataset == "tsai":
            command.extend([
                "--metadata-csv", env["TSAI_METADATA_CSV"],
                "--derived-batches-csv", env["TSAI_DERIVED_BATCHES_CSV"],
            ])
        else:
            command.extend(["--derived-batches-csv", env["DEJAGER_DERIVED_BATCHES_CSV"]])
        if record.get("skip_harmony") or args.skip_harmony:
            command.append("--skip-harmony")
        else:
            command.extend(["--harmony-batch-key", str(record["harmony_batch_key"])])
        parameters = {**integration, **record}
    else:
        raise ValueError(f"Unsupported stage: {stage}")

    if args.sample_ids:
        command.extend(["--sample-ids", args.sample_ids])
    if args.aggregate_only:
        command.append("--aggregate-only")
    if args.overwrite:
        command.append("--overwrite")
    command.extend(passthrough)

    inputs = {"input_dir": str(input_dir)}
    outputs = {"output_dir": str(output_dir)}
    return command, parameters, inputs, outputs, output_dir, variant_id


def parse_args(argv: list[str] | None = None) -> tuple[argparse.Namespace, list[str]]:
    parser = argparse.ArgumentParser(description="ROSMAP transcriptomics processing launcher.")
    parser.add_argument("--dataset", required=True, choices=["tsai", "dejager"])
    parser.add_argument("--stage", required=True, choices=["1", "2", "3"])
    parser.add_argument("--variant", default="canonical")
    parser.add_argument("--sample-ids", default="")
    parser.add_argument("--threads", type=int)
    parser.add_argument("--list-samples", action="store_true")
    parser.add_argument("--aggregate-only", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--skip-harmony", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_known_args(argv)


def main(argv: list[str] | None = None) -> int:
    args, passthrough = parse_args(argv)
    env = paths_env()
    command, parameters, inputs, outputs, output_dir, variant_id = build_command(args, passthrough, env)

    if args.dry_run:
        print(command_string(command))
        return 0

    completed = subprocess.run(command)
    if completed.returncode != 0:
        return completed.returncode

    if not args.list_samples:
        manifest = write_run_manifest(
            output_dir,
            dataset=args.dataset,
            stage=args.stage,
            variant=variant_id,
            command=command,
            parameters=parameters,
            inputs=inputs,
            outputs=outputs,
            env=env,
        )
        print(f"[manifest] wrote {manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

