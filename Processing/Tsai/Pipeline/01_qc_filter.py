#!/usr/bin/env python3
from __future__ import annotations

import argparse
import fcntl
import os
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import median_abs_deviation


os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
sc.settings.verbosity = 0


def default_paths() -> dict[str, Path]:
    script_path = Path(__file__).resolve()
    repo_root = script_path.parents[3]
    workspace_root = repo_root.parent
    return {
        "repo_root": repo_root,
        "workspace_root": workspace_root,
        "input_dir": workspace_root / "Tsai_Data" / "Cellbender_Outputs",
        "output_dir": workspace_root / "Tsai_Data" / "Processing_Outputs" / "01_QC_Filtered",
        "metadata_csv": repo_root
        / "Preprocessing"
        / "Tsai"
        / "02_Cellranger_Counts"
        / "Tracking"
        / "patient_metadata.csv",
    }


def parse_args() -> argparse.Namespace:
    paths = default_paths()
    parser = argparse.ArgumentParser(
        description="Stage 1 Tsai snRNA-seq QC filtering from CellBender outputs."
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=paths["input_dir"],
        help="Directory containing per-sample CellBender output folders.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=paths["output_dir"],
        help="Directory to write per-sample QC-filtered h5ad files.",
    )
    parser.add_argument(
        "--metadata-csv",
        type=Path,
        default=paths["metadata_csv"],
        help="Metadata CSV used to preserve sample ordering.",
    )
    parser.add_argument(
        "--sample-ids",
        type=str,
        default="",
        help="Comma-separated list of sample IDs to process. Defaults to all complete samples.",
    )
    parser.add_argument(
        "--list-samples",
        action="store_true",
        help="Print discoverable complete sample IDs and exit.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing outputs instead of skipping them.",
    )
    parser.add_argument(
        "--counts-lower-percentile",
        type=float,
        default=4.5,
        help="Lower percentile threshold for log1p total counts.",
    )
    parser.add_argument(
        "--counts-upper-percentile",
        type=float,
        default=96.0,
        help="Upper percentile threshold for log1p total counts.",
    )
    parser.add_argument(
        "--genes-lower-percentile",
        type=float,
        default=5.0,
        help="Lower percentile threshold for log1p genes by counts.",
    )
    parser.add_argument(
        "--genes-upper-percentile",
        type=float,
        default=100.0,
        help="Upper percentile threshold for log1p genes by counts.",
    )
    parser.add_argument(
        "--max-mt",
        type=float,
        default=10.0,
        help="Maximum allowed mitochondrial percentage.",
    )
    parser.add_argument(
        "--top20-nmads",
        type=float,
        default=4.0,
        help="Number of MADs for pct_counts_in_top_20_genes outlier detection. Set to 0 to disable.",
    )
    return parser.parse_args()


def discover_complete_samples(input_dir: Path, metadata_csv: Path | None) -> list[str]:
    complete_samples = sorted(
        path.name
        for path in input_dir.iterdir()
        if path.is_dir() and (path / "processed_feature_bc_matrix_filtered.h5").exists()
    )

    if metadata_csv is None or not metadata_csv.exists():
        return complete_samples

    metadata = pd.read_csv(metadata_csv, dtype={"projid": str})
    metadata_samples = metadata["projid"].dropna().astype(str).tolist()
    ordered_samples = [sample for sample in metadata_samples if sample in set(complete_samples)]
    remaining_samples = sorted(set(complete_samples) - set(ordered_samples))
    return ordered_samples + remaining_samples


def parse_requested_samples(raw_sample_ids: str, available_samples: list[str]) -> list[str]:
    if not raw_sample_ids.strip():
        return available_samples

    requested = [sample.strip() for sample in raw_sample_ids.split(",") if sample.strip()]
    missing = [sample for sample in requested if sample not in set(available_samples)]
    if missing:
        missing_str = ", ".join(missing)
        raise ValueError(f"Requested sample IDs are not available as complete CellBender outputs: {missing_str}")
    return requested


def is_outlier_percentile(
    adata: sc.AnnData, metric: str, lower_percentile: float, upper_percentile: float
) -> pd.Series:
    values = adata.obs[metric]
    lower_threshold = np.percentile(values, lower_percentile)
    upper_threshold = np.percentile(values, upper_percentile)
    return (values < lower_threshold) | (values > upper_threshold)


def is_outlier_mad(adata: sc.AnnData, metric: str, nmads: int) -> pd.Series:
    values = adata.obs[metric]
    deviation = median_abs_deviation(values)
    median_value = np.median(values)
    return (values < median_value - nmads * deviation) | (values > median_value + nmads * deviation)


def append_qc_summary(
    summary_path: Path,
    sample_id: str,
    n_cells_before: int,
    n_cells_after: int,
    n_outlier: int,
    n_mt_outlier: int,
    median_genes: float,
    median_counts: float,
    median_pct_mt: float,
) -> None:
    """Append a single row to the shared QC summary CSV (file-lock safe)."""
    header = "sample_id,n_cells_before,n_cells_after,n_removed,pct_retained,n_outlier,n_mt_outlier,median_genes,median_counts,median_pct_mt\n"
    row = (
        f"{sample_id},{n_cells_before},{n_cells_after},{n_cells_before - n_cells_after},"
        f"{100 * n_cells_after / max(n_cells_before, 1):.1f},"
        f"{n_outlier},{n_mt_outlier},{median_genes:.1f},{median_counts:.1f},{median_pct_mt:.2f}\n"
    )
    with open(summary_path, "a") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        if fh.tell() == 0:
            fh.write(header)
        fh.write(row)
        fcntl.flock(fh, fcntl.LOCK_UN)


def run_qc_filter(sample_id: str, args: argparse.Namespace) -> Path:
    input_path = args.input_dir / sample_id / "processed_feature_bc_matrix_filtered.h5"
    output_path = args.output_dir / f"{sample_id}_qc.h5ad"

    if output_path.exists() and not args.overwrite:
        print(f"[skip] {sample_id}: {output_path} already exists")
        return output_path

    adata = sc.read_10x_h5(input_path)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    adata.obs["sample_id"] = sample_id

    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.contains(r"^HB[^P]", regex=True)

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo", "hb"],
        inplace=True,
        percent_top=[20],
        log1p=True,
    )

    adata.obs["outlier"] = (
        is_outlier_percentile(
            adata,
            "log1p_total_counts",
            args.counts_lower_percentile,
            args.counts_upper_percentile,
        )
        | is_outlier_percentile(
            adata,
            "log1p_n_genes_by_counts",
            args.genes_lower_percentile,
            args.genes_upper_percentile,
        )
    )
    adata.obs["mt_outlier"] = adata.obs["pct_counts_mt"] > args.max_mt

    if args.top20_nmads > 0:
        adata.obs["top20_outlier"] = is_outlier_mad(
            adata, "pct_counts_in_top_20_genes", int(args.top20_nmads)
        )
    else:
        adata.obs["top20_outlier"] = False

    adata.obs["mad_total_counts_outlier"] = is_outlier_mad(adata, "log1p_total_counts", 5)
    adata.obs["mad_n_genes_outlier"] = is_outlier_mad(adata, "log1p_n_genes_by_counts", 5)

    pre_filter_n_cells = adata.n_obs
    filtered = adata[
        (~adata.obs["outlier"])
        & (~adata.obs["mt_outlier"])
        & (~adata.obs["top20_outlier"])
    ].copy()
    filtered.uns["qc_filtering"] = {
        "counts_lower_percentile": args.counts_lower_percentile,
        "counts_upper_percentile": args.counts_upper_percentile,
        "genes_lower_percentile": args.genes_lower_percentile,
        "genes_upper_percentile": args.genes_upper_percentile,
        "max_mt": args.max_mt,
        "top20_nmads": args.top20_nmads,
        "input_path": str(input_path),
    }

    filtered.write_h5ad(output_path)

    summary_path = args.output_dir / "qc_summary.csv"
    append_qc_summary(
        summary_path,
        sample_id=sample_id,
        n_cells_before=pre_filter_n_cells,
        n_cells_after=filtered.n_obs,
        n_outlier=int(adata.obs["outlier"].sum()),
        n_mt_outlier=int(adata.obs["mt_outlier"].sum()),
        median_genes=float(np.median(adata.obs["n_genes_by_counts"])),
        median_counts=float(np.median(adata.obs["total_counts"])),
        median_pct_mt=float(np.median(adata.obs["pct_counts_mt"])),
    )

    print(f"[done] {sample_id}: kept {filtered.n_obs}/{pre_filter_n_cells} cells -> {output_path}")
    return output_path


def main() -> None:
    args = parse_args()
    args.input_dir = args.input_dir.resolve()
    args.output_dir = args.output_dir.resolve()
    args.metadata_csv = args.metadata_csv.resolve()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    available_samples = discover_complete_samples(args.input_dir, args.metadata_csv)
    if args.list_samples:
        for sample_id in available_samples:
            print(sample_id)
        return

    selected_samples = parse_requested_samples(args.sample_ids, available_samples)
    if not selected_samples:
        raise SystemExit("No complete CellBender outputs were found to process.")

    for sample_id in selected_samples:
        run_qc_filter(sample_id, args)


if __name__ == "__main__":
    main()
