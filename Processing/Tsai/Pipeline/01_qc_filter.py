#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import time
from pathlib import Path

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
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
        "--counts-low-pct",
        type=float,
        default=4.5,
        help="Lower percentile threshold for log1p_total_counts outlier detection.",
    )
    parser.add_argument(
        "--counts-high-pct",
        type=float,
        default=96.0,
        help="Upper percentile threshold for log1p_total_counts outlier detection.",
    )
    parser.add_argument(
        "--genes-low-pct",
        type=float,
        default=5.0,
        help="Lower percentile threshold for log1p_n_genes_by_counts outlier detection.",
    )
    parser.add_argument(
        "--mt-pct-threshold",
        type=float,
        default=10.0,
        help="Hard upper threshold for mitochondrial percentage.",
    )
    parser.add_argument(
        "--aggregate-only",
        action="store_true",
        help="Skip per-sample processing; generate aggregate plots from qc_summary.csv.",
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
    adata: sc.AnnData, metric: str, lower_pct: float | None, upper_pct: float | None
) -> pd.Series:
    values = adata.obs[metric]
    mask = pd.Series(False, index=adata.obs.index)
    if lower_pct is not None:
        mask = mask | (values < np.percentile(values, lower_pct))
    if upper_pct is not None:
        mask = mask | (values > np.percentile(values, upper_pct))
    return mask


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
    """Append a single row to the shared QC summary CSV (mkdir-lock safe)."""
    header = "sample_id,n_cells_before,n_cells_after,n_removed,pct_retained,n_outlier,n_mt_outlier,median_genes,median_counts,median_pct_mt\n"
    row = (
        f"{sample_id},{n_cells_before},{n_cells_after},{n_cells_before - n_cells_after},"
        f"{100 * n_cells_after / max(n_cells_before, 1):.1f},"
        f"{n_outlier},{n_mt_outlier},{median_genes:.1f},{median_counts:.1f},{median_pct_mt:.2f}\n"
    )
    lock_dir = str(summary_path) + ".lock"
    acquired = False
    for _ in range(300):
        try:
            os.mkdir(lock_dir)
            acquired = True
            break
        except FileExistsError:
            time.sleep(0.1)
    if not acquired:
        print(f"[WARN] Lock timeout for {summary_path}, proceeding without lock")
    try:
        write_header = not summary_path.exists() or summary_path.stat().st_size == 0
        with open(summary_path, "a") as fh:
            if write_header:
                fh.write(header)
            fh.write(row)
    finally:
        if acquired:
            try:
                os.rmdir(lock_dir)
            except OSError:
                pass


def plot_qc_persample(
    adata: sc.AnnData,
    sample_id: str,
    figures_dir: Path,
    args: argparse.Namespace,
) -> None:
    """Generate per-sample QC visualizations (3 PNGs)."""
    figures_dir.mkdir(parents=True, exist_ok=True)

    # --- 1. QC metric violins with percentile threshold lines ---
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    metric_specs = [
        ("total_counts", "Total Counts"),
        ("n_genes_by_counts", "Genes Detected"),
        ("pct_counts_mt", "MT %"),
    ]
    for ax, (metric, label) in zip(axes, metric_specs):
        values = adata.obs[metric].values
        ax.violinplot(values, showmedians=True)
        if metric == "pct_counts_mt":
            ax.axhline(args.mt_pct_threshold, color="red", linestyle="--", alpha=0.7,
                        label=f"threshold ({args.mt_pct_threshold}%)")
        elif metric == "total_counts":
            log_metric = "log1p_total_counts"
            if log_metric in adata.obs.columns:
                log_vals = adata.obs[log_metric].values
                lo_log = np.percentile(log_vals, args.counts_low_pct)
                hi_log = np.percentile(log_vals, args.counts_high_pct)
                lo_orig = np.expm1(lo_log)
                hi_orig = np.expm1(hi_log)
                ax.axhline(lo_orig, color="red", linestyle="--", alpha=0.7, label=f"P{args.counts_low_pct} ({lo_orig:.0f})")
                ax.axhline(hi_orig, color="red", linestyle="--", alpha=0.7, label=f"P{args.counts_high_pct} ({hi_orig:.0f})")
        elif metric == "n_genes_by_counts":
            log_metric = "log1p_n_genes_by_counts"
            if log_metric in adata.obs.columns:
                log_vals = adata.obs[log_metric].values
                lo_log = np.percentile(log_vals, args.genes_low_pct)
                lo_orig = np.expm1(lo_log)
                ax.axhline(lo_orig, color="red", linestyle="--", alpha=0.7, label=f"P{args.genes_low_pct} ({lo_orig:.0f})")
        ax.set_title(label)
        ax.set_xticks([])
        ax.legend(fontsize=7, loc="upper right")
    fig.suptitle(f"{sample_id} — QC Metrics (n={adata.n_obs})", fontsize=12)
    fig.tight_layout()
    fig.savefig(figures_dir / f"{sample_id}_qc_violins.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- 2. Genes vs counts scatter, colored by MT% ---
    fig, ax = plt.subplots(figsize=(7, 6))
    sc_plot = ax.scatter(
        adata.obs["total_counts"],
        adata.obs["n_genes_by_counts"],
        c=adata.obs["pct_counts_mt"],
        cmap="viridis",
        s=3,
        alpha=0.5,
        rasterized=True,
    )
    fig.colorbar(sc_plot, ax=ax, label="MT %")
    ax.set_xlabel("Total Counts")
    ax.set_ylabel("Genes Detected")
    ax.set_title(f"{sample_id} — Genes vs Counts")
    fig.savefig(figures_dir / f"{sample_id}_genes_vs_counts.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- 3. Filter overlay — kept vs removed cells ---
    kept = (~adata.obs["outlier"]) & (~adata.obs["mt_outlier"])
    outlier_only = adata.obs["outlier"] & (~adata.obs["mt_outlier"])
    mt_only = adata.obs["mt_outlier"] & (~adata.obs["outlier"])
    both = adata.obs["outlier"] & adata.obs["mt_outlier"]

    fig, ax = plt.subplots(figsize=(7, 6))
    groups = [
        (kept, "Kept", "#999999", 0.3),
        (outlier_only, "Outlier", "#e41a1c", 0.6),
        (mt_only, "MT outlier", "#ff7f00", 0.6),
        (both, "Both", "#800000", 0.7),
    ]
    for mask, label, color, alpha in groups:
        if mask.sum() == 0:
            continue
        ax.scatter(
            adata.obs.loc[mask, "total_counts"],
            adata.obs.loc[mask, "n_genes_by_counts"],
            c=color, s=3, alpha=alpha, label=f"{label} ({mask.sum()})",
            rasterized=True,
        )
    ax.set_xlabel("Total Counts")
    ax.set_ylabel("Genes Detected")
    ax.set_title(f"{sample_id} — Filtering Overlay")
    ax.legend(fontsize=8, markerscale=3)
    fig.savefig(figures_dir / f"{sample_id}_filter_overlay.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_qc_aggregate(args: argparse.Namespace) -> None:
    """Generate aggregate QC plots from qc_summary.csv."""
    summary_path = args.output_dir / "qc_summary.csv"
    if not summary_path.exists():
        print(f"[WARN] {summary_path} not found, skipping aggregate plots")
        return

    df = pd.read_csv(summary_path, dtype={"sample_id": str})
    figures_dir = args.output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)

    # --- 1. Cells retained bar plot ---
    df_sorted = df.sort_values("pct_retained")
    fig, ax = plt.subplots(figsize=(8, max(12, len(df_sorted) * 0.06)))
    colors = plt.cm.RdYlGn(df_sorted["pct_retained"].values / 100)
    ax.barh(range(len(df_sorted)), df_sorted["pct_retained"], color=colors, edgecolor="none")
    ax.axvline(df_sorted["pct_retained"].median(), color="black", linestyle="--", alpha=0.7,
               label=f"Median ({df_sorted['pct_retained'].median():.1f}%)")
    ax.set_yticks(range(len(df_sorted)))
    ax.set_yticklabels(df_sorted["sample_id"], fontsize=3)
    ax.set_xlabel("% Cells Retained")
    ax.set_title("QC: Cells Retained per Sample")
    ax.legend(fontsize=9)
    fig.savefig(figures_dir / "aggregate_cells_retained.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- 2. QC metric distributions across samples ---
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    for ax, (col, label) in zip(axes, [
        ("median_genes", "Median Genes per Cell"),
        ("median_counts", "Median Counts per Cell"),
        ("median_pct_mt", "Median MT %"),
    ]):
        bp = ax.boxplot(df[col].dropna(), vert=True, patch_artist=True)
        bp["boxes"][0].set_facecolor("#4292c6")
        ax.set_title(label)
        ax.set_xticks([])
        median_val = df[col].median()
        ax.axhline(median_val, color="red", linestyle="--", alpha=0.5)
        ax.text(1.15, median_val, f"{median_val:.1f}", va="center", fontsize=8, color="red")
    fig.suptitle(f"QC Metric Distributions (n={len(df)} samples)", fontsize=12)
    fig.tight_layout()
    fig.savefig(figures_dir / "aggregate_qc_distributions.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- 3. Attrition breakdown ---
    df_sorted = df.sort_values("n_removed", ascending=True)
    fig, ax = plt.subplots(figsize=(8, max(12, len(df_sorted) * 0.06)))
    ax.barh(range(len(df_sorted)), df_sorted["n_outlier"], color="#e41a1c", alpha=0.8, label="Outlier")
    ax.barh(range(len(df_sorted)), df_sorted["n_mt_outlier"], left=df_sorted["n_outlier"],
            color="#ff7f00", alpha=0.8, label="MT outlier")
    ax.set_yticks(range(len(df_sorted)))
    ax.set_yticklabels(df_sorted["sample_id"], fontsize=3)
    ax.set_xlabel("Cells Removed")
    ax.set_title("QC: Attrition Breakdown per Sample")
    ax.legend(fontsize=9)
    fig.savefig(figures_dir / "aggregate_attrition_breakdown.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"[done] Aggregate QC plots saved to {figures_dir}")


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
        is_outlier_percentile(adata, "log1p_total_counts", args.counts_low_pct, args.counts_high_pct)
        | is_outlier_percentile(adata, "log1p_n_genes_by_counts", args.genes_low_pct, None)
    )
    adata.obs["mt_outlier"] = adata.obs["pct_counts_mt"] > args.mt_pct_threshold

    figures_dir = args.output_dir / "figures"
    plot_qc_persample(adata, sample_id, figures_dir, args)

    pre_filter_n_cells = adata.n_obs
    filtered = adata[
        (~adata.obs["outlier"])
        & (~adata.obs["mt_outlier"])
    ].copy()
    filtered.uns["qc_filtering"] = {
        "method": "percentile",
        "counts_low_pct": args.counts_low_pct,
        "counts_high_pct": args.counts_high_pct,
        "genes_low_pct": args.genes_low_pct,
        "mt_pct_threshold": args.mt_pct_threshold,
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

    if args.aggregate_only:
        plot_qc_aggregate(args)
        return

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
