#!/usr/bin/env python3
"""
Evaluate batch correction quality for Stage 3 integration outputs.

Computes metrics that measure:
  - Batch mixing (iLISI): Are batches well-mixed within cell types?
  - Cell type separation (cLISI): Are cell types still distinct after correction?
  - Silhouette score: How well-separated are cell types in the corrected embedding?
  - Cluster-batch independence: Are clusters drawn proportionally from all batches?

Usage:
    python 03b_evaluate_correction.py --input path/to/tsai_integrated.h5ad
    python 03b_evaluate_correction.py --input path/to/tsai_annotated.h5ad

    # Compare multiple correction approaches:
    python 03b_evaluate_correction.py \\
        --input run_A/tsai_integrated.h5ad run_B/tsai_integrated.h5ad \\
        --labels "derived_batch" "projid"
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Evaluate batch correction quality.")
    parser.add_argument(
        "--input",
        type=Path,
        nargs="+",
        required=True,
        help="One or more integrated h5ad files to evaluate.",
    )
    parser.add_argument(
        "--labels",
        type=str,
        nargs="*",
        default=None,
        help="Labels for each input (for comparison plots). Defaults to filenames.",
    )
    parser.add_argument(
        "--batch-key",
        type=str,
        default=None,
        help="obs column to evaluate as batch. Auto-detected from harmony_params if not set.",
    )
    parser.add_argument(
        "--cell-type-key",
        type=str,
        default="cell_type",
        help="obs column for cell type labels (default: 'cell_type').",
    )
    parser.add_argument(
        "--cluster-key",
        type=str,
        default="leiden_res0_5",
        help="obs column for cluster labels (default: 'leiden_res0_5').",
    )
    parser.add_argument(
        "--embedding-key",
        type=str,
        default=None,
        help="obsm key for corrected embedding. Auto-detected from harmony_params if not set.",
    )
    parser.add_argument(
        "--subsample",
        type=int,
        default=50000,
        help="Subsample to N cells for LISI/silhouette computation (default: 50000). "
        "Set to 0 for no subsampling.",
    )
    parser.add_argument("--output-dir", type=Path, default=None)
    return parser.parse_args()


def detect_batch_and_embedding(adata: ad.AnnData, args: argparse.Namespace) -> tuple[str, str]:
    """Auto-detect batch key and embedding key from harmony_params or fallbacks."""
    harmony_params = adata.uns.get("harmony_params", {})

    batch_key = args.batch_key
    if batch_key is None:
        batch_key = harmony_params.get("batch_key", "batch")
        if batch_key == "SKIPPED":
            batch_key = "batch"
    if batch_key not in adata.obs.columns:
        # Fallback chain
        for fallback in ["derived_batch", "batch", "projid", "sample_id"]:
            if fallback in adata.obs.columns:
                batch_key = fallback
                break

    embedding_key = args.embedding_key
    if embedding_key is None:
        embedding_key = harmony_params.get("neighbor_rep", "X_harmony")
    if embedding_key not in adata.obsm:
        embedding_key = "X_pca"

    return batch_key, embedding_key


def subsample_adata(adata: ad.AnnData, n: int, seed: int = 42) -> ad.AnnData:
    """Subsample to n cells, preserving obs/obsm."""
    if n <= 0 or adata.n_obs <= n:
        return adata
    rng = np.random.default_rng(seed)
    idx = rng.choice(adata.n_obs, size=n, replace=False)
    idx.sort()
    return adata[idx].copy()


def compute_silhouette(adata: ad.AnnData, embedding_key: str, label_key: str) -> float:
    """Compute Average Silhouette Width for cell type labels on the corrected embedding."""
    from sklearn.metrics import silhouette_score

    X = adata.obsm[embedding_key]
    labels = adata.obs[label_key].values

    # Need at least 2 unique labels
    unique = np.unique(labels)
    if len(unique) < 2:
        print(f"  [WARN] Only {len(unique)} unique label(s) for '{label_key}' — skipping silhouette")
        return float("nan")

    score = silhouette_score(X, labels, metric="euclidean", sample_size=min(10000, len(labels)))
    return score


def compute_batch_silhouette(adata: ad.AnnData, embedding_key: str, batch_key: str) -> float:
    """Compute batch ASW — should be close to 0 (well-mixed) after correction."""
    from sklearn.metrics import silhouette_score

    X = adata.obsm[embedding_key]
    labels = adata.obs[batch_key].values

    unique = np.unique(labels)
    if len(unique) < 2:
        return float("nan")

    score = silhouette_score(X, labels, metric="euclidean", sample_size=min(10000, len(labels)))
    return score


def compute_batch_entropy(adata: ad.AnnData, batch_key: str, cluster_key: str) -> float:
    """Compute mean per-cluster batch entropy — higher = better mixing."""
    from scipy.stats import entropy

    ct = pd.crosstab(adata.obs[cluster_key], adata.obs[batch_key])
    # Normalize each cluster to get batch proportions
    proportions = ct.div(ct.sum(axis=1), axis=0)
    # Compute entropy for each cluster
    cluster_entropies = proportions.apply(entropy, axis=1)
    return cluster_entropies.mean()


def evaluate_single(
    adata: ad.AnnData,
    batch_key: str,
    embedding_key: str,
    cell_type_key: str,
    cluster_key: str,
    label: str,
) -> dict[str, float]:
    """Compute all evaluation metrics for a single h5ad."""
    metrics = {"label": label, "n_cells": adata.n_obs}

    print(f"\n{'='*60}")
    print(f"Evaluating: {label}")
    print(f"  Cells: {adata.n_obs}")
    print(f"  Batch key: {batch_key} ({adata.obs[batch_key].nunique()} groups)")
    print(f"  Embedding: {embedding_key}")
    print(f"{'='*60}")

    # 1. Cell type silhouette (higher = better separation)
    if cell_type_key in adata.obs.columns:
        ct_sil = compute_silhouette(adata, embedding_key, cell_type_key)
        metrics["celltype_silhouette"] = ct_sil
        print(f"  Cell type silhouette (ASW):  {ct_sil:.4f}  (higher = better separation)")
    else:
        print(f"  [SKIP] '{cell_type_key}' not in obs — skipping cell type silhouette")

    # 2. Batch silhouette (closer to 0 = better mixing)
    batch_sil = compute_batch_silhouette(adata, embedding_key, batch_key)
    metrics["batch_silhouette"] = batch_sil
    print(f"  Batch silhouette (ASW):      {batch_sil:.4f}  (closer to 0 = better mixing)")

    # 3. Batch entropy per cluster (higher = better mixing)
    if cluster_key in adata.obs.columns:
        batch_ent = compute_batch_entropy(adata, batch_key, cluster_key)
        metrics["batch_entropy"] = batch_ent
        print(f"  Batch entropy (per cluster): {batch_ent:.4f}  (higher = better mixing)")
    else:
        print(f"  [SKIP] '{cluster_key}' not in obs — skipping batch entropy")

    # 4. Try LISI if scib is available
    try:
        from scib.metrics import ilisi_graph, clisi_graph

        ilisi = ilisi_graph(adata, batch_key=batch_key, type_="embed", use_rep=embedding_key)
        metrics["iLISI"] = float(ilisi)
        print(f"  iLISI (integration):         {ilisi:.4f}  (higher = better mixing)")

        if cell_type_key in adata.obs.columns:
            clisi = clisi_graph(adata, label_key=cell_type_key, type_="embed", use_rep=embedding_key)
            metrics["cLISI"] = float(clisi)
            print(f"  cLISI (cell type):           {clisi:.4f}  (lower = better separation)")
    except ImportError:
        print("  [SKIP] scib not installed — skipping LISI metrics")
        print("         Install with: pip install scib")
    except Exception as exc:
        print(f"  [WARN] LISI computation failed: {exc}")

    return metrics


def save_comparison(
    all_metrics: list[dict],
    output_dir: Path,
) -> None:
    """Save comparison table and plots."""
    output_dir.mkdir(parents=True, exist_ok=True)

    df = pd.DataFrame(all_metrics)
    csv_path = output_dir / "batch_correction_comparison.csv"
    df.to_csv(csv_path, index=False)
    print(f"\n[write] Comparison table -> {csv_path}")

    # Plot comparison if multiple runs
    if len(all_metrics) > 1:
        metric_cols = [c for c in df.columns if c not in ("label", "n_cells")]
        n_metrics = len(metric_cols)
        if n_metrics > 0:
            fig, axes = plt.subplots(1, n_metrics, figsize=(4 * n_metrics, 5))
            if n_metrics == 1:
                axes = [axes]
            for ax, col in zip(axes, metric_cols):
                vals = df[col].values
                ax.bar(range(len(vals)), vals, color="#4292c6")
                ax.set_xticks(range(len(vals)))
                ax.set_xticklabels(df["label"].values, rotation=45, ha="right")
                ax.set_title(col)
                for i, v in enumerate(vals):
                    if not np.isnan(v):
                        ax.text(i, v, f"{v:.3f}", ha="center", va="bottom", fontsize=8)
            fig.suptitle("Batch Correction Comparison", fontsize=13)
            fig.tight_layout()
            fig.savefig(output_dir / "batch_correction_comparison.png", dpi=150, bbox_inches="tight")
            plt.close(fig)
            print(f"[write] Comparison plot -> {output_dir / 'batch_correction_comparison.png'}")


def main() -> None:
    args = parse_args()

    labels = args.labels or [p.parent.name for p in args.input]
    if len(labels) != len(args.input):
        print("[ERROR] Number of --labels must match number of --input files")
        sys.exit(1)

    all_metrics = []
    for input_path, label in zip(args.input, labels):
        if not input_path.exists():
            print(f"[ERROR] File not found: {input_path}")
            continue

        print(f"\nLoading {input_path}...")
        adata = ad.read_h5ad(input_path)

        batch_key, embedding_key = detect_batch_and_embedding(adata, args)

        # Subsample for efficiency
        adata_sub = subsample_adata(adata, args.subsample)
        if adata_sub.n_obs < adata.n_obs:
            print(f"  Subsampled {adata.n_obs} -> {adata_sub.n_obs} cells")

        metrics = evaluate_single(
            adata_sub,
            batch_key=batch_key,
            embedding_key=embedding_key,
            cell_type_key=args.cell_type_key,
            cluster_key=args.cluster_key,
            label=label,
        )
        all_metrics.append(metrics)

    output_dir = args.output_dir or args.input[0].parent
    save_comparison(all_metrics, output_dir)

    print("\nDone.")


if __name__ == "__main__":
    main()
