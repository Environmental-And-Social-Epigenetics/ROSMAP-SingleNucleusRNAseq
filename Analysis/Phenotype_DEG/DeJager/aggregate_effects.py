#!/usr/bin/env python3
"""
aggregate_effects.py -- Phase D3: aggregate per-(integration × phenotype × celltype)
DESeq2 results into a single comparison of batch correction quality based on
truth-set vs negative-control gene effect sizes.

Inputs:
  Analysis_Outputs/Phenotype_DEG/DeJager/results_<integration>/<phenotype>/deseq_<phenotype>_<celltype>.csv

Outputs (under --output-dir):
  effect_size_summary.csv       Long-format: integration × phenotype × celltype × gene_set × stat
  effect_size_ranking.csv       One row per integration with overall SNR + ranks
  effect_size_heatmap.png       Heatmap of mean truth-set |log2FC| (integration × celltype × phenotype)
  truth_vs_noise_scatter.png    Scatter of truth vs noise effect, one point per (integration, celltype, phenotype)

Interpretation:
  - High mean(|log2FC|_truth) means the batch correction preserves real biological signal
  - High mean(|log2FC|_negative_control) means the batch correction inflates spurious noise
  - SNR = truth / noise; the best batch correction maximizes this ratio
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Truth-set gene lists (relative import — same package directory)
sys.path.insert(0, str(Path(__file__).resolve().parent))
from truth_genes import PHENOTYPE_TRUTH, NEGATIVE_CONTROL  # noqa: E402

INTEGRATIONS = ["library_id", "patient_id", "pool_batch", "derived_batch", "sequencing_date"]
PHENOTYPES = ["msex", "cogdx_binary"]
CELL_TYPES = ["Exc", "Inh", "Ast", "Mic", "Oli", "OPC"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results-root", type=Path, required=True,
                   help="Path to Analysis_Outputs/Phenotype_DEG/DeJager/")
    p.add_argument("--output-dir", type=Path, required=True,
                   help="Where to write the comparison CSV + plots")
    return p.parse_args()


def load_one_result(results_root: Path, integration: str, phenotype: str, celltype: str) -> pd.DataFrame | None:
    """Load one DEG result CSV. Returns None if file is missing."""
    path = results_root / f"results_{integration}" / phenotype / f"deseq_{phenotype}_{celltype}.csv"
    if not path.exists():
        # Try broad_ prefix variant
        path_alt = results_root / f"results_{integration}" / phenotype / f"deseq_{phenotype}_broad_{celltype}.csv"
        if path_alt.exists():
            path = path_alt
        else:
            print(f"  [MISSING] {path}")
            return None
    try:
        df = pd.read_csv(path)
        return df
    except Exception as e:
        print(f"  [ERROR] {path}: {e}")
        return None


def gene_set_stats(deg_df: pd.DataFrame, gene_set: list[str]) -> dict:
    """Compute summary stats for the DEG result restricted to gene_set."""
    sub = deg_df[deg_df["gene"].isin(gene_set)].copy()
    n_present = len(sub)
    if n_present == 0:
        return {
            "n_genes_present": 0,
            "n_genes_set": len(gene_set),
            "mean_abs_lfc": np.nan,
            "median_abs_lfc": np.nan,
            "max_abs_lfc": np.nan,
            "mean_neg_log10_padj": np.nan,
            "n_sig_padj05": 0,
            "frac_sig_padj05": np.nan,
        }
    abs_lfc = sub["log2FoldChange"].abs()
    padj = sub["padj"].fillna(1.0).clip(lower=1e-300)
    return {
        "n_genes_present": int(n_present),
        "n_genes_set": len(gene_set),
        "mean_abs_lfc": float(abs_lfc.mean()),
        "median_abs_lfc": float(abs_lfc.median()),
        "max_abs_lfc": float(abs_lfc.max()),
        "mean_neg_log10_padj": float((-np.log10(padj)).mean()),
        "n_sig_padj05": int((padj < 0.05).sum()),
        "frac_sig_padj05": float((padj < 0.05).mean()),
    }


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # ----- Build long-format summary -----
    rows: list[dict] = []
    for integration in INTEGRATIONS:
        for phenotype in PHENOTYPES:
            truth_set = PHENOTYPE_TRUTH[phenotype]
            for celltype in CELL_TYPES:
                df = load_one_result(args.results_root, integration, phenotype, celltype)
                if df is None:
                    continue

                truth_stats = gene_set_stats(df, truth_set)
                noise_stats = gene_set_stats(df, NEGATIVE_CONTROL)

                base = {
                    "integration": integration,
                    "phenotype": phenotype,
                    "cell_type": celltype,
                    "n_total_genes": len(df),
                    "n_total_sig_padj05": int((df["padj"].fillna(1.0) < 0.05).sum()),
                }
                for prefix, stats in (("truth", truth_stats), ("neg_ctrl", noise_stats)):
                    for k, v in stats.items():
                        base[f"{prefix}_{k}"] = v

                # SNR = truth mean |lfc| / noise mean |lfc|
                t = truth_stats["mean_abs_lfc"]
                n = noise_stats["mean_abs_lfc"]
                base["snr_mean_abs_lfc"] = (t / n) if (n and n > 0) else np.nan

                rows.append(base)

    if not rows:
        print("[ERROR] No DEG result files found. Did Phase D2 complete?")
        sys.exit(1)

    summary = pd.DataFrame(rows)
    summary_path = args.output_dir / "effect_size_summary.csv"
    summary.to_csv(summary_path, index=False)
    print(f"[write] {summary_path} ({len(summary)} rows)")

    # ----- Per-integration ranking (collapsed across cell types) -----
    rank_rows: list[dict] = []
    for integration in INTEGRATIONS:
        sub = summary[summary["integration"] == integration]
        if sub.empty:
            continue
        for phenotype in PHENOTYPES:
            ssub = sub[sub["phenotype"] == phenotype]
            if ssub.empty:
                continue
            t_lfc = ssub["truth_mean_abs_lfc"].mean(skipna=True)
            n_lfc = ssub["neg_ctrl_mean_abs_lfc"].mean(skipna=True)
            t_sig = ssub["truth_n_sig_padj05"].sum(skipna=True)
            rank_rows.append({
                "integration": integration,
                "phenotype": phenotype,
                "n_celltypes": len(ssub),
                "truth_mean_abs_lfc": float(t_lfc),
                "neg_ctrl_mean_abs_lfc": float(n_lfc),
                "snr_mean_abs_lfc": float(t_lfc / n_lfc) if n_lfc and n_lfc > 0 else float("nan"),
                "truth_total_sig_padj05": int(t_sig),
            })

    rank = pd.DataFrame(rank_rows)
    if len(rank) > 0:
        # Add ordinal ranks per phenotype (1 = best SNR)
        rank["rank_within_phenotype"] = (
            rank.groupby("phenotype")["snr_mean_abs_lfc"]
            .rank(ascending=False, method="dense")
            .astype("Int64")
        )
        # Combined score = average SNR across phenotypes
        combined = rank.groupby("integration")["snr_mean_abs_lfc"].mean().rename("combined_snr")
        rank = rank.merge(combined, on="integration")
        rank = rank.sort_values(["phenotype", "rank_within_phenotype"])
        rank_path = args.output_dir / "effect_size_ranking.csv"
        rank.to_csv(rank_path, index=False)
        print(f"[write] {rank_path}")

    # ----- Heatmap: integration × cell_type, faceted by phenotype, colored by truth_mean_abs_lfc -----
    fig, axes = plt.subplots(1, len(PHENOTYPES), figsize=(5 * len(PHENOTYPES), 5), squeeze=False)
    for j, phenotype in enumerate(PHENOTYPES):
        ax = axes[0][j]
        sub = summary[summary["phenotype"] == phenotype]
        if sub.empty:
            ax.set_visible(False)
            continue
        pivot = sub.pivot_table(index="integration", columns="cell_type",
                                values="truth_mean_abs_lfc", aggfunc="mean")
        pivot = pivot.reindex(index=[i for i in INTEGRATIONS if i in pivot.index],
                              columns=[c for c in CELL_TYPES if c in pivot.columns])
        im = ax.imshow(pivot.values, aspect="auto", cmap="viridis")
        ax.set_xticks(range(len(pivot.columns)))
        ax.set_xticklabels(pivot.columns, rotation=45, ha="right")
        ax.set_yticks(range(len(pivot.index)))
        ax.set_yticklabels(pivot.index)
        ax.set_title(f"{phenotype}\nmean |log2FC| on truth set")
        for r in range(pivot.shape[0]):
            for c in range(pivot.shape[1]):
                v = pivot.values[r, c]
                if np.isfinite(v):
                    ax.text(c, r, f"{v:.2f}", ha="center", va="center",
                            color="white" if v < np.nanmean(pivot.values) else "black", fontsize=8)
        fig.colorbar(im, ax=ax, shrink=0.7)
    fig.suptitle("Truth-set effect size by batch correction × cell type", fontsize=12)
    fig.tight_layout()
    heatmap_path = args.output_dir / "effect_size_heatmap.png"
    fig.savefig(heatmap_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[write] {heatmap_path}")

    # ----- Scatter: truth vs noise effect, one point per (integration, celltype, phenotype) -----
    fig, axes = plt.subplots(1, len(PHENOTYPES), figsize=(6 * len(PHENOTYPES), 5), squeeze=False)
    cmap = plt.colormaps.get_cmap("tab10")
    integ_color = {integ: cmap(i) for i, integ in enumerate(INTEGRATIONS)}
    for j, phenotype in enumerate(PHENOTYPES):
        ax = axes[0][j]
        sub = summary[summary["phenotype"] == phenotype]
        if sub.empty:
            ax.set_visible(False)
            continue
        for integration in INTEGRATIONS:
            ssub = sub[sub["integration"] == integration]
            if ssub.empty:
                continue
            ax.scatter(ssub["neg_ctrl_mean_abs_lfc"], ssub["truth_mean_abs_lfc"],
                       label=integration, color=integ_color[integration],
                       s=70, alpha=0.8, edgecolors="black")
            for _, r in ssub.iterrows():
                ax.annotate(r["cell_type"],
                            (r["neg_ctrl_mean_abs_lfc"], r["truth_mean_abs_lfc"]),
                            fontsize=7, alpha=0.7)
        # Diagonal y=x line for reference
        all_vals = np.concatenate([sub["neg_ctrl_mean_abs_lfc"].dropna().values,
                                   sub["truth_mean_abs_lfc"].dropna().values])
        if len(all_vals):
            lim = float(np.nanmax(all_vals)) * 1.1
            ax.plot([0, lim], [0, lim], "k--", alpha=0.3, label="y=x (no signal)")
            ax.set_xlim(0, lim)
            ax.set_ylim(0, lim)
        ax.set_xlabel("Mean |log2FC| on housekeeping genes (noise)")
        ax.set_ylabel("Mean |log2FC| on truth set (signal)")
        ax.set_title(f"{phenotype}: above the line = real biology")
        ax.legend(fontsize=8, loc="best")
    fig.suptitle("Signal vs. noise per cell type for each batch correction", fontsize=12)
    fig.tight_layout()
    scatter_path = args.output_dir / "truth_vs_noise_scatter.png"
    fig.savefig(scatter_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"[write] {scatter_path}")

    print("\nDone.")


if __name__ == "__main__":
    main()
