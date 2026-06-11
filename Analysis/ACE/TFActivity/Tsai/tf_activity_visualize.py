#!/usr/bin/env python3
"""
ACE TF & Pathway Activity Visualization

Generates publication-quality figures from the tf_summary.csv and
pathway_summary.csv produced by tf_activity_analysis.py.

Outputs (all PDF):
  1. Cross-cell-type TF heatmap (significant TFs only)
  2. Per-cell-type volcano plots
  3. Epigenetic TF highlight panel
  4. Hub TF barplot (TFs significant in 3+ cell types)
  5. Pathway heatmap (if pathway_summary.csv exists)

Usage:
    python tf_activity_visualize.py \
        --results-dir ${ACE_OUTPUT_ROOT}/TFActivity/Tsai/results_derived_batch/tot_adverse_exp \
        --output-dir  ${ACE_OUTPUT_ROOT}/TFActivity/Tsai/figures_derived_batch/tot_adverse_exp \
        --phenotype tot_adverse_exp \
        --sex Male
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path

sys.stdout.reconfigure(line_buffering=True)
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    stream=sys.stdout,
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Known epigenetic regulators for the highlight panel
# ---------------------------------------------------------------------------

EPIGENETIC_TFS = {
    "HDAC1", "HDAC2", "HDAC3", "HDAC4", "HDAC5", "HDAC6", "HDAC7",
    "HDAC8", "HDAC9", "HDAC10", "HDAC11",
    "DNMT1", "DNMT3A", "DNMT3B",
    "EZH2", "SUZ12", "REST", "SIN3A", "CTCF",
    "KDM5A", "KDM5B",
    "KAT2A", "KAT2B",
    "EP300", "CREBBP",
}


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------


def cross_celltype_heatmap(
    tf_df: pd.DataFrame, output_dir: Path, phenotype: str, sex: str
) -> None:
    """Heatmap of TF coefficients: rows = TFs significant in >= 1 cell type."""
    sig = tf_df[tf_df["padj"] < 0.05]
    if sig.empty:
        log.warning("No significant TFs for heatmap (sex=%s)", sex)
        return

    sig_names = sig["name"].unique()
    subset = tf_df[tf_df["name"].isin(sig_names)]

    pivot_coef = subset.pivot_table(
        index="name", columns="cell_type", values="coef", aggfunc="first"
    )
    pivot_padj = subset.pivot_table(
        index="name", columns="cell_type", values="padj", aggfunc="first"
    )

    # Annotation: * for padj < 0.05
    annot = pivot_padj.map(lambda x: "*" if pd.notna(x) and x < 0.05 else "")

    # Limit to top 50 TFs by minimum p-value across cell types for readability
    if len(pivot_coef) > 50:
        min_padj = pivot_padj.min(axis=1).sort_values()
        top_names = min_padj.head(50).index
        pivot_coef = pivot_coef.loc[top_names]
        annot = annot.loc[top_names]

    fig_height = max(8, len(pivot_coef) * 0.3)
    fig_width = max(6, len(pivot_coef.columns) * 1.2 + 2)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    sns.heatmap(
        pivot_coef,
        cmap="RdBu_r",
        center=0,
        annot=annot,
        fmt="",
        linewidths=0.5,
        ax=ax,
        cbar_kws={"label": f"Coefficient ({phenotype})"},
    )
    ax.set_title(f"TF Activity vs {phenotype} ({sex})\n* = FDR < 0.05")
    ax.set_ylabel("Transcription Factor")
    ax.set_xlabel("Cell Type")
    plt.tight_layout()

    out = output_dir / f"tf_heatmap_{sex}.pdf"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    log.info("Saved heatmap: %s", out)


def per_celltype_volcano(
    tf_df: pd.DataFrame, output_dir: Path, phenotype: str, sex: str
) -> None:
    """Volcano plot per cell type: coefficient vs -log10(padj)."""
    cell_types = tf_df["cell_type"].unique()
    n = len(cell_types)
    if n == 0:
        return

    ncols = min(3, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 5 * nrows), squeeze=False)

    for idx, ct in enumerate(sorted(cell_types)):
        ax = axes[idx // ncols, idx % ncols]
        ct_data = tf_df[tf_df["cell_type"] == ct].copy()
        ct_data["neg_log10_padj"] = -np.log10(ct_data["padj"].clip(lower=1e-300))

        colors = []
        for _, row in ct_data.iterrows():
            if row["padj"] < 0.05 and row["coef"] > 0:
                colors.append("red")
            elif row["padj"] < 0.05 and row["coef"] < 0:
                colors.append("blue")
            else:
                colors.append("grey")

        ax.scatter(ct_data["coef"], ct_data["neg_log10_padj"], c=colors, alpha=0.6, s=15)
        ax.axhline(-np.log10(0.05), linestyle="--", color="black", linewidth=0.8)
        ax.axvline(0, linestyle="--", color="black", linewidth=0.8)

        # Label top 10 TFs by significance
        top = ct_data.nsmallest(10, "padj")
        for _, row in top.iterrows():
            if row["padj"] < 0.05:
                ax.text(
                    row["coef"],
                    row["neg_log10_padj"],
                    row["name"],
                    fontsize=6,
                    ha="right",
                    va="bottom",
                )

        ax.set_title(ct, fontsize=10)
        ax.set_xlabel("Coefficient")
        ax.set_ylabel("-log10(FDR)")

    # Hide unused axes
    for idx in range(n, nrows * ncols):
        axes[idx // ncols, idx % ncols].set_visible(False)

    fig.suptitle(f"TF Activity Volcanos - {phenotype} ({sex})", fontsize=14, y=1.02)
    plt.tight_layout()

    out = output_dir / f"tf_volcano_{sex}.pdf"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved volcanos: %s", out)


def epigenetic_highlight(
    tf_df: pd.DataFrame, output_dir: Path, phenotype: str, sex: str
) -> None:
    """Highlight panel for known epigenetic regulators."""
    epi = tf_df[tf_df["name"].isin(EPIGENETIC_TFS)]
    if epi.empty:
        log.warning("No epigenetic TFs found in results (sex=%s)", sex)
        return

    pivot_coef = epi.pivot_table(
        index="name", columns="cell_type", values="coef", aggfunc="first"
    )
    pivot_padj = epi.pivot_table(
        index="name", columns="cell_type", values="padj", aggfunc="first"
    )
    annot = pivot_padj.map(lambda x: "*" if pd.notna(x) and x < 0.05 else "")

    fig_height = max(4, len(pivot_coef) * 0.4)
    fig_width = max(6, len(pivot_coef.columns) * 1.2 + 2)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    sns.heatmap(
        pivot_coef,
        cmap="RdBu_r",
        center=0,
        annot=annot,
        fmt="",
        linewidths=0.5,
        ax=ax,
        cbar_kws={"label": f"Coefficient ({phenotype})"},
    )
    ax.set_title(f"Epigenetic TF Activity vs {phenotype} ({sex})\n* = FDR < 0.05")
    ax.set_ylabel("Epigenetic Regulator")
    ax.set_xlabel("Cell Type")
    plt.tight_layout()

    out = output_dir / f"epigenetic_tf_{sex}.pdf"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    log.info("Saved epigenetic panel: %s", out)


def hub_tf_barplot(
    tf_df: pd.DataFrame, output_dir: Path, phenotype: str, sex: str
) -> None:
    """Barplot of TFs significant in 3+ cell types."""
    sig = tf_df[tf_df["padj"] < 0.05]
    if sig.empty:
        log.warning("No significant TFs for hub barplot (sex=%s)", sex)
        return

    counts = sig.groupby("name")["cell_type"].nunique()
    hubs = counts[counts >= 3].sort_values(ascending=False)

    if hubs.empty:
        log.info("No hub TFs (significant in 3+ cell types) for sex=%s", sex)
        return

    fig, ax = plt.subplots(figsize=(max(6, len(hubs) * 0.5), 5))
    hubs.plot(kind="bar", ax=ax, color="steelblue", edgecolor="black")
    ax.set_ylabel("Number of Cell Types (FDR < 0.05)")
    ax.set_xlabel("Transcription Factor")
    ax.set_title(f"Hub TFs: Significant in 3+ Cell Types\n{phenotype} ({sex})")
    ax.axhline(3, linestyle="--", color="red", linewidth=0.8, label="Threshold = 3")
    ax.legend()
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()

    out = output_dir / f"hub_tf_barplot_{sex}.pdf"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    log.info("Saved hub TF barplot: %s", out)


def pathway_heatmap(
    pw_df: pd.DataFrame, output_dir: Path, phenotype: str, sex: str
) -> None:
    """Heatmap of PROGENy pathway coefficients across cell types."""
    if pw_df.empty:
        return

    pivot_coef = pw_df.pivot_table(
        index="name", columns="cell_type", values="coef", aggfunc="first"
    )
    pivot_padj = pw_df.pivot_table(
        index="name", columns="cell_type", values="padj", aggfunc="first"
    )
    annot = pivot_padj.map(lambda x: "*" if pd.notna(x) and x < 0.05 else "")

    fig_height = max(5, len(pivot_coef) * 0.5)
    fig_width = max(6, len(pivot_coef.columns) * 1.2 + 2)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    sns.heatmap(
        pivot_coef,
        cmap="RdBu_r",
        center=0,
        annot=annot,
        fmt="",
        linewidths=0.5,
        ax=ax,
        cbar_kws={"label": f"Coefficient ({phenotype})"},
    )
    ax.set_title(f"Pathway Activity (PROGENy) vs {phenotype} ({sex})\n* = FDR < 0.05")
    ax.set_ylabel("Pathway")
    ax.set_xlabel("Cell Type")
    plt.tight_layout()

    out = output_dir / f"pathway_heatmap_{sex}.pdf"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    log.info("Saved pathway heatmap: %s", out)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Visualize ACE TF/pathway activity results."
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        required=True,
        help="Directory containing tf_summary.csv and optionally pathway_summary.csv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory for output figures",
    )
    parser.add_argument(
        "--phenotype",
        type=str,
        required=True,
        help="ACE phenotype label for titles",
    )
    parser.add_argument(
        "--sex",
        type=str,
        choices=["Male", "Female", "Both"],
        default="Both",
        help="Which sex to visualize (default: Both, generating separate plots per sex)",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    log.info("=" * 60)
    log.info("ACE TF Activity Visualization")
    log.info("  Results dir: %s", args.results_dir)
    log.info("  Output dir:  %s", args.output_dir)
    log.info("  Phenotype:   %s", args.phenotype)
    log.info("  Sex filter:  %s", args.sex)
    log.info("=" * 60)

    # Load TF summary
    tf_path = args.results_dir / "tf_summary.csv"
    if not tf_path.exists():
        log.error("tf_summary.csv not found in %s", args.results_dir)
        sys.exit(1)

    tf_df = pd.read_csv(tf_path)
    log.info("Loaded TF summary: %d rows", len(tf_df))

    # Load pathway summary if present
    pw_path = args.results_dir / "pathway_summary.csv"
    pw_df = pd.read_csv(pw_path) if pw_path.exists() else pd.DataFrame()
    if not pw_df.empty:
        log.info("Loaded pathway summary: %d rows", len(pw_df))

    # Determine which sexes to plot
    if args.sex == "Both":
        sexes = [s for s in ["Male", "Female"] if s in tf_df["sex"].unique()]
    else:
        sexes = [args.sex]

    for sex in sexes:
        log.info("Generating figures for %s...", sex)
        tf_sex = tf_df[tf_df["sex"] == sex]
        pw_sex = pw_df[pw_df["sex"] == sex] if not pw_df.empty else pd.DataFrame()

        if tf_sex.empty:
            log.warning("No TF data for sex=%s, skipping", sex)
            continue

        cross_celltype_heatmap(tf_sex, args.output_dir, args.phenotype, sex)
        per_celltype_volcano(tf_sex, args.output_dir, args.phenotype, sex)
        epigenetic_highlight(tf_sex, args.output_dir, args.phenotype, sex)
        hub_tf_barplot(tf_sex, args.output_dir, args.phenotype, sex)

        if not pw_sex.empty:
            pathway_heatmap(pw_sex, args.output_dir, args.phenotype, sex)

    log.info("Visualization complete.")


if __name__ == "__main__":
    main()
