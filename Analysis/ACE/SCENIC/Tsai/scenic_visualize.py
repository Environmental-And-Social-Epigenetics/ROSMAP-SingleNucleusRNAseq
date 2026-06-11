#!/usr/bin/env python3
"""
ACE SCENIC Visualization

Reads regression_results.csv files produced by scenic_analysis.py across
all cell-type x sex combinations and generates summary figures:

1. Volcano plot per cell type (coef vs -log10 padj)
2. Cross-cell-type regulon heatmap (significant regulons x cell types)
3. Bar plot of significant regulon counts per cell type

Usage:
    python scenic_visualize.py \
        --results-dir /path/to/results_derived_batch/tot_adverse_exp \
        --output-dir /path/to/figures \
        --phenotype tot_adverse_exp
"""

import argparse
import logging
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

PALETTE = {"Up": "red", "Down": "blue", "Not Significant": "gray"}


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def discover_results(results_dir):
    """
    Walk ``results_dir`` looking for regression_results.csv files.

    Expected layout::

        results_dir/
            Male_Mic/regression_results.csv
            Female_Ast/regression_results.csv
            ...

    Returns a list of (sex, cell_type, DataFrame) tuples.
    """
    records = []
    if not os.path.isdir(results_dir):
        log.error("Results directory not found: %s", results_dir)
        return records

    for entry in sorted(os.listdir(results_dir)):
        csv_path = os.path.join(results_dir, entry, "regression_results.csv")
        if not os.path.isfile(csv_path):
            continue

        parts = entry.split("_", 1)
        if len(parts) != 2:
            log.warning("Skipping unexpected directory name: %s", entry)
            continue

        sex, cell_type = parts
        df = pd.read_csv(csv_path)
        if df.empty:
            log.warning("Empty results for %s/%s -- skipping", sex, cell_type)
            continue

        df["sex"] = sex
        df["cell_type"] = cell_type
        records.append((sex, cell_type, df))

    log.info("Found %d result files in %s", len(records), results_dir)
    return records


# ---------------------------------------------------------------------------
# Volcano plot
# ---------------------------------------------------------------------------

def plot_volcano(df, sex, cell_type, phenotype, ax=None):
    """
    Volcano plot: coefficient vs -log10(padj).
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(10, 6))

    df = df.copy()
    df["neg_log10_padj"] = -np.log10(df["padj"].clip(lower=1e-300))
    df["significance"] = "Not Significant"
    df.loc[(df["padj"] < 0.05) & (df["coef"] > 0), "significance"] = "Up"
    df.loc[(df["padj"] < 0.05) & (df["coef"] < 0), "significance"] = "Down"

    sns.scatterplot(
        data=df, x="coef", y="neg_log10_padj", hue="significance",
        palette=PALETTE, edgecolor=None, ax=ax, legend=True,
    )

    # Label significant regulons
    sig = df[df["significance"] != "Not Significant"]
    for _, row in sig.iterrows():
        ax.text(
            row["coef"], row["neg_log10_padj"], row["regulon"],
            fontsize=6, ha="right", va="bottom", alpha=0.8,
        )

    ax.axhline(-np.log10(0.05), ls="--", color="black", lw=0.8)
    ax.axvline(0, ls="--", color="black", lw=0.8)
    ax.set_xlabel("Effect size (coefficient)")
    ax.set_ylabel("-log10(FDR)")
    ax.set_title(f"{cell_type} ({sex}) -- ACE ({phenotype}) TF Activity")
    ax.legend(title="Significance", loc="upper right", fontsize=7)

    return ax


# ---------------------------------------------------------------------------
# Cross-cell-type heatmap
# ---------------------------------------------------------------------------

def plot_heatmap(all_df, phenotype, ax=None):
    """
    Heatmap of coefficients for regulons significant in at least one cell type.
    Rows = regulons, columns = sex_celltype.
    """
    sig = all_df[all_df["padj"] < 0.05]
    if sig.empty:
        log.warning("No significant regulons for heatmap.")
        return None

    sig_regulons = sig["regulon"].unique()
    subset = all_df[all_df["regulon"].isin(sig_regulons)].copy()
    subset["label"] = subset["sex"] + "_" + subset["cell_type"]

    pivot = subset.pivot_table(
        index="regulon", columns="label", values="coef", aggfunc="first"
    )

    if pivot.empty:
        return None

    if ax is None:
        height = max(4, len(pivot) * 0.3)
        width = max(6, len(pivot.columns) * 0.8)
        _, ax = plt.subplots(figsize=(width, height))

    sns.heatmap(
        pivot, cmap="RdBu_r", center=0, ax=ax,
        linewidths=0.5, linecolor="white",
        cbar_kws={"label": "Coefficient"},
    )
    ax.set_title(f"ACE ({phenotype}) -- Significant Regulon Coefficients")
    ax.set_ylabel("Regulon")
    ax.set_xlabel("")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=8)
    plt.setp(ax.get_yticklabels(), fontsize=7)

    return ax


# ---------------------------------------------------------------------------
# Bar plot: count of significant regulons
# ---------------------------------------------------------------------------

def plot_sig_counts(all_df, phenotype, ax=None):
    """
    Bar plot showing number of significant regulons per cell type x sex.
    """
    sig = all_df[all_df["padj"] < 0.05].copy()
    sig["label"] = sig["sex"] + "_" + sig["cell_type"]

    counts = sig.groupby("label").size().sort_values(ascending=False)

    if counts.empty:
        log.warning("No significant regulons for bar plot.")
        return None

    if ax is None:
        _, ax = plt.subplots(figsize=(max(6, len(counts) * 0.6), 5))

    counts.plot.bar(ax=ax, color="steelblue", edgecolor="black")
    ax.set_xlabel("")
    ax.set_ylabel("Number of significant regulons (padj < 0.05)")
    ax.set_title(f"ACE ({phenotype}) -- Significant Regulon Counts")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", fontsize=9)

    return ax


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description="ACE SCENIC visualization")
    p.add_argument("--results-dir", required=True,
                   help="Root results directory containing Sex_CellType/ subdirs")
    p.add_argument("--output-dir", required=True,
                   help="Directory to save figures")
    p.add_argument("--phenotype", default="tot_adverse_exp",
                   help="Phenotype name (for labels, default: tot_adverse_exp)")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    records = discover_results(args.results_dir)
    if not records:
        log.error("No result files found in %s", args.results_dir)
        sys.exit(1)

    # Combine all results into one DataFrame
    all_df = pd.concat([r[2] for r in records], ignore_index=True)

    pdf_path = os.path.join(
        args.output_dir,
        f"scenic_summary_{args.phenotype}.pdf",
    )
    log.info("Generating summary PDF: %s", pdf_path)

    with PdfPages(pdf_path) as pdf:
        # -- Volcano plots per cell type x sex --
        for sex, cell_type, df in records:
            fig, ax = plt.subplots(figsize=(10, 6))
            plot_volcano(df, sex, cell_type, args.phenotype, ax=ax)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

        # -- Cross-cell-type heatmap --
        sig_regulons = all_df[all_df["padj"] < 0.05]["regulon"].unique()
        if len(sig_regulons) > 0:
            height = max(4, len(sig_regulons) * 0.3)
            labels = (all_df["sex"] + "_" + all_df["cell_type"]).nunique()
            width = max(6, labels * 0.8)
            fig, ax = plt.subplots(figsize=(width, height))
            result = plot_heatmap(all_df, args.phenotype, ax=ax)
            if result is not None:
                fig.tight_layout()
                pdf.savefig(fig)
            plt.close(fig)

        # -- Bar plot of sig counts --
        fig, ax = plt.subplots(figsize=(10, 5))
        result = plot_sig_counts(all_df, args.phenotype, ax=ax)
        if result is not None:
            fig.tight_layout()
            pdf.savefig(fig)
        plt.close(fig)

    log.info("Saved summary PDF -> %s", pdf_path)

    # Also save a combined CSV for convenience
    combined_path = os.path.join(
        args.output_dir,
        f"all_regression_results_{args.phenotype}.csv",
    )
    all_df.to_csv(combined_path, index=False)
    log.info("Saved combined results -> %s", combined_path)

    log.info("Visualization complete.")


if __name__ == "__main__":
    main()
