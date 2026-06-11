#!/usr/bin/env python3
"""
ACE TF Convergence Analysis

Lightweight post-hoc script to identify hub transcription factors that are
consistently associated with ACE phenotypes across multiple cell types.

Identifies:
  - Hub TFs: significant (FDR < 0.05) in 3+ cell types with consistent
    direction of effect
  - Cross-database validation: TFs significant in both DoRothEA and CollecTRI
  - Optional SCENIC cross-reference if results are available

Outputs:
  - hub_tf_table.csv: hub TFs with cell type counts and mean coefficients
  - convergence_heatmap.pdf: heatmap of hub TF coefficients across cell types

Designed for interactive use (no SLURM needed).

Usage:
    python tf_convergence_analysis.py \
        --results-dir ${ACE_OUTPUT_ROOT}/TFActivity/Tsai/results_derived_batch/tot_adverse_exp \
        --output-dir  ${ACE_OUTPUT_ROOT}/TFActivity/Tsai/convergence_derived_batch/tot_adverse_exp \
        --phenotype tot_adverse_exp
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

# Minimum number of cell types for a TF to be considered a hub
HUB_THRESHOLD = 3


def identify_hubs(
    tf_df: pd.DataFrame,
    threshold: int = HUB_THRESHOLD,
) -> pd.DataFrame:
    """Identify hub TFs: significant in *threshold*+ cell types, same direction.

    Returns a DataFrame with one row per hub TF/sex/source_db combination.
    """
    sig = tf_df[tf_df["padj"] < 0.05].copy()
    if sig.empty:
        return pd.DataFrame()

    results = []
    for (name, sex, source_db), grp in sig.groupby(["name", "sex", "source_db"]):
        n_ct = grp["cell_type"].nunique()
        if n_ct < threshold:
            continue

        # Check direction consistency
        directions = grp["coef"].apply(lambda x: "up" if x > 0 else "down")
        dominant_dir = directions.mode().iloc[0]
        consistent = (directions == dominant_dir).sum()
        consistency_pct = consistent / len(directions) * 100

        results.append(
            {
                "name": name,
                "sex": sex,
                "source_db": source_db,
                "n_celltypes_sig": n_ct,
                "cell_types": ", ".join(sorted(grp["cell_type"].unique())),
                "mean_coef": grp["coef"].mean(),
                "median_coef": grp["coef"].median(),
                "min_padj": grp["padj"].min(),
                "direction": dominant_dir,
                "direction_consistency_pct": round(consistency_pct, 1),
            }
        )

    if not results:
        return pd.DataFrame()

    return (
        pd.DataFrame(results)
        .sort_values(["n_celltypes_sig", "min_padj"], ascending=[False, True])
        .reset_index(drop=True)
    )


def cross_database_validation(tf_df: pd.DataFrame) -> pd.DataFrame:
    """Find TFs significant in both DoRothEA and CollecTRI for same cell type/sex."""
    sig = tf_df[tf_df["padj"] < 0.05].copy()
    if sig.empty:
        return pd.DataFrame()

    dorothea = sig[sig["source_db"] == "DoRothEA"][["name", "cell_type", "sex", "coef", "padj"]]
    collectri = sig[sig["source_db"] == "CollecTRI"][["name", "cell_type", "sex", "coef", "padj"]]

    if dorothea.empty or collectri.empty:
        return pd.DataFrame()

    merged = dorothea.merge(
        collectri,
        on=["name", "cell_type", "sex"],
        suffixes=("_dorothea", "_collectri"),
    )
    merged["same_direction"] = (merged["coef_dorothea"] * merged["coef_collectri"]) > 0
    return merged.sort_values("padj_dorothea")


def load_scenic_results(scenic_dir: Path) -> pd.DataFrame | None:
    """Try to load SCENIC AUCell results for cross-reference."""
    if not scenic_dir.exists():
        return None

    import glob

    csvs = sorted(scenic_dir.glob("*pVals*.csv"))
    if not csvs:
        return None

    dfs = []
    for f in csvs:
        try:
            df = pd.read_csv(f)
            if "regulon" in df.columns:
                # Extract TF name from regulon (e.g. "STAT3(+)" -> "STAT3")
                df["tf_name"] = df["regulon"].str.replace(r"\([+-]\)$", "", regex=True)
                df["source_file"] = f.name
                dfs.append(df)
        except Exception:
            continue

    if not dfs:
        return None

    return pd.concat(dfs, ignore_index=True)


def convergence_heatmap(
    hub_df: pd.DataFrame,
    tf_df: pd.DataFrame,
    output_dir: Path,
    phenotype: str,
) -> None:
    """Generate convergence heatmap for hub TFs."""
    if hub_df.empty:
        log.info("No hub TFs to plot.")
        return

    for sex in hub_df["sex"].unique():
        sex_hubs = hub_df[hub_df["sex"] == sex]
        hub_names = sex_hubs["name"].unique()

        # Get coefficients for hub TFs across all cell types
        tf_sex = tf_df[(tf_df["sex"] == sex) & (tf_df["name"].isin(hub_names))]
        if tf_sex.empty:
            continue

        # Use DoRothEA results preferentially for the heatmap
        dorothea_data = tf_sex[tf_sex["source_db"] == "DoRothEA"]
        if dorothea_data.empty:
            dorothea_data = tf_sex

        pivot_coef = dorothea_data.pivot_table(
            index="name", columns="cell_type", values="coef", aggfunc="first"
        )
        pivot_padj = dorothea_data.pivot_table(
            index="name", columns="cell_type", values="padj", aggfunc="first"
        )
        annot = pivot_padj.map(lambda x: "*" if pd.notna(x) and x < 0.05 else "")

        fig_height = max(5, len(pivot_coef) * 0.5)
        fig_width = max(6, len(pivot_coef.columns) * 1.2 + 3)

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
        ax.set_title(
            f"Hub TFs (sig. in {HUB_THRESHOLD}+ cell types) vs {phenotype} ({sex})\n* = FDR < 0.05"
        )
        ax.set_ylabel("Transcription Factor")
        ax.set_xlabel("Cell Type")
        plt.tight_layout()

        out = output_dir / f"convergence_heatmap_{sex}.pdf"
        fig.savefig(out, dpi=150)
        plt.close(fig)
        log.info("Saved convergence heatmap: %s", out)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="ACE TF convergence analysis across cell types."
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        required=True,
        help="Directory containing tf_summary.csv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory for convergence outputs",
    )
    parser.add_argument(
        "--phenotype",
        type=str,
        required=True,
        help="ACE phenotype label (for titles and file names)",
    )
    parser.add_argument(
        "--scenic-dir",
        type=Path,
        default=None,
        help="Optional: path to SCENIC results for cross-reference",
    )
    parser.add_argument(
        "--hub-threshold",
        type=int,
        default=HUB_THRESHOLD,
        help=f"Minimum cell types for hub TF (default: {HUB_THRESHOLD})",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    log.info("=" * 60)
    log.info("ACE TF Convergence Analysis")
    log.info("  Results dir:   %s", args.results_dir)
    log.info("  Output dir:    %s", args.output_dir)
    log.info("  Phenotype:     %s", args.phenotype)
    log.info("  Hub threshold: %d", args.hub_threshold)
    log.info("=" * 60)

    # Load TF summary
    tf_path = args.results_dir / "tf_summary.csv"
    if not tf_path.exists():
        log.error("tf_summary.csv not found in %s", args.results_dir)
        sys.exit(1)

    tf_df = pd.read_csv(tf_path)
    log.info("Loaded TF summary: %d rows", len(tf_df))

    # Identify hub TFs
    hub_df = identify_hubs(tf_df, threshold=args.hub_threshold)
    if hub_df.empty:
        log.info("No hub TFs found at threshold=%d.", args.hub_threshold)
    else:
        hub_path = args.output_dir / "hub_tf_table.csv"
        hub_df.to_csv(hub_path, index=False)
        log.info("Hub TFs: %d entries -> %s", len(hub_df), hub_path)

        # Print summary
        for sex in hub_df["sex"].unique():
            sex_hubs = hub_df[hub_df["sex"] == sex]
            log.info("  %s: %d hub TF/db combos", sex, len(sex_hubs))
            for _, row in sex_hubs.head(10).iterrows():
                log.info(
                    "    %s (%s): %d cell types, mean coef=%.4f, dir=%s (%.0f%% consistent)",
                    row["name"],
                    row["source_db"],
                    row["n_celltypes_sig"],
                    row["mean_coef"],
                    row["direction"],
                    row["direction_consistency_pct"],
                )

    # Cross-database validation
    xdb = cross_database_validation(tf_df)
    if not xdb.empty:
        xdb_path = args.output_dir / "cross_database_validation.csv"
        xdb.to_csv(xdb_path, index=False)
        n_consistent = xdb["same_direction"].sum()
        log.info(
            "Cross-database validation: %d TF-celltype pairs validated, "
            "%d (%.0f%%) same direction",
            len(xdb),
            n_consistent,
            n_consistent / len(xdb) * 100 if len(xdb) > 0 else 0,
        )

    # SCENIC cross-reference
    scenic_dir = args.scenic_dir
    if scenic_dir is None:
        # Try default SCENIC location
        # Infer from results-dir: .../TFActivity/Tsai/... -> .../SCENIC/Tsai/
        try:
            base = args.results_dir
            while base.name and base.name != "TFActivity":
                base = base.parent
            scenic_dir = base.parent / "SCENIC" / "Tsai"
        except Exception:
            scenic_dir = None

    if scenic_dir is not None:
        scenic_df = load_scenic_results(scenic_dir)
        if scenic_df is not None and not hub_df.empty:
            hub_names = set(hub_df["name"].unique())
            scenic_overlap = scenic_df[scenic_df["tf_name"].isin(hub_names)]
            if not scenic_overlap.empty:
                overlap_path = args.output_dir / "scenic_overlap.csv"
                scenic_overlap.to_csv(overlap_path, index=False)
                log.info(
                    "SCENIC overlap: %d rows for %d hub TFs",
                    len(scenic_overlap),
                    scenic_overlap["tf_name"].nunique(),
                )
        elif scenic_df is None:
            log.info("No SCENIC results found at %s", scenic_dir)

    # Generate convergence heatmap
    if not hub_df.empty:
        convergence_heatmap(hub_df, tf_df, args.output_dir, args.phenotype)

    log.info("Convergence analysis complete.")


if __name__ == "__main__":
    main()
