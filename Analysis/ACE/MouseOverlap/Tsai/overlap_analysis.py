#!/usr/bin/env python3
"""Compute direction-aware mouse LNB vs human ACE DEG overlap."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats


MOUSE_DIRECTIONS = ("LNB_up", "LNB_down")
HUMAN_DIRECTIONS = ("ACE_up", "ACE_down")
DEFAULT_HUMAN_TARGETS = ("Exc", "Inh", "In-PV_Basket", "In-PV_Chandelier")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mouse-deg", type=Path, required=True)
    parser.add_argument("--human-deg", type=Path, required=True)
    parser.add_argument("--ortholog-table", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--figures-dir", type=Path, required=True)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--human-targets", default=",".join(DEFAULT_HUMAN_TARGETS))
    parser.add_argument("--one-to-one-only", action="store_true", default=True)
    return parser.parse_args()


def split_csv(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def bh_adjust(pvalues: pd.Series) -> pd.Series:
    p = pd.to_numeric(pvalues, errors="coerce").to_numpy(dtype=float)
    out = np.full(len(p), np.nan)
    valid = np.isfinite(p)
    if not valid.any():
        return pd.Series(out, index=pvalues.index)
    idx = np.where(valid)[0]
    order = idx[np.argsort(p[idx])]
    ranked = p[order]
    n = len(ranked)
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    out[order] = adjusted
    return pd.Series(out, index=pvalues.index)


def deduplicate_by_best_padj(df: pd.DataFrame, group_cols: list[str], padj_col: str, effect_col: str) -> pd.DataFrame:
    work = df.copy()
    work["_padj_sort"] = pd.to_numeric(work[padj_col], errors="coerce").fillna(np.inf)
    work["_abs_effect_sort"] = pd.to_numeric(work[effect_col], errors="coerce").abs().fillna(-np.inf)
    work = work.sort_values(group_cols + ["_padj_sort", "_abs_effect_sort"], ascending=[True] * len(group_cols) + [True, False])
    work = work.drop_duplicates(group_cols, keep="first")
    return work.drop(columns=["_padj_sort", "_abs_effect_sort"])


def comparison_role(mouse_population: str, human_cell_type: str) -> str:
    if mouse_population == "NeuN" and human_cell_type == "Exc":
        return "primary"
    if mouse_population == "PV" and human_cell_type == "Inh":
        return "primary"
    if mouse_population == "PV" and human_cell_type in {"In-PV_Basket", "In-PV_Chandelier"}:
        return "pv_sensitivity"
    return "exploratory_off_target"


def signed_concordance(x: pd.Series, y: pd.Series) -> dict[str, float]:
    data = pd.DataFrame({"x": pd.to_numeric(x, errors="coerce"), "y": pd.to_numeric(y, errors="coerce")}).dropna()
    if len(data) < 3:
        return {
            "n_ranked": len(data),
            "spearman_r": np.nan,
            "spearman_pvalue": np.nan,
            "sign_concordance": np.nan,
            "signed_product_mean": np.nan,
            "signed_product_wilcoxon_pvalue": np.nan,
        }

    spearman = stats.spearmanr(data["x"], data["y"])
    sx = np.sign(data["x"])
    sy = np.sign(data["y"])
    nonzero = (sx != 0) & (sy != 0)
    sign_conc = float((sx[nonzero] == sy[nonzero]).mean()) if nonzero.any() else np.nan

    xz = (data["x"] - data["x"].mean()) / data["x"].std(ddof=0)
    yz = (data["y"] - data["y"].mean()) / data["y"].std(ddof=0)
    product = (xz * yz).replace([np.inf, -np.inf], np.nan).dropna()
    if len(product) >= 3 and product.abs().sum() > 0:
        try:
            wilcoxon_p = stats.wilcoxon(product, alternative="greater", zero_method="wilcox").pvalue
        except ValueError:
            wilcoxon_p = np.nan
    else:
        wilcoxon_p = np.nan

    return {
        "n_ranked": int(len(data)),
        "spearman_r": float(spearman.statistic),
        "spearman_pvalue": float(spearman.pvalue),
        "sign_concordance": sign_conc,
        "signed_product_mean": float(product.mean()) if len(product) else np.nan,
        "signed_product_wilcoxon_pvalue": wilcoxon_p,
    }


def load_inputs(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    mouse = pd.read_csv(args.mouse_deg)
    human = pd.read_csv(args.human_deg)
    orth = pd.read_csv(args.ortholog_table)

    if args.one_to_one_only and "is_one_to_one" in orth.columns:
        orth = orth[orth["is_one_to_one"].astype(str).str.lower().isin(["true", "1"])]

    mouse = mouse.merge(orth, on="mouse_symbol", how="inner")
    mouse = mouse.rename(columns={"human_symbol": "ortholog_human_symbol"})
    mouse = deduplicate_by_best_padj(
        mouse,
        ["mouse_dataset", "ortholog_human_symbol"],
        "mouse_padj",
        "mouse_lnb_log2fc",
    )

    human = human.rename(columns={"gene_symbol": "ortholog_human_symbol"})
    human = deduplicate_by_best_padj(
        human,
        ["integration", "phenotype", "cell_type", "sex", "ortholog_human_symbol"],
        "padj",
        "log2FoldChange",
    )

    return mouse, human, orth


def compute_overlap(args: argparse.Namespace, mouse: pd.DataFrame, human: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    targets = set(split_csv(args.human_targets))
    human = human[human["cell_type"].isin(targets)].copy()

    summary_rows: list[dict[str, object]] = []
    gene_rows: list[pd.DataFrame] = []

    human_groups = ["integration", "phenotype", "cell_type", "sex"]
    for mouse_dataset, mouse_df in mouse.groupby("mouse_dataset", sort=True):
        mouse_population = str(mouse_df["mouse_population"].iloc[0])
        for human_key, human_df in human.groupby(human_groups, sort=True):
            integration, phenotype, cell_type, sex = human_key
            role = comparison_role(mouse_population, str(cell_type))

            merged = mouse_df.merge(
                human_df,
                on="ortholog_human_symbol",
                how="inner",
                suffixes=("_mouse", "_human"),
            )
            universe = merged.dropna(subset=["mouse_lnb_log2fc", "log2FoldChange"])
            universe_n = int(len(universe))
            concordance = signed_concordance(universe["mouse_lnb_log2fc"], universe["log2FoldChange"])

            for mouse_dir in MOUSE_DIRECTIONS:
                mouse_set_df = universe[
                    (universe["mouse_sig"].astype(bool)) & (universe["mouse_direction"] == mouse_dir)
                ]
                mouse_genes = set(mouse_set_df["ortholog_human_symbol"])
                for human_dir in HUMAN_DIRECTIONS:
                    human_set_df = universe[
                        (universe["human_sig"].astype(bool)) & (universe["human_direction"] == human_dir)
                    ]
                    human_genes = set(human_set_df["ortholog_human_symbol"])
                    overlap = mouse_genes & human_genes
                    union = mouse_genes | human_genes
                    a = len(overlap)
                    b = len(mouse_genes - human_genes)
                    c = len(human_genes - mouse_genes)
                    d = max(universe_n - a - b - c, 0)

                    if universe_n > 0:
                        oddsratio, fisher_p = stats.fisher_exact([[a, b], [c, d]], alternative="greater")
                        hypergeom_p = stats.hypergeom.sf(a - 1, universe_n, len(human_genes), len(mouse_genes))
                    else:
                        oddsratio, fisher_p, hypergeom_p = np.nan, np.nan, np.nan

                    comparison_id = "|".join(
                        [
                            mouse_dataset,
                            integration,
                            phenotype,
                            str(cell_type),
                            str(sex),
                            mouse_dir,
                            human_dir,
                        ]
                    )
                    row = {
                        "comparison_id": comparison_id,
                        "comparison_role": role,
                        "mouse_dataset": mouse_dataset,
                        "mouse_population": mouse_population,
                        "mouse_condition": mouse_df["condition"].iloc[0],
                        "human_integration": integration,
                        "human_phenotype": phenotype,
                        "human_cell_type": cell_type,
                        "human_sex": sex,
                        "mouse_direction": mouse_dir,
                        "human_direction": human_dir,
                        "universe_n": universe_n,
                        "mouse_set_n": len(mouse_genes),
                        "human_set_n": len(human_genes),
                        "overlap_n": a,
                        "fisher_odds_ratio": oddsratio,
                        "fisher_pvalue": fisher_p,
                        "hypergeom_pvalue": hypergeom_p,
                        "jaccard": a / len(union) if union else 0.0,
                        **concordance,
                    }
                    summary_rows.append(row)

                    if overlap:
                        genes = universe[universe["ortholog_human_symbol"].isin(overlap)].copy()
                        genes["comparison_id"] = comparison_id
                        genes["mouse_overlap_direction"] = mouse_dir
                        genes["human_overlap_direction"] = human_dir
                        gene_rows.append(
                            genes[
                                [
                                    "comparison_id",
                                    "ortholog_human_symbol",
                                    "mouse_symbol",
                                    "mouse_ensembl_id",
                                    "mouse_dataset",
                                    "mouse_population",
                                    "condition",
                                    "mouse_lnb_log2fc",
                                    "mouse_padj",
                                    "mouse_direction",
                                    "integration",
                                    "phenotype",
                                    "cell_type",
                                    "sex",
                                    "log2FoldChange",
                                    "padj",
                                    "human_direction",
                                ]
                            ].rename(
                                columns={
                                    "ortholog_human_symbol": "human_symbol",
                                    "condition": "mouse_condition",
                                    "integration": "human_integration",
                                    "phenotype": "human_phenotype",
                                    "cell_type": "human_cell_type",
                                    "sex": "human_sex",
                                    "log2FoldChange": "human_log2FoldChange",
                                    "padj": "human_padj",
                                }
                            )
                        )

    summary = pd.DataFrame(summary_rows)
    if not summary.empty:
        summary["fisher_padj"] = bh_adjust(summary["fisher_pvalue"])
        summary["hypergeom_padj"] = bh_adjust(summary["hypergeom_pvalue"])
        summary["spearman_padj"] = bh_adjust(summary["spearman_pvalue"])
        summary["signed_product_wilcoxon_padj"] = bh_adjust(summary["signed_product_wilcoxon_pvalue"])
        summary = summary.sort_values(
            ["comparison_role", "human_integration", "human_phenotype", "fisher_pvalue", "comparison_id"],
            na_position="last",
        )

    genes_out = pd.concat(gene_rows, ignore_index=True) if gene_rows else pd.DataFrame()
    return summary, genes_out


def plot_heatmap(summary: pd.DataFrame, figures_dir: Path) -> None:
    primary = summary[
        (summary["human_integration"] == "derived_batch")
        & (summary["human_phenotype"] == "tot_adverse_exp")
    ].copy()
    if primary.empty:
        return
    primary["target"] = primary["human_cell_type"] + "_" + primary["human_sex"].astype(str)
    primary["direction_pair"] = primary["mouse_direction"] + "_vs_" + primary["human_direction"]
    primary["score"] = -np.log10(primary["fisher_padj"].clip(lower=1e-300))

    direction_pairs = [f"{m}_vs_{h}" for m in MOUSE_DIRECTIONS for h in HUMAN_DIRECTIONS]
    fig, axes = plt.subplots(2, 2, figsize=(14, 8), constrained_layout=True)
    for ax, direction_pair in zip(axes.flat, direction_pairs):
        sub = primary[primary["direction_pair"] == direction_pair]
        mat = sub.pivot_table(index="mouse_dataset", columns="target", values="score", aggfunc="max").fillna(0)
        sns.heatmap(mat, ax=ax, cmap="viridis", cbar=True, linewidths=0.3, linecolor="white")
        ax.set_title(direction_pair)
        ax.set_xlabel("")
        ax.set_ylabel("")
    fig.suptitle("Directional overlap strength, derived_batch tot_adverse_exp")
    for ext in ("png", "pdf"):
        fig.savefig(figures_dir / f"directional_overlap_heatmap.{ext}", dpi=300)
    plt.close(fig)


def plot_overlap_bars(summary: pd.DataFrame, figures_dir: Path) -> None:
    if summary.empty:
        return
    sig = summary[summary["fisher_padj"] < 0.1].copy()
    if sig.empty:
        sig = summary.nsmallest(20, "fisher_pvalue").copy()
    sig["label"] = (
        sig["mouse_dataset"]
        + " | "
        + sig["human_integration"]
        + " | "
        + sig["human_phenotype"]
        + " | "
        + sig["human_cell_type"]
        + " "
        + sig["human_sex"].astype(str)
        + " | "
        + sig["mouse_direction"]
        + "/"
        + sig["human_direction"]
    )
    sig = sig.sort_values("overlap_n", ascending=True).tail(25)
    height = max(4, 0.28 * len(sig))
    fig, ax = plt.subplots(figsize=(12, height))
    ax.barh(sig["label"], sig["overlap_n"], color="#4C78A8")
    ax.set_xlabel("Overlapping genes")
    ax.set_ylabel("")
    ax.set_title("Top directional overlap counts")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(figures_dir / f"significant_overlap_barplot.{ext}", dpi=300)
    plt.close(fig)


def plot_scatters(summary: pd.DataFrame, mouse: pd.DataFrame, human: pd.DataFrame, figures_dir: Path) -> None:
    candidates = summary[
        (summary["human_integration"] == "derived_batch")
        & (summary["human_phenotype"] == "tot_adverse_exp")
        & (summary["comparison_role"].isin(["primary", "pv_sensitivity"]))
    ].copy()
    if candidates.empty:
        candidates = summary.copy()
    candidates = candidates.sort_values(["fisher_pvalue", "comparison_id"], na_position="last").head(6)
    if candidates.empty:
        return

    n = len(candidates)
    ncols = 2
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(11, 4.5 * nrows), squeeze=False)
    human_targets = set(candidates["human_cell_type"])
    human_sub = human[human["cell_type"].isin(human_targets)].copy()

    for ax, (_, row) in zip(axes.flat, candidates.iterrows()):
        m = mouse[mouse["mouse_dataset"] == row["mouse_dataset"]]
        h = human_sub[
            (human_sub["integration"] == row["human_integration"])
            & (human_sub["phenotype"] == row["human_phenotype"])
            & (human_sub["cell_type"] == row["human_cell_type"])
            & (human_sub["sex"] == row["human_sex"])
        ]
        merged = m.merge(h, on="ortholog_human_symbol", how="inner")
        ax.scatter(merged["mouse_lnb_log2fc"], merged["log2FoldChange"], s=5, alpha=0.25, color="#4C78A8")
        ax.axhline(0, color="black", linewidth=0.6)
        ax.axvline(0, color="black", linewidth=0.6)
        title = (
            f"{row['mouse_dataset']} vs {row['human_cell_type']} {row['human_sex']}\n"
            f"{row['human_phenotype']} rho={row['spearman_r']:.3f}"
        )
        ax.set_title(title)
        ax.set_xlabel("Mouse LNB log2FC")
        ax.set_ylabel("Human ACE log2FC")
    for ax in axes.flat[n:]:
        ax.axis("off")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(figures_dir / f"signed_effect_scatter_top.{ext}", dpi=300)
    plt.close(fig)


def write_coverage(mouse: pd.DataFrame, human: pd.DataFrame, orth: pd.DataFrame, output_dir: Path) -> None:
    rows = []
    for dataset, sub in mouse.groupby("mouse_dataset", sort=True):
        rows.append(
            {
                "entity": "mouse_dataset",
                "label": dataset,
                "mapped_one_to_one_genes": int(sub["ortholog_human_symbol"].nunique()),
                "rows": int(len(sub)),
            }
        )
    for key, sub in human.groupby(["integration", "phenotype", "cell_type", "sex"], sort=True):
        rows.append(
            {
                "entity": "human_contrast",
                "label": "|".join(map(str, key)),
                "mapped_one_to_one_genes": int(sub["ortholog_human_symbol"].nunique()),
                "rows": int(len(sub)),
            }
        )
    rows.append(
        {
            "entity": "ortholog_table",
            "label": "one_to_one_used",
            "mapped_one_to_one_genes": int(orth[orth["is_one_to_one"].astype(str).str.lower().isin(["true", "1"])]["human_symbol"].nunique())
            if "is_one_to_one" in orth.columns
            else int(orth["human_symbol"].nunique()),
            "rows": int(len(orth)),
        }
    )
    pd.DataFrame(rows).to_csv(output_dir / "ortholog_coverage_summary.csv", index=False)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.figures_dir.mkdir(parents=True, exist_ok=True)

    mouse, human, orth = load_inputs(args)
    summary, genes = compute_overlap(args, mouse, human)

    summary_path = args.output_dir / "overlap_summary.csv"
    genes_path = args.output_dir / "overlap_genes.csv"
    summary.to_csv(summary_path, index=False)
    genes.to_csv(genes_path, index=False)
    write_coverage(mouse, human, orth, args.output_dir)

    plot_heatmap(summary, args.figures_dir)
    plot_overlap_bars(summary, args.figures_dir)
    plot_scatters(summary, mouse, human, args.figures_dir)

    print(f"Wrote: {summary_path}")
    print(f"Wrote: {genes_path}")
    print(f"Summary rows: {len(summary)}")
    print(f"Overlap gene rows: {len(genes)}")


if __name__ == "__main__":
    main()

