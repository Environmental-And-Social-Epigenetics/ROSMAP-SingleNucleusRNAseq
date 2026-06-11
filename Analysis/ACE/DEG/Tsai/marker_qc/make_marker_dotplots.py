#!/usr/bin/env python3
"""
Canonical-marker QC dotplots for the Tsai integrated single-nucleus object.

Loads tsai_annotated.h5ad, visualizes expression of canonical lineage markers
across cell-type categories, and writes:
  - figures/dotplot_all_categories.{pdf,png}  (one dotplot, brackets per category)
  - figures/dotplot_<category>.pdf            (per-category dotplots)
  - figures/stacked_violin_<category>.pdf
  - figures/stacked_violin_all_categories.pdf
  - tables/marker_fraction_expressing.csv
  - tables/marker_mean_expression.csv
  - tables/marker_flags.csv                    (missing-own + unexpected-other)
  - tables/markers_present.csv                 (canonical markers present in raw)
  - tables/run_metadata.json
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime
from pathlib import Path

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
sys.stdout.reconfigure(line_buffering=True)

import anndata as ad  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
import scipy.sparse as sp  # noqa: E402

# Mentor's canonical marker panel. Keep gene order intra-category as given.
MARKERS = {
    "PV_inhibitory": ["PVALB", "GAD1", "GAD2", "SLC32A1", "SOX6", "ERBB4", "KCNS3"],
    "Microglia": ["CX3CR1", "P2RY12", "TMEM119", "AIF1", "CSF1R", "C3", "TYROBP", "TREM2", "APOE"],
    "Oligo_myelin": ["PLP1", "MBP", "MOG", "MAG", "MOBP", "CNP", "MYRF", "CLDN11"],
    "Astrocyte": ["AQP4", "GFAP", "SLC1A2", "SLC1A3", "ALDH1L1", "ATP1A2", "KCNJ10"],
    "Excitatory": ["SLC17A7", "CAMK2A", "SATB2", "CUX2", "RORB", "TBR1"],
}

# Deterministic cell-type column order for plots. Includes both the raw H5AD
# labels (slashes/parens, used in the integrated object) and the cleaned forms
# (underscores, used downstream by the DEG pipeline) so the script works
# against either annotation.
CELLTYPE_ORDER = [
    "Ex-L2/3", "Ex-L2_3",
    "Ex-L4",
    "Ex-L4/5", "Ex-L4_5",
    "Ex-L5",
    "Ex-L5/6", "Ex-L5_6",
    "Ex-L5/6-CC", "Ex-L5_6-CC",
    "Ex-NRGN",
    "In-PV (Basket)", "In-PV_Basket",
    "In-PV (Chandelier)", "In-PV_Chandelier",
    "In-Rosehip", "In-SST", "In-VIP",
    "Ast", "Oli", "OPC", "Mic", "Endo",
]

# Which canonical-marker category each cell-type label "should" match. Keyed by
# both raw and cleaned forms.
EXPECTED_CATEGORY = {
    "Ex-L2/3": "Excitatory", "Ex-L2_3": "Excitatory",
    "Ex-L4": "Excitatory",
    "Ex-L4/5": "Excitatory", "Ex-L4_5": "Excitatory",
    "Ex-L5": "Excitatory",
    "Ex-L5/6": "Excitatory", "Ex-L5_6": "Excitatory",
    "Ex-L5/6-CC": "Excitatory", "Ex-L5_6-CC": "Excitatory",
    "Ex-NRGN": "Excitatory",
    "In-PV (Basket)": "PV_inhibitory", "In-PV_Basket": "PV_inhibitory",
    "In-PV (Chandelier)": "PV_inhibitory", "In-PV_Chandelier": "PV_inhibitory",
    "In-Rosehip": "PV_inhibitory", "In-SST": "PV_inhibitory", "In-VIP": "PV_inhibitory",
    "Mic": "Microglia",
    "Ast": "Astrocyte",
    "Oli": "Oligo_myelin", "OPC": "Oligo_myelin",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Marker-QC dotplots for Tsai snRNA.")
    parser.add_argument("--input-h5ad", type=Path, default=None,
                        help="Annotated H5AD. Defaults to ${TSAI_INTEGRATED}/tsai_annotated.h5ad.")
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent,
                        help="Output root containing figures/ and tables/.")
    parser.add_argument("--per-celltype-cap", type=int, default=20_000,
                        help="Max cells per cell type to retain for plotting.")
    parser.add_argument("--threshold", type=float, default=0.30,
                        help="Fraction-expressing threshold for marker flags.")
    parser.add_argument("--seed", type=int, default=0)
    return parser.parse_args()


def resolve_input(args: argparse.Namespace) -> Path:
    if args.input_h5ad is not None:
        return args.input_h5ad.resolve()
    tsai_integrated = os.environ.get("TSAI_INTEGRATED")
    if not tsai_integrated:
        raise SystemExit("ERROR: --input-h5ad not provided and TSAI_INTEGRATED env not set.")
    return (Path(tsai_integrated) / "tsai_annotated.h5ad").resolve()


def per_celltype_subsample(adata: ad.AnnData, cap: int, seed: int) -> ad.AnnData:
    rng = np.random.default_rng(seed)
    keep_idx: list[int] = []
    cell_types = adata.obs["cell_type"].astype(str)
    for ct, idx in cell_types.groupby(cell_types).groups.items():
        idx_arr = np.asarray(idx)
        n_take = min(cap, len(idx_arr))
        if n_take < len(idx_arr):
            chosen = rng.choice(idx_arr, size=n_take, replace=False)
        else:
            chosen = idx_arr
        keep_idx.extend(chosen.tolist())
        print(f"  {ct}: {len(idx_arr)} -> {n_take}")
    keep_mask = np.zeros(adata.n_obs, dtype=bool)
    pos_lookup = pd.Series(np.arange(adata.n_obs), index=adata.obs_names)
    keep_mask[pos_lookup.loc[keep_idx].values] = True
    return adata[keep_mask].copy()


def fraction_expressing_matrix(adata_full: ad.AnnData,
                               genes_present: list[str]) -> pd.DataFrame:
    """Per-cell-type fraction of cells with expression > 0 for each marker."""
    if adata_full.raw is None:
        raise SystemExit("ERROR: adata.raw is None; cannot compute fractions on full genes.")
    raw = adata_full.raw
    gene_idx = pd.Index(raw.var_names).get_indexer(genes_present)
    if (gene_idx < 0).any():
        missing = [g for g, i in zip(genes_present, gene_idx) if i < 0]
        raise SystemExit(f"ERROR: genes missing from adata.raw.var_names: {missing}")
    mat = raw.X[:, gene_idx]
    if sp.issparse(mat):
        nonzero_per_cell = mat > 0
    else:
        nonzero_per_cell = mat > 0
    ct = adata_full.obs["cell_type"].astype(str).values
    rows = []
    for celltype in pd.Index(ct).unique():
        mask = ct == celltype
        if not mask.any():
            continue
        if sp.issparse(nonzero_per_cell):
            frac = np.asarray(nonzero_per_cell[mask].mean(axis=0)).ravel()
        else:
            frac = nonzero_per_cell[mask].mean(axis=0)
        rows.append(pd.Series(frac, index=genes_present, name=celltype))
    return pd.DataFrame(rows)


def mean_expression_matrix(adata_full: ad.AnnData,
                           genes_present: list[str]) -> pd.DataFrame:
    raw = adata_full.raw
    gene_idx = pd.Index(raw.var_names).get_indexer(genes_present)
    mat = raw.X[:, gene_idx]
    ct = adata_full.obs["cell_type"].astype(str).values
    rows = []
    for celltype in pd.Index(ct).unique():
        mask = ct == celltype
        if not mask.any():
            continue
        if sp.issparse(mat):
            mean_vec = np.asarray(mat[mask].mean(axis=0)).ravel()
        else:
            mean_vec = mat[mask].mean(axis=0)
        rows.append(pd.Series(mean_vec, index=genes_present, name=celltype))
    return pd.DataFrame(rows)


def flag_table(frac_df: pd.DataFrame, threshold: float) -> pd.DataFrame:
    gene_to_cat: dict[str, str] = {}
    for cat, genes in MARKERS.items():
        for g in genes:
            gene_to_cat[g] = cat

    rows: list[dict] = []
    for celltype in frac_df.index:
        expected = EXPECTED_CATEGORY.get(celltype)
        for gene in frac_df.columns:
            frac = float(frac_df.loc[celltype, gene])
            cat = gene_to_cat.get(gene, "unknown")
            if expected is not None and cat == expected and frac < threshold:
                rows.append({
                    "cell_type": celltype, "gene": gene, "category": cat,
                    "fraction_expressing": frac, "flag_type": "MISSING_OWN",
                })
            if expected is not None and cat != expected and frac > threshold:
                rows.append({
                    "cell_type": celltype, "gene": gene, "category": cat,
                    "fraction_expressing": frac, "flag_type": "UNEXPECTED_OTHER",
                })
    if not rows:
        return pd.DataFrame(columns=[
            "cell_type", "gene", "category", "fraction_expressing", "flag_type",
        ])
    return pd.DataFrame(rows).sort_values(
        ["flag_type", "cell_type", "fraction_expressing"],
        ascending=[True, True, False],
    ).reset_index(drop=True)


def main() -> None:
    args = parse_args()
    out_dir = args.output_dir.resolve()
    fig_dir = out_dir / "figures"
    tab_dir = out_dir / "tables"
    fig_dir.mkdir(parents=True, exist_ok=True)
    tab_dir.mkdir(parents=True, exist_ok=True)

    input_path = resolve_input(args)
    print(f"Input H5AD: {input_path}")
    if not input_path.exists():
        raise SystemExit(f"ERROR: input H5AD not found: {input_path}")

    print("Loading adata...")
    adata = ad.read_h5ad(input_path)
    print(f"  n_obs={adata.n_obs:,}  n_vars={adata.n_vars:,}")
    if "cell_type" not in adata.obs.columns:
        raise SystemExit("ERROR: obs['cell_type'] missing.")
    print("  cell_type counts:")
    print(adata.obs["cell_type"].value_counts().to_string())
    print(f"  adata.raw is not None: {adata.raw is not None}")
    if adata.raw is not None:
        print(f"  adata.raw.n_vars={adata.raw.n_vars:,}")

    # Intersect markers with raw var_names.
    marker_var = adata.raw.var_names if adata.raw is not None else adata.var_names
    present_rows: list[dict] = []
    markers_present: dict[str, list[str]] = {}
    for cat, genes in MARKERS.items():
        present = [g for g in genes if g in marker_var]
        absent = [g for g in genes if g not in marker_var]
        markers_present[cat] = present
        for g in genes:
            present_rows.append({
                "category": cat, "gene": g, "present_in_raw": g in marker_var,
            })
        if absent:
            print(f"  [{cat}] absent from raw: {absent}")
    pd.DataFrame(present_rows).to_csv(tab_dir / "markers_present.csv", index=False)

    # Constrain cell_type to known categories and the deterministic order.
    obs_ct = adata.obs["cell_type"].astype(str)
    present_categories = [ct for ct in CELLTYPE_ORDER if ct in obs_ct.unique()]
    extra = sorted(set(obs_ct.unique()) - set(CELLTYPE_ORDER))
    final_order = present_categories + extra
    if extra:
        print(f"  WARN: extra cell types not in CELLTYPE_ORDER appended: {extra}")
    adata.obs["cell_type"] = pd.Categorical(obs_ct, categories=final_order, ordered=True)

    # Subsample per cell type for plotting.
    print(f"Subsampling per cell type (cap={args.per_celltype_cap}, seed={args.seed})...")
    adata_sub = per_celltype_subsample(adata, args.per_celltype_cap, args.seed)
    print(f"  subsampled n_obs={adata_sub.n_obs:,}")

    # Primary dotplot (dict form -> auto-brackets per category).
    print("Rendering primary dotplot...")
    dp_all = sc.pl.dotplot(
        adata_sub,
        var_names=markers_present,
        groupby="cell_type",
        use_raw=True,
        standard_scale="var",
        dendrogram=False,
        return_fig=True,
        show=False,
    )
    dp_all.savefig(fig_dir / "dotplot_all_categories.pdf", bbox_inches="tight")
    dp_all.savefig(fig_dir / "dotplot_all_categories.png", bbox_inches="tight", dpi=200)
    plt.close("all")

    dp_all_swap = sc.pl.dotplot(
        adata_sub,
        var_names=markers_present,
        groupby="cell_type",
        use_raw=True,
        standard_scale="var",
        swap_axes=True,
        dendrogram=False,
        return_fig=True,
        show=False,
    )
    dp_all_swap.savefig(fig_dir / "dotplot_all_categories_swapped.pdf", bbox_inches="tight")
    plt.close("all")

    # Per-category dotplots.
    for cat, genes_present in markers_present.items():
        if not genes_present:
            print(f"  [{cat}] no genes present in raw; skipping per-category plots.")
            continue
        print(f"Rendering per-category dotplot: {cat}")
        dp = sc.pl.dotplot(
            adata_sub,
            var_names=genes_present,
            groupby="cell_type",
            use_raw=True,
            standard_scale="var",
            dendrogram=False,
            return_fig=True,
            show=False,
        )
        dp.savefig(fig_dir / f"dotplot_{cat}.pdf", bbox_inches="tight")
        plt.close("all")

        print(f"Rendering stacked-violin: {cat}")
        sv = sc.pl.stacked_violin(
            adata_sub,
            var_names=genes_present,
            groupby="cell_type",
            use_raw=True,
            dendrogram=False,
            return_fig=True,
            show=False,
        )
        sv.savefig(fig_dir / f"stacked_violin_{cat}.pdf", bbox_inches="tight")
        plt.close("all")

    print("Rendering combined stacked-violin...")
    sv_all = sc.pl.stacked_violin(
        adata_sub,
        var_names=markers_present,
        groupby="cell_type",
        use_raw=True,
        dendrogram=False,
        return_fig=True,
        show=False,
    )
    sv_all.savefig(fig_dir / "stacked_violin_all_categories.pdf", bbox_inches="tight")
    plt.close("all")

    # Fraction-expressing + mean-expression matrices (use FULL adata).
    # Preserve the panel's intra-category gene order across categories.
    ordered_present: list[str] = []
    for cat in MARKERS:
        for g in markers_present[cat]:
            if g not in ordered_present:
                ordered_present.append(g)

    print("Computing fraction-expressing matrix on full adata...")
    frac_df = fraction_expressing_matrix(adata, ordered_present)
    frac_df = frac_df.reindex([ct for ct in final_order if ct in frac_df.index])
    frac_df.to_csv(tab_dir / "marker_fraction_expressing.csv")

    print("Computing mean-expression matrix on full adata...")
    mean_df = mean_expression_matrix(adata, ordered_present)
    mean_df = mean_df.reindex([ct for ct in final_order if ct in mean_df.index])
    mean_df.to_csv(tab_dir / "marker_mean_expression.csv")

    print(f"Building flag table (threshold={args.threshold})...")
    flags_df = flag_table(frac_df, args.threshold)
    flags_df.to_csv(tab_dir / "marker_flags.csv", index=False)
    print(f"  n_flags total: {len(flags_df)}")
    if len(flags_df) > 0:
        print(flags_df.groupby("flag_type").size().to_string())

    # Provenance.
    meta = {
        "input_h5ad": str(input_path),
        "output_dir": str(out_dir),
        "n_cells_full": int(adata.n_obs),
        "n_cells_subsample": int(adata_sub.n_obs),
        "per_celltype_cap": args.per_celltype_cap,
        "threshold": args.threshold,
        "seed": args.seed,
        "scanpy_version": sc.__version__,
        "anndata_version": ad.__version__,
        "celltype_order": final_order,
        "markers_present": markers_present,
        "timestamp": datetime.utcnow().isoformat() + "Z",
    }
    (tab_dir / "run_metadata.json").write_text(json.dumps(meta, indent=2))

    print(f"\nDone. Figures: {fig_dir}\n      Tables:  {tab_dir}")


if __name__ == "__main__":
    main()
