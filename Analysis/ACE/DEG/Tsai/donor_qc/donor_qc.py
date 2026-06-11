#!/usr/bin/env python3
"""
Donor-level + annotation QC for Tsai snRNA cell-type labels.

Loads tsai_annotated.h5ad and walks through six diagnostic sections:
  A. Cohort-level per-donor QC baseline.
  B. Cluster x donor composition.
  C. Annotation confidence (ORA gap top1 vs runner-up).
  D. Pipeline-panel (Mohammadi) marker QC.
  E. Per-donor drill-down on problematic cell types.
  F. Doublet drill-down.

Writes CSVs to tables/ and PDFs/PNGs to figures/. Section G assembles the
markdown report (see build_report.py).
"""

from __future__ import annotations

import argparse
import gc
import json
import os
import sys
import warnings
from datetime import datetime
from pathlib import Path

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
sys.stdout.reconfigure(line_buffering=True)

import anndata as ad  # noqa: E402
import h5py  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
import scipy.sparse as sp  # noqa: E402

# Problematic cell types flagged by the canonical-panel QC. Section E drills
# into each, computing per-donor marker fractions for the markers that failed
# (MISSING_OWN) in the canonical-panel run.
PROBLEM_CELLTYPES = {
    "Mic": ["P2RY12", "CX3CR1", "TYROBP", "AIF1", "TREM2", "TMEM119", "APOE"],
    "OPC": ["PLP1", "MOG", "MYRF", "CLDN11", "MOBP", "MAG"],
    "Ex-NRGN": ["RORB", "SATB2", "TBR1", "CUX2"],
    "Ex-L4": ["CUX2"],
    "Ex-L4/5": ["CUX2"],
    "Ex-L5": ["RORB", "CUX2"],
    "Ex-L5/6": ["RORB", "CUX2"],
    "Ex-L5/6-CC": ["CUX2"],
}

# Canonical-panel set (used for the panel comparison in Section D3).
CANONICAL_MARKERS = {
    "PV_inhibitory": ["PVALB", "GAD1", "GAD2", "SLC32A1", "SOX6", "ERBB4", "KCNS3"],
    "Microglia": ["CX3CR1", "P2RY12", "TMEM119", "AIF1", "CSF1R", "C3", "TYROBP", "TREM2", "APOE"],
    "Oligo_myelin": ["PLP1", "MBP", "MOG", "MAG", "MOBP", "CNP", "MYRF", "CLDN11"],
    "Astrocyte": ["AQP4", "GFAP", "SLC1A2", "SLC1A3", "ALDH1L1", "ATP1A2", "KCNJ10"],
    "Excitatory": ["SLC17A7", "CAMK2A", "SATB2", "CUX2", "RORB", "TBR1"],
}

# Map cell_type -> canonical category (for the panel-comparison summary).
EXPECTED_CATEGORY = {
    "Ex-L2/3": "Excitatory", "Ex-L4": "Excitatory", "Ex-L4/5": "Excitatory",
    "Ex-L5": "Excitatory", "Ex-L5/6": "Excitatory", "Ex-L5/6-CC": "Excitatory",
    "Ex-NRGN": "Excitatory",
    "In-PV (Basket)": "PV_inhibitory", "In-PV (Chandelier)": "PV_inhibitory",
    "In-Rosehip": "PV_inhibitory", "In-SST": "PV_inhibitory", "In-VIP": "PV_inhibitory",
    "Mic": "Microglia",
    "Ast": "Astrocyte",
    "Oli": "Oligo_myelin", "OPC": "Oligo_myelin",
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Donor-level + annotation QC")
    p.add_argument("--input-h5ad", type=Path, default=None,
                   help="Default: ${TSAI_INTEGRATED}/tsai_annotated.h5ad")
    p.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    p.add_argument("--qc-summary-csv", type=Path, default=None)
    p.add_argument("--doublet-summary-csv", type=Path, default=None)
    p.add_argument("--cluster-rankings-csv", type=Path, default=None,
                   help="cluster_annotation_rankings.csv from the integration stage.")
    p.add_argument("--mohammadi-rds", type=Path, default=None)
    p.add_argument("--smoke", action="store_true",
                   help="Subsample to 50k cells and run only sections A1/A2.")
    p.add_argument("--celltype", type=str, default=None,
                   help="If set, run E + F drill-down on this cell type only.")
    p.add_argument("--per-celltype-cap", type=int, default=20_000,
                   help="Per-celltype subsample cap for dotplots.")
    p.add_argument("--drop-step", type=int, default=1,
                   help="Section E3: how many donors to drop per step.")
    p.add_argument("--drop-max", type=int, default=20,
                   help="Section E3: max donors to drop.")
    p.add_argument("--seed", type=int, default=0)
    return p.parse_args()


def resolve_integrated_dir(configured: str) -> Path:
    """Locate the integration output dir containing tsai_annotated.h5ad.

    The repo convention is 03_Integrated/<variant>/ (e.g. .../derived_batch/),
    but existing data may still use the flat 03_Integrated/ layout. Probe, in
    order: the configured path, its parent (flat layout), and a derived_batch
    child (in case the configured path is the flat parent). Returns the first
    candidate that contains tsai_annotated.h5ad; otherwise the configured path
    unchanged (so the caller's not-found error names the configured location).
    """
    cfg = Path(configured)
    candidates = [cfg, cfg.parent, cfg / "derived_batch"]
    seen: set[Path] = set()
    for cand in candidates:
        if cand in seen:
            continue
        seen.add(cand)
        if (cand / "tsai_annotated.h5ad").exists():
            if cand != cfg:
                log(f"  [layout] integrated dir resolved to {cand} (configured: {cfg})")
            return cand
    return cfg


def resolve_paths(args: argparse.Namespace) -> dict:
    repo = Path(__file__).resolve().parents[5]
    tsai_processing = repo / "Data" / "Transcriptomics" / "Tsai" / "Processing_Outputs"
    integrated_cfg = os.environ.get("TSAI_INTEGRATED", str(tsai_processing / "03_Integrated" / "derived_batch"))
    integrated = resolve_integrated_dir(integrated_cfg)
    qc_filtered = os.environ.get("TSAI_QC_FILTERED", str(tsai_processing / "01_QC_Filtered"))
    doublet_dir = os.environ.get("TSAI_DOUBLET_REMOVED", str(tsai_processing / "02_Doublet_Removed"))
    markers = os.environ.get(
        "TSAI_MARKERS_RDS",
        str(repo / "Processing" / "Tsai" / "Pipeline" / "Resources" / "Brain_Human_PFC_Markers_Mohammadi2020.rds"),
    )
    return {
        "input_h5ad": args.input_h5ad or integrated / "tsai_annotated.h5ad",
        "qc_summary": args.qc_summary_csv or Path(qc_filtered) / "qc_summary.csv",
        "doublet_summary": args.doublet_summary_csv or Path(doublet_dir) / "doublet_summary.csv",
        "cluster_rankings": (args.cluster_rankings_csv
                             or integrated / "cluster_annotation_rankings.csv"),
        "mohammadi_rds": args.mohammadi_rds or Path(markers),
    }


def log(msg: str) -> None:
    print(msg, flush=True)


# ---------- Section A: cohort-level per-donor QC ----------

def section_A(adata: ad.AnnData, tab_dir: Path, fig_dir: Path,
              qc_summary_path: Path, doublet_summary_path: Path) -> dict:
    log("\n=== Section A: cohort-level per-donor QC ===")
    obs = adata.obs
    grouped = obs.groupby("projid", observed=True)
    donor_qc = pd.DataFrame({
        "n_cells": grouped.size(),
        "median_total_counts": grouped["total_counts"].median(),
        "median_n_genes_by_counts": grouped["n_genes_by_counts"].median(),
        "median_pct_counts_mt": grouped["pct_counts_mt"].median(),
        "median_pct_counts_ribo": grouped["pct_counts_ribo"].median(),
        "median_scDblFinder_score": grouped["scDblFinder.score"].median(),
        "p95_scDblFinder_score": grouped["scDblFinder.score"].quantile(0.95),
        "frac_doublet_class": grouped["scDblFinder.class"].apply(
            lambda s: float((s.astype(str) == "doublet").mean())
        ),
    })

    if qc_summary_path.exists():
        qc_in = pd.read_csv(qc_summary_path)
        keep = [c for c in ["sample_id", "n_cells_before", "n_cells_after",
                            "pct_retained", "median_pct_mt"] if c in qc_in.columns]
        qc_in = qc_in[keep].copy()
        qc_in["sample_id"] = qc_in["sample_id"].astype(str)
        donor_qc = donor_qc.reset_index()
        donor_qc["projid"] = donor_qc["projid"].astype(str)
        donor_qc = donor_qc.merge(qc_in, left_on="projid", right_on="sample_id", how="left")
        donor_qc = donor_qc.drop(columns=["sample_id"]).set_index("projid")
        log(f"  merged qc_summary on {donor_qc['n_cells_before'].notna().sum()} donors")
    if doublet_summary_path.exists():
        dbl_in = pd.read_csv(doublet_summary_path)
        keep = [c for c in ["sample_id", "n_cells_before", "n_singlets",
                            "n_doublets", "pct_doublets"] if c in dbl_in.columns]
        dbl_in = dbl_in[keep].copy()
        dbl_in.columns = [c if c == "sample_id" else f"dbl_{c}" for c in dbl_in.columns]
        dbl_in["sample_id"] = dbl_in["sample_id"].astype(str)
        donor_qc = donor_qc.reset_index()
        donor_qc["projid"] = donor_qc["projid"].astype(str)
        donor_qc = donor_qc.merge(dbl_in, left_on="projid", right_on="sample_id", how="left")
        donor_qc = donor_qc.drop(columns=["sample_id"]).set_index("projid")
        log(f"  merged doublet_summary on {donor_qc['dbl_n_singlets'].notna().sum()} donors")

    donor_qc.to_csv(tab_dir / "donor_qc.csv")
    log(f"  wrote tables/donor_qc.csv  ({len(donor_qc)} donors)")

    # Outlier flagging via robust MAD.
    def mad_flag(s: pd.Series, n_mad: float = 3.0, direction: str = "high") -> pd.Series:
        med = s.median()
        mad = (s - med).abs().median() * 1.4826
        if mad == 0:
            return pd.Series(False, index=s.index)
        z = (s - med) / mad
        return z > n_mad if direction == "high" else z < -n_mad

    flags = pd.DataFrame({
        "high_mt": mad_flag(donor_qc["median_pct_counts_mt"], 3.0, "high"),
        "low_n_genes": mad_flag(donor_qc["median_n_genes_by_counts"], 3.0, "low"),
        "high_doublet": mad_flag(donor_qc["frac_doublet_class"], 3.0, "high"),
    })
    flags["n_flags"] = flags.sum(axis=1)
    outliers = flags[flags["n_flags"] > 0].sort_values("n_flags", ascending=False)
    outliers = outliers.join(donor_qc[[
        "median_pct_counts_mt", "median_n_genes_by_counts", "frac_doublet_class", "n_cells",
    ]])
    outliers.to_csv(tab_dir / "donor_outliers.csv")
    log(f"  wrote tables/donor_outliers.csv  ({len(outliers)} donors flagged on >=1 metric, "
        f"{(outliers['n_flags'] >= 2).sum()} on >=2)")

    # Violin/box panel.
    metrics = [
        ("median_pct_counts_mt", "median % counts MT"),
        ("median_n_genes_by_counts", "median n_genes / cell"),
        ("frac_doublet_class", "doublet fraction (scDblFinder)"),
        ("median_total_counts", "median total counts / cell"),
    ]
    fig, axes = plt.subplots(1, len(metrics), figsize=(4 * len(metrics), 5))
    for ax, (col, title) in zip(axes, metrics):
        vals = donor_qc[col].dropna()
        ax.boxplot(vals, vert=True, showfliers=False, widths=0.6)
        ax.scatter(np.full(len(vals), 1) + (np.random.default_rng(0).normal(0, 0.04, len(vals))),
                   vals, s=8, alpha=0.4, color="steelblue")
        med = vals.median()
        ax.axhline(med, color="grey", linestyle="--", linewidth=0.8)
        ax.set_title(title, fontsize=10)
        ax.set_xticks([])
    fig.suptitle("Per-donor QC distributions (each point = one donor)", fontsize=12)
    fig.tight_layout()
    fig.savefig(fig_dir / "cohort_qc_violins.pdf", bbox_inches="tight")
    fig.savefig(fig_dir / "cohort_qc_violins.png", bbox_inches="tight", dpi=180)
    plt.close(fig)
    log("  wrote figures/cohort_qc_violins.{pdf,png}")

    return {"donor_qc": donor_qc, "donor_outliers": outliers}


# ---------- Section B: cluster x donor composition ----------

def section_B(adata: ad.AnnData, tab_dir: Path, fig_dir: Path,
              donor_outliers: pd.DataFrame) -> dict:
    log("\n=== Section B: cluster x donor composition ===")
    obs = adata.obs
    cross = pd.crosstab(obs["leiden_res0_5"], obs["projid"])
    cross_frac = cross.div(cross.sum(axis=1), axis=0)
    cross.to_csv(tab_dir / "cluster_donor_counts.csv")
    cross_frac.to_csv(tab_dir / "cluster_donor_composition.csv")
    log(f"  wrote cluster_donor_{{counts,composition}}.csv  ({cross.shape[0]} clusters x "
        f"{cross.shape[1]} donors)")

    # Inverse Simpson per cluster.
    inv_simpson = 1.0 / (cross_frac ** 2).sum(axis=1).replace(0, np.nan)
    # Map cluster -> dominant cell_type.
    ct_per_cluster = obs.groupby("leiden_res0_5", observed=True)["cell_type"].agg(
        lambda s: s.value_counts().idxmax() if len(s) else None
    )
    skew = pd.DataFrame({
        "cell_type": ct_per_cluster,
        "n_cells": cross.sum(axis=1),
        "n_donors_contributing": (cross > 0).sum(axis=1),
        "inverse_simpson": inv_simpson,
        "top_donor_frac": cross_frac.max(axis=1),
        "top_donor": cross_frac.idxmax(axis=1),
    }).sort_values("inverse_simpson")
    skew.to_csv(tab_dir / "cluster_donor_skew.csv")
    log("  wrote cluster_donor_skew.csv")

    # Heatmap: clusters x donors. Order rows by cell_type alphabetical, then by
    # cluster id within type. Order columns by outlier-rank (worst first).
    skew_sorted = skew.sort_values(["cell_type", "inverse_simpson"])
    donor_order = (donor_outliers.sort_values("n_flags", ascending=False).index.tolist()
                   if len(donor_outliers) else [])
    other_donors = [d for d in cross_frac.columns if d not in donor_order]
    full_donor_order = donor_order + other_donors
    mat = cross_frac.loc[skew_sorted.index, full_donor_order].values
    fig, ax = plt.subplots(figsize=(min(0.05 * len(full_donor_order) + 6, 30),
                                    min(0.25 * len(skew_sorted) + 3, 25)))
    im = ax.imshow(mat, aspect="auto", cmap="viridis", interpolation="nearest", vmin=0, vmax=0.2)
    ax.set_yticks(np.arange(len(skew_sorted)))
    ax.set_yticklabels([f"{cl} ({ct})" for cl, ct in
                        zip(skew_sorted.index, skew_sorted["cell_type"])],
                       fontsize=7)
    ax.set_xticks([])
    ax.set_xlabel(f"{len(full_donor_order)} donors (worst-QC first)")
    ax.set_ylabel("Leiden res0.5 cluster (cell_type)")
    ax.set_title("Cluster x donor composition (row-normalized; capped at 0.2)")
    plt.colorbar(im, ax=ax, label="frac of cluster from donor")
    fig.tight_layout()
    fig.savefig(fig_dir / "cluster_donor_heatmap.pdf", bbox_inches="tight")
    fig.savefig(fig_dir / "cluster_donor_heatmap.png", bbox_inches="tight", dpi=180)
    plt.close(fig)
    log("  wrote cluster_donor_heatmap.{pdf,png}")

    return {"cluster_donor_counts": cross, "cluster_donor_skew": skew}


# ---------- Section C: annotation confidence ----------

def read_ora_estimate(h5ad_path: Path) -> pd.DataFrame:
    """Read obsm['ora_estimate'] as a (cells x signatures) DataFrame."""
    with h5py.File(h5ad_path, "r") as f:
        g = f["obsm/ora_estimate"]
        col_order = list(g.attrs.get("column-order", []))
        cols = [c.decode() if isinstance(c, bytes) else str(c) for c in col_order]
        data = {c: g[c][:] for c in cols}
    return pd.DataFrame(data)


def section_C(adata: ad.AnnData, cluster_rankings_csv: Path,
              tab_dir: Path, fig_dir: Path,
              ora_df: pd.DataFrame | None = None) -> dict:
    log("\n=== Section C: annotation confidence ===")
    out = {}

    if cluster_rankings_csv.exists():
        rk = pd.read_csv(cluster_rankings_csv)
        # decoupler.rank_sources_groups long-format:
        # columns = group, reference, names (= source), statistic, meanchange, pvals, pvals_adj
        if {"group", "names", "meanchange"}.issubset(rk.columns):
            score_col = "meanchange" if "meanchange" in rk.columns else "statistic"
            rk = rk.sort_values(["group", score_col], ascending=[True, False])
            top2 = rk.groupby("group").head(2).copy()
            top2["rk"] = top2.groupby("group").cumcount() + 1
            piv = top2.pivot(index="group", columns="rk", values=["names", score_col])
            gap = pd.DataFrame({
                "top1_source": piv[("names", 1)],
                "top1_score": piv[(score_col, 1)],
                "top2_source": piv[("names", 2)],
                "top2_score": piv[(score_col, 2)],
            })
            gap["gap"] = gap["top1_score"] - gap["top2_score"]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                gap["ratio"] = gap["top1_score"] / gap["top2_score"].replace(0, np.nan)
            ct_per_cluster = adata.obs.groupby(
                "leiden_res0_5", observed=True)["cell_type"].agg(
                lambda s: s.value_counts().idxmax() if len(s) else None
            )
            ct_per_cluster.index = ct_per_cluster.index.astype(str)
            gap.index = gap.index.astype(str)
            gap["assigned_cell_type"] = ct_per_cluster.reindex(gap.index).values
            gap = gap.sort_values("gap")
            gap.to_csv(tab_dir / "cluster_annotation_gap.csv")
            log(f"  wrote cluster_annotation_gap.csv ({len(gap)} clusters; score={score_col})")
            out["cluster_gap"] = gap
        else:
            log(f"  WARN: {cluster_rankings_csv.name} schema unexpected; columns={list(rk.columns)}")
    else:
        log(f"  WARN: {cluster_rankings_csv} not found; skipping C1")

    # Per-cell ORA gap by cell_type.
    if ora_df is None:
        ora_df = read_ora_estimate(adata.uns.get("_input_h5ad_path", None)
                                   or Path(adata.uns.get("_input", "")))
    # Normalize column names: ORA has Ex-L2|3, obs has Ex-L2/3. Build a map.
    ora_cols = {c: c.replace("|", "/") for c in ora_df.columns}
    ora_df = ora_df.rename(columns=ora_cols)

    log(f"  per-cell ORA matrix: {ora_df.shape}")
    sorted_vals = np.sort(ora_df.values, axis=1)
    top1 = sorted_vals[:, -1]
    top2 = sorted_vals[:, -2]
    gap_per_cell = top1 - top2
    df = pd.DataFrame({
        "cell_type": adata.obs["cell_type"].astype(str).values,
        "ora_gap": gap_per_cell,
        "ora_top1": top1,
        "ora_top2": top2,
    })
    quant = df.groupby("cell_type")["ora_gap"].quantile([0.05, 0.25, 0.5, 0.75, 0.95]).unstack()
    quant["n_cells"] = df.groupby("cell_type").size()
    quant.to_csv(tab_dir / "per_cell_ora_gap_quantiles.csv")
    log("  wrote per_cell_ora_gap_quantiles.csv")

    # Violin per cell_type.
    order = sorted(df["cell_type"].unique())
    fig, ax = plt.subplots(figsize=(max(8, 0.5 * len(order)), 5))
    data = [df.loc[df["cell_type"] == ct, "ora_gap"].values for ct in order]
    bp = ax.boxplot(data, labels=order, showfliers=False, widths=0.6)
    ax.set_ylabel("ORA score: top1 - top2 (per cell)")
    ax.set_title("Annotation confidence by cell_type")
    plt.setp(ax.get_xticklabels(), rotation=60, ha="right", fontsize=8)
    fig.tight_layout()
    fig.savefig(fig_dir / "per_cell_ora_gap_by_celltype.pdf", bbox_inches="tight")
    fig.savefig(fig_dir / "per_cell_ora_gap_by_celltype.png", bbox_inches="tight", dpi=180)
    plt.close(fig)
    log("  wrote per_cell_ora_gap_by_celltype.{pdf,png}")
    out["per_cell_gap"] = quant
    return out


# ---------- Section D: pipeline-panel marker QC ----------

def load_mohammadi_markers(rds_path: Path, csv_fallback: Path | None = None) -> dict:
    """Returns {source: [genes...]} for the Mohammadi 2020 PFC panel.

    Reads from a sibling CSV (`tables/mohammadi_markers.csv`, columns
    `source,gene`) if present. The CSV is produced once by an Rscript helper:

        Rscript -e 'lst <- readRDS("<rds>"); df <- do.call(rbind,
          lapply(names(lst), function(n) data.frame(source=n, gene=lst[[n]],
          stringsAsFactors=FALSE))); write.csv(df, "<csv>", row.names=FALSE)'

    Falls back to rpy2 if the CSV is absent and rpy2 happens to be installed.
    """
    if csv_fallback is None:
        csv_fallback = Path(__file__).resolve().parent / "tables" / "mohammadi_markers.csv"
    if csv_fallback.exists():
        log(f"  loading Mohammadi markers from {csv_fallback}")
        df = pd.read_csv(csv_fallback)
        return {src: grp["gene"].tolist() for src, grp in df.groupby("source", sort=False)}

    log(f"  CSV not found; trying rpy2 against {rds_path}")
    try:
        from rpy2.robjects import r  # noqa: F401
    except ImportError:
        raise SystemExit(
            f"Mohammadi marker CSV not found at {csv_fallback} and rpy2 is not "
            "available. Generate the CSV with:\n  Rscript -e 'lst <- readRDS(\""
            f"{rds_path}\"); df <- do.call(rbind, lapply(names(lst), function(n) "
            "data.frame(source=n, gene=lst[[n]], stringsAsFactors=FALSE))); "
            f"write.csv(df, \"{csv_fallback}\", row.names=FALSE)'"
        )
    obj = r["readRDS"](str(rds_path))
    names = list(r["names"](obj))
    out = {}
    for i, src in enumerate(names):
        vec = obj[i]
        out[str(src)] = [str(g) for g in vec]
    return out


def fraction_expressing_by_celltype(adata: ad.AnnData, genes: list[str]) -> pd.DataFrame:
    """Per-cell-type fraction of cells with raw expression > 0 for each gene.
    Uses adata.raw to access the full gene matrix."""
    raw = adata.raw
    var_names = pd.Index(raw.var_names)
    gene_idx = var_names.get_indexer(genes)
    present = np.array(genes)[gene_idx >= 0]
    valid = gene_idx[gene_idx >= 0]
    if len(valid) == 0:
        return pd.DataFrame(index=adata.obs["cell_type"].astype(str).unique())
    mat = raw.X[:, valid]
    nz = mat > 0
    ct = adata.obs["cell_type"].astype(str).values
    rows = []
    for celltype in pd.Index(ct).unique():
        mask = ct == celltype
        if not mask.any():
            continue
        if sp.issparse(nz):
            frac = np.asarray(nz[mask].mean(axis=0)).ravel()
        else:
            frac = nz[mask].mean(axis=0)
        rows.append(pd.Series(frac, index=present, name=celltype))
    return pd.DataFrame(rows)


def section_D(adata: ad.AnnData, mohammadi_rds: Path,
              tab_dir: Path, fig_dir: Path) -> dict:
    log("\n=== Section D: pipeline-panel (Mohammadi) marker QC ===")
    try:
        moh = load_mohammadi_markers(mohammadi_rds,
                                     csv_fallback=tab_dir / "mohammadi_markers.csv")
    except SystemExit as e:
        log(f"  ERROR: {e}")
        return {}
    summary = []
    for src, genes in moh.items():
        summary.append({"source": src, "n_genes": len(genes),
                        "n_in_raw": int(sum(1 for g in genes if g in adata.raw.var_names))})
    pd.DataFrame(summary).to_csv(tab_dir / "mohammadi_markers_summary.csv", index=False)
    # mohammadi_markers.csv is the canonical input; do not overwrite.
    log(f"  loaded {len(moh)} Mohammadi sources, "
        f"{sum(s['n_genes'] for s in summary)} marker entries total")

    # Per cell_type: fraction expressing the Mohammadi genes of its own
    # assigned signature.
    all_pipeline_genes = sorted({g for genes in moh.values() for g in genes
                                 if g in adata.raw.var_names})
    log(f"  computing fraction-expressing on {len(all_pipeline_genes)} Mohammadi genes "
        f"(present in raw)")
    frac_all = fraction_expressing_by_celltype(adata, all_pipeline_genes)
    frac_all.to_csv(tab_dir / "marker_fraction_expressing_mohammadi_all.csv")
    log("  wrote marker_fraction_expressing_mohammadi_all.csv")

    # Per cell_type vs its own Mohammadi source: mean fraction-expressing of
    # that source's marker genes.
    own_rows = []
    for ct in adata.obs["cell_type"].astype(str).unique():
        # The Mohammadi source name matches the cell_type label (with |->/).
        # Find the Mohammadi key that matches.
        candidates = [src for src in moh if src.replace("|", "/") == ct]
        if not candidates:
            own_rows.append({"cell_type": ct, "mohammadi_source": None,
                             "n_markers_in_raw": 0, "mean_frac_expressing": np.nan,
                             "median_frac_expressing": np.nan})
            continue
        src = candidates[0]
        genes_in_raw = [g for g in moh[src] if g in frac_all.columns]
        if not genes_in_raw:
            own_rows.append({"cell_type": ct, "mohammadi_source": src,
                             "n_markers_in_raw": 0, "mean_frac_expressing": np.nan,
                             "median_frac_expressing": np.nan})
            continue
        fr = frac_all.loc[ct, genes_in_raw]
        own_rows.append({
            "cell_type": ct, "mohammadi_source": src,
            "n_markers_in_raw": len(genes_in_raw),
            "mean_frac_expressing": float(fr.mean()),
            "median_frac_expressing": float(fr.median()),
            "frac_markers_above_30pct": float((fr > 0.30).mean()),
        })
    own_df = pd.DataFrame(own_rows).sort_values("mean_frac_expressing")
    own_df.to_csv(tab_dir / "mohammadi_own_source_fractions.csv", index=False)
    log("  wrote mohammadi_own_source_fractions.csv")

    # Panel comparison: own canonical-panel mean fraction vs own Mohammadi
    # mean fraction.
    canon_rows = []
    for ct, cat in EXPECTED_CATEGORY.items():
        canon_genes = [g for g in CANONICAL_MARKERS[cat] if g in adata.raw.var_names]
        moh_match = [src for src in moh if src.replace("|", "/") == ct]
        moh_genes = ([g for g in moh[moh_match[0]] if g in adata.raw.var_names]
                     if moh_match else [])
        canon_frac = (fraction_expressing_by_celltype(adata, canon_genes)
                      if canon_genes else None)
        moh_frac = fraction_expressing_by_celltype(adata, moh_genes) if moh_genes else None
        canon_rows.append({
            "cell_type": ct,
            "canonical_category": cat,
            "n_canonical_in_raw": len(canon_genes),
            "n_mohammadi_in_raw": len(moh_genes),
            "canonical_mean_frac": (float(canon_frac.loc[ct].mean())
                                    if canon_frac is not None and ct in canon_frac.index
                                    else np.nan),
            "mohammadi_mean_frac": (float(moh_frac.loc[ct].mean())
                                    if moh_frac is not None and ct in moh_frac.index
                                    else np.nan),
        })
    panel_df = pd.DataFrame(canon_rows)
    panel_df.to_csv(tab_dir / "panel_comparison.csv", index=False)
    log("  wrote panel_comparison.csv")
    return {"moh_markers": moh, "own_df": own_df, "panel_df": panel_df}


# ---------- Section E: per-donor drill-down on problematic cell types ----------

def section_E(adata: ad.AnnData, problem_celltypes: dict,
              donor_qc: pd.DataFrame, tab_dir: Path, fig_dir: Path) -> dict:
    log("\n=== Section E: per-donor drill-down ===")
    out = {}
    raw = adata.raw
    var_names = pd.Index(raw.var_names)

    drop_curves = []
    for ct, failing_genes in problem_celltypes.items():
        if ct not in adata.obs["cell_type"].astype(str).unique():
            log(f"  skip {ct}: not in cell_type categories")
            continue
        log(f"  -- {ct} --")
        cell_mask = (adata.obs["cell_type"].astype(str).values == ct)
        if cell_mask.sum() == 0:
            continue
        donor_arr = adata.obs["projid"].astype(str).values[cell_mask]
        gene_idx = var_names.get_indexer(failing_genes)
        present_genes = np.array(failing_genes)[gene_idx >= 0]
        valid = gene_idx[gene_idx >= 0]
        if len(valid) == 0:
            log(f"    no failing genes found in raw; skipping")
            continue
        mat = raw.X[cell_mask][:, valid]
        nz = (mat > 0)
        if sp.issparse(nz):
            nz = np.asarray(nz.todense())
        else:
            nz = np.asarray(nz)

        # E1: per-donor x failing-gene fraction-expressing.
        donor_df = pd.DataFrame(nz, index=donor_arr, columns=present_genes)
        donor_frac = donor_df.groupby(level=0).mean()
        donor_frac["n_cells"] = donor_df.groupby(level=0).size()
        donor_frac.index.name = "projid"
        donor_frac = donor_frac.sort_values(list(present_genes))  # worst-overall first
        donor_frac.to_csv(tab_dir / f"donor_problem_markers_{ct.replace('/', '_')}.csv")

        # Heatmap (donors x genes).
        plot_donors = donor_frac.index[donor_frac["n_cells"] >= 5][:80]
        if len(plot_donors) > 0:
            mat_plot = donor_frac.loc[plot_donors, list(present_genes)].values
            fig, ax = plt.subplots(figsize=(max(5, 0.8 * len(present_genes)),
                                            min(0.18 * len(plot_donors) + 2, 16)))
            im = ax.imshow(mat_plot, aspect="auto", cmap="magma", vmin=0, vmax=1)
            ax.set_yticks(np.arange(len(plot_donors)))
            ax.set_yticklabels(plot_donors, fontsize=6)
            ax.set_xticks(np.arange(len(present_genes)))
            ax.set_xticklabels(present_genes, rotation=60, ha="right", fontsize=8)
            ax.set_title(f"{ct}: per-donor fraction expressing failing markers")
            plt.colorbar(im, ax=ax, label="fraction expressing")
            fig.tight_layout()
            fig.savefig(fig_dir / f"donor_problem_markers_{ct.replace('/', '_')}_heatmap.pdf",
                        bbox_inches="tight")
            plt.close(fig)

        # E2: per-donor QC vs per-donor mean failing-marker fraction.
        donor_frac["mean_frac"] = donor_frac[list(present_genes)].mean(axis=1)
        merged = donor_frac[["mean_frac", "n_cells"]].join(
            donor_qc[["median_pct_counts_mt", "median_n_genes_by_counts",
                      "frac_doublet_class"]], how="left"
        )
        merged = merged[merged["n_cells"] >= 5]
        if len(merged) > 0:
            fig, axes = plt.subplots(1, 3, figsize=(15, 5))
            for ax, (xcol, xlab) in zip(axes, [
                ("median_pct_counts_mt", "donor median % MT"),
                ("median_n_genes_by_counts", "donor median n_genes"),
                ("frac_doublet_class", "donor doublet fraction"),
            ]):
                ax.scatter(merged[xcol], merged["mean_frac"],
                           s=20, alpha=0.6, color="steelblue")
                ax.set_xlabel(xlab)
                ax.set_ylabel(f"mean frac expressing failing markers in {ct}")
                # Spearman.
                v = merged[[xcol, "mean_frac"]].dropna()
                if len(v) >= 5:
                    rho = v.corr(method="spearman").iloc[0, 1]
                    ax.set_title(f"spearman={rho:.2f}", fontsize=9)
            fig.suptitle(f"{ct}: donor QC vs failing-marker fraction "
                         f"(donors with >=5 {ct} cells)")
            fig.tight_layout()
            fig.savefig(fig_dir / f"donor_qc_vs_problem_marker_{ct.replace('/', '_')}.pdf",
                        bbox_inches="tight")
            plt.close(fig)

        # E3: drop-one-donor sensitivity. Iteratively drop donors with the
        # LOWEST mean_frac (i.e. donors most pulling the cell type down),
        # weighted by their n_cells, and recompute cohort-wide mean fraction.
        sens_rows = [{"donors_dropped": 0,
                      "n_cells_used": int(donor_frac["n_cells"].sum()),
                      "cohort_mean_frac": float(
                          (donor_frac[list(present_genes)] *
                           donor_frac[["n_cells"]].values).sum().sum()
                          / (donor_frac["n_cells"].sum() * len(present_genes))
                      )}]
        ranked = donor_frac.sort_values("mean_frac")  # worst-first
        for n_drop in range(1, min(20, len(ranked)) + 1):
            keep = ranked.iloc[n_drop:]
            n_cells_used = int(keep["n_cells"].sum())
            if n_cells_used == 0:
                break
            cohort_mean = float(
                (keep[list(present_genes)] * keep[["n_cells"]].values).sum().sum()
                / (n_cells_used * len(present_genes))
            )
            sens_rows.append({"donors_dropped": n_drop, "n_cells_used": n_cells_used,
                              "cohort_mean_frac": cohort_mean})
        sens_df = pd.DataFrame(sens_rows)
        sens_df["cell_type"] = ct
        drop_curves.append(sens_df)

        fig, ax = plt.subplots(figsize=(7, 5))
        ax.plot(sens_df["donors_dropped"], sens_df["cohort_mean_frac"], marker="o")
        ax.set_xlabel("donors dropped (worst-mean-fraction first)")
        ax.set_ylabel("cohort mean frac expressing failing markers")
        ax.set_title(f"{ct}: drop-one-donor sensitivity")
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(fig_dir / f"drop_one_donor_sensitivity_{ct.replace('/', '_')}.pdf",
                    bbox_inches="tight")
        fig.savefig(fig_dir / f"drop_one_donor_sensitivity_{ct.replace('/', '_')}.png",
                    bbox_inches="tight", dpi=180)
        plt.close(fig)

    if drop_curves:
        all_sens = pd.concat(drop_curves, ignore_index=True)
        all_sens.to_csv(tab_dir / "drop_one_donor_sensitivity.csv", index=False)
        log("  wrote drop_one_donor_sensitivity.csv (combined)")
        out["drop_curves"] = all_sens
    return out


# ---------- Section F: doublet drill-down ----------

def section_F(adata: ad.AnnData, tab_dir: Path, fig_dir: Path) -> dict:
    log("\n=== Section F: doublet drill-down ===")
    obs = adata.obs.copy()
    obs["scDblFinder.score"] = pd.to_numeric(obs["scDblFinder.score"], errors="coerce")
    grouped = obs.groupby("cell_type", observed=True)["scDblFinder.score"]
    summary = grouped.agg(["median", lambda s: s.quantile(0.9), "mean", "count"])
    summary.columns = ["median", "p90", "mean", "n_cells"]
    summary["frac_doublet_class"] = obs.groupby(
        "cell_type", observed=True
    )["scDblFinder.class"].apply(lambda s: float((s.astype(str) == "doublet").mean()))
    summary = summary.sort_values("median", ascending=False)
    summary.to_csv(tab_dir / "doublet_score_by_celltype.csv")
    log("  wrote doublet_score_by_celltype.csv")

    order = summary.index.tolist()
    fig, ax = plt.subplots(figsize=(max(8, 0.55 * len(order)), 5))
    data = [obs.loc[obs["cell_type"].astype(str) == ct, "scDblFinder.score"].dropna().values
            for ct in order]
    ax.boxplot(data, labels=order, showfliers=False, widths=0.6)
    ax.set_ylabel("scDblFinder.score")
    ax.set_title("Doublet score distribution by cell_type (ordered by median)")
    plt.setp(ax.get_xticklabels(), rotation=60, ha="right", fontsize=8)
    fig.tight_layout()
    fig.savefig(fig_dir / "doublet_score_by_celltype.pdf", bbox_inches="tight")
    fig.savefig(fig_dir / "doublet_score_by_celltype.png", bbox_inches="tight", dpi=180)
    plt.close(fig)
    log("  wrote doublet_score_by_celltype.{pdf,png}")
    return {"doublet_by_celltype": summary}


# ---------- main ----------

def main() -> None:
    args = parse_args()
    paths = resolve_paths(args)
    out_dir = args.output_dir.resolve()
    fig_dir = out_dir / "figures"
    tab_dir = out_dir / "tables"
    fig_dir.mkdir(parents=True, exist_ok=True)
    tab_dir.mkdir(parents=True, exist_ok=True)

    log(f"Input H5AD: {paths['input_h5ad']}")
    log(f"Output dir: {out_dir}")
    if not Path(paths["input_h5ad"]).exists():
        raise SystemExit(f"ERROR: input H5AD not found: {paths['input_h5ad']}")

    log("Loading adata...")
    adata = ad.read_h5ad(paths["input_h5ad"])
    log(f"  n_obs={adata.n_obs:,}  n_vars={adata.n_vars:,}  raw_vars={adata.raw.n_vars:,}")
    log(f"  cell_type categories: {sorted(adata.obs['cell_type'].astype(str).unique())}")

    # Smoke mode: subsample, run A only.
    if args.smoke:
        log("[smoke] subsampling to 50k cells, running Section A only")
        sc.pp.subsample(adata, n_obs=50_000, random_state=args.seed)
        section_A(adata, tab_dir, fig_dir, Path(paths["qc_summary"]),
                  Path(paths["doublet_summary"]))
        log("[smoke] done")
        return

    # Section A (always).
    a_out = section_A(adata, tab_dir, fig_dir, Path(paths["qc_summary"]),
                      Path(paths["doublet_summary"]))
    donor_qc = a_out["donor_qc"]
    donor_outliers = a_out["donor_outliers"]

    # If --celltype specified, restrict E/F to that one and skip B/C/D.
    if args.celltype is not None:
        problem = {args.celltype: PROBLEM_CELLTYPES.get(args.celltype, [])}
        if not problem[args.celltype]:
            log(f"  WARN: no PROBLEM_CELLTYPES entry for {args.celltype}; using canonical own panel")
            cat = EXPECTED_CATEGORY.get(args.celltype)
            problem[args.celltype] = CANONICAL_MARKERS.get(cat, [])
        section_E(adata, problem, donor_qc, tab_dir, fig_dir)
        section_F(adata, tab_dir, fig_dir)
        return

    # Full pass: B, C, D, E (all problem types), F.
    section_B(adata, tab_dir, fig_dir, donor_outliers)
    # ORA dataframe needed for C.
    ora_df = read_ora_estimate(Path(paths["input_h5ad"]))
    section_C(adata, Path(paths["cluster_rankings"]), tab_dir, fig_dir, ora_df=ora_df)
    del ora_df
    gc.collect()
    section_D(adata, Path(paths["mohammadi_rds"]), tab_dir, fig_dir)
    section_E(adata, PROBLEM_CELLTYPES, donor_qc, tab_dir, fig_dir)
    section_F(adata, tab_dir, fig_dir)

    # Provenance.
    meta = {
        "input_h5ad": str(paths["input_h5ad"]),
        "output_dir": str(out_dir),
        "n_cells": int(adata.n_obs),
        "n_vars_raw": int(adata.raw.n_vars),
        "scanpy_version": sc.__version__,
        "anndata_version": ad.__version__,
        "timestamp": datetime.utcnow().isoformat() + "Z",
    }
    (tab_dir / "run_metadata.json").write_text(json.dumps(meta, indent=2))
    log(f"\nDone. Figures: {fig_dir}\n      Tables:  {tab_dir}")


if __name__ == "__main__":
    main()
