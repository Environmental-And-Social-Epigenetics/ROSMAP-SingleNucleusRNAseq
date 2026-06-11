#!/usr/bin/env python3
"""
Stage 3b: specificity-based cell-type annotation for the Tsai integrated object.

Rebuilds the cell-type labels around the existing Mohammadi marker lists, fixing
two flaws of the ORA top-1-per-cluster approach in 03_integration_annotation.py:

  1. ORA is biased by marker-set size (Mohammadi signatures span 73-453 genes).
  2. The signatures overlap heavily, so labels can rest on non-specific genes
     (e.g. microglia markers appearing among PV-interneuron top genes).

Method (three levers):
  - trim_markers(): keep only each cell type's most SPECIFIC Mohammadi genes
    (tau specificity + argmax-belonging + detection), then truncate all types to
    the same top-N. Equal N removes size bias; specificity removes overlap.
  - score(): AUCell (rank-based, size-robust; primary) + ULM (specificity-weighted
    cross-check) via decoupler 2.x (dc.mt.aucell / dc.mt.ulm).
  - assign_labels(): per-cell argmax -> per-cluster majority vote + confidence,
    preserving the 17 cluster-level label names for downstream consumers.

NON-DESTRUCTIVE: writes tsai_annotated_v2.h5ad and a stage3b/ side directory;
never modifies tsai_annotated.h5ad. Old labels are carried as cell_type_ora_v1.

decoupler 2.x API verified (2.1.5): dc.mt.aucell(data, net, tmin, raw, empty,
bsize, verbose) -> obsm['score_aucell']; dc.mt.ulm(...) -> obsm['score_ulm'] +
obsm['padj_ulm']. Net columns: source, target, weight. raw=True selects adata.raw.
"""

from __future__ import annotations

import argparse
import gc
import json
import os
import sys
from datetime import datetime
from pathlib import Path

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
sys.stdout.reconfigure(line_buffering=True)

import anndata as ad  # noqa: E402
import decoupler as dc  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scanpy as sc  # noqa: E402
import scipy.sparse as sp  # noqa: E402
from scipy import stats  # noqa: E402

# The 17 canonical label names (slash form, matching obs["cell_type"]).
CELLTYPE_ORDER = [
    "Ast", "Endo", "Ex-L2/3", "Ex-L4", "Ex-L4/5", "Ex-L5", "Ex-L5/6",
    "Ex-L5/6-CC", "Ex-NRGN", "In-PV (Basket)", "In-PV (Chandelier)",
    "In-Rosehip", "In-SST", "In-VIP", "Mic", "OPC", "Oli",
]


def log(msg: str) -> None:
    print(msg, flush=True)


# ---------- path resolution (mirrors donor_qc.resolve_integrated_dir) ----------

def resolve_integrated_dir(configured: str) -> Path:
    """Locate the dir containing tsai_annotated.h5ad. The repo convention is
    03_Integrated/<variant>/, but existing data uses the flat 03_Integrated/.
    Probe configured -> parent -> derived_batch child."""
    cfg = Path(configured)
    seen: set[Path] = set()
    for cand in (cfg, cfg.parent, cfg / "derived_batch"):
        if cand in seen:
            continue
        seen.add(cand)
        if (cand / "tsai_annotated.h5ad").exists():
            if cand != cfg:
                log(f"  [layout] integrated dir resolved to {cand} (configured: {cfg})")
            return cand
    return cfg


def default_paths() -> dict:
    repo_root = Path(__file__).resolve().parents[3]
    tsai_processing = repo_root / "Data" / "Transcriptomics" / "Tsai" / "Processing_Outputs"
    integrated_cfg = os.environ.get("TSAI_INTEGRATED", str(tsai_processing / "03_Integrated" / "derived_batch"))
    markers = os.environ.get(
        "TSAI_MARKERS_RDS",
        str(repo_root / "Processing" / "Tsai" / "Pipeline" / "Resources"
            / "Brain_Human_PFC_Markers_Mohammadi2020.rds"),
    )
    return {"repo_root": repo_root, "integrated_cfg": integrated_cfg, "markers_rds": markers}


# ---------- marker loading (reuse Stage 3 net shape) ----------

def load_markers(markers_rds: Path) -> pd.DataFrame:
    """Read Mohammadi RDS -> long source,target,weight (weight=1.0). Uses rpy2."""
    from rpy2.robjects import r
    marker_list = r["readRDS"](str(markers_rds))
    names = list(marker_list.names)
    records = []
    for cell_type, genes in zip(names, list(marker_list)):
        for gene in list(genes):
            records.append({"source": str(cell_type), "target": str(gene), "weight": 1.0})
    df = pd.DataFrame(records).drop_duplicates()
    if df.empty:
        raise ValueError(f"No markers loaded from {markers_rds}")
    return df


# ---------- Step 1: specificity trimming ----------

def cluster_mean_matrix(adata: ad.AnnData, genes: list[str], cluster_key: str) -> pd.DataFrame:
    """Mean log-norm expression per leiden cluster (label-independent), over `genes`.
    Returns clusters x genes DataFrame from adata.raw.

    Implementation: build a sparse (n_clusters x n_obs) indicator matrix and
    compute `means = (I @ X[:, valid]) / cluster_sizes` in one shot. This
    avoids fancy CSR row-slicing per cluster (which would materialize many
    transient copies of the column-sliced submatrix) and keeps memory bounded
    to the n_clusters x len(genes) result.
    """
    raw = adata.raw
    idx = pd.Index(raw.var_names).get_indexer(genes)
    present = [g for g, i in zip(genes, idx) if i >= 0]
    valid = idx[idx >= 0]

    clusters = adata.obs[cluster_key].astype(str).values
    cluster_ids, inverse = np.unique(clusters, return_inverse=True)
    n_clusters = len(cluster_ids)
    n_obs = adata.n_obs
    counts = np.bincount(inverse, minlength=n_clusters).astype(np.float64)

    # (n_clusters x n_obs) sparse indicator. CSR for fast right-multiplication.
    ind = sp.csr_matrix(
        (np.ones(n_obs, dtype=np.float32),
         (inverse.astype(np.int64), np.arange(n_obs, dtype=np.int64))),
        shape=(n_clusters, n_obs),
    )

    X = raw.X
    if sp.issparse(X):
        # ind @ X is (n_clusters x n_obs) @ (n_obs x n_var_raw) = small dense.
        # This avoids materializing a column slice of the 6.6B-nnz raw matrix.
        X_csr = X if isinstance(X, sp.csr_matrix) else X.tocsr()
        sums_full = (ind @ X_csr).toarray()       # (n_clusters x n_var_raw)
        sums = sums_full[:, valid]
    else:
        sums = (ind @ X)[:, valid]

    means = sums / counts[:, None]
    df = pd.DataFrame(means, index=cluster_ids, columns=present)
    return df


def fraction_expressing(adata: ad.AnnData, genes: list[str]) -> pd.Series:
    """Cohort-wide fraction of cells with raw expression > 0, per gene.
    Counts CSC column nonzeros directly to avoid materializing a boolean copy
    of a ~2.4M-row CSR submatrix."""
    raw = adata.raw
    idx = pd.Index(raw.var_names).get_indexer(genes)
    present = [g for g, i in zip(genes, idx) if i >= 0]
    valid = idx[idx >= 0]
    X = raw.X
    n_obs = adata.n_obs
    if sp.issparse(X):
        X_csr = X if isinstance(X, sp.csr_matrix) else X.tocsr()
        # CSR stores one entry per stored nonzero (no duplicate (row,col)),
        # and log-normalized expression is nonneg, so binning the column-index
        # array directly gives the number of cells with value > 0 per gene.
        counts_all = np.bincount(X_csr.indices, minlength=X_csr.shape[1])
        frac = counts_all[valid].astype(np.float64) / n_obs
    else:
        block = X[:, valid] > 0
        frac = np.asarray(block).mean(axis=0)
    return pd.Series(frac, index=present)


def compute_tau(type_means: pd.DataFrame) -> pd.Series:
    """tau specificity per gene (cols) across cell types (rows). type_means is
    cell_type x gene, non-negative. tau in [0,1]; 1 = perfectly specific."""
    X = type_means.to_numpy(dtype=float)            # types x genes
    n = X.shape[0]
    col_max = X.max(axis=0)                          # per gene
    tau = np.full(X.shape[1], np.nan)
    nz = col_max > 0
    Xhat = X[:, nz] / col_max[nz]
    tau[nz] = (1.0 - Xhat).sum(axis=0) / (n - 1)
    return pd.Series(tau, index=type_means.columns)


def trim_markers(
    adata: ad.AnnData,
    markers_df: pd.DataFrame,
    cluster_key: str,
    top_n: int,
    tau_min: float,
    min_detect: float,
    allow_tau_fallback: bool,
    tau_floor: float,
    tab_dir: Path,
) -> pd.DataFrame:
    """Build an equal-size, high-specificity marker panel from the Mohammadi sets.

    Returns the trimmed long net (source, target, weight=tau, plus diagnostic cols).
    Writes trimmed_markers.csv, trimmed_markers_summary.csv,
    marker_overlap_before_after.csv.
    """
    log("\n=== Step 1: marker specificity trimming ===")
    moh_genes = sorted(set(markers_df["target"]) & set(adata.raw.var_names))
    log(f"  Mohammadi genes present in raw: {len(moh_genes)} / "
        f"{markers_df['target'].nunique()} unique")

    # Per-cluster means -> map clusters to provisional type via current cell_type
    # majority -> aggregate to a 17 x gene matrix (label-independent signal).
    cl_means = cluster_mean_matrix(adata, moh_genes, cluster_key)            # clusters x genes
    cl_to_type = (adata.obs.groupby(cluster_key, observed=True)["cell_type"]
                  .agg(lambda s: s.astype(str).value_counts().idxmax()))
    cl_to_type.index = cl_to_type.index.astype(str)
    cl_means["__type__"] = cl_to_type.reindex(cl_means.index).values
    type_means = cl_means.groupby("__type__").mean()                        # types x genes
    type_means = type_means.reindex([t for t in CELLTYPE_ORDER if t in type_means.index])
    log(f"  provisional type-mean matrix: {type_means.shape[0]} types x "
        f"{type_means.shape[1]} genes")

    # Detection + tau.
    frac = fraction_expressing(adata, moh_genes)
    tau = compute_tau(type_means)
    argmax_type = type_means.idxmax(axis=0)        # per gene -> which type is highest

    # Assemble per (type, gene) candidate table restricted to that type's own
    # Mohammadi membership.
    membership = markers_df[["source", "target"]].drop_duplicates()
    membership = membership[membership["target"].isin(moh_genes)]
    cand = membership.rename(columns={"source": "cell_type", "target": "gene"}).copy()
    cand["tau"] = cand["gene"].map(tau)
    cand["frac_expressing"] = cand["gene"].map(frac)
    cand["argmax_type"] = cand["gene"].map(argmax_type)
    cand["mean_expr_in_type"] = [
        type_means.loc[ct, g] if (ct in type_means.index and g in type_means.columns) else np.nan
        for ct, g in zip(cand["cell_type"], cand["gene"])
    ]

    # Keep criteria: detected, belongs (argmax==type), tau>=tau_min.
    detected = cand["frac_expressing"] >= min_detect
    belongs = cand["argmax_type"] == cand["cell_type"]
    spec = cand["tau"] >= tau_min
    cand["keep"] = detected & belongs & spec
    cand["spec_score"] = cand["tau"] * cand["mean_expr_in_type"]

    # Optional fallback for under-filled types: relax tau to tau_floor (keep
    # detection + belonging) only for types with < top_n survivors.
    if allow_tau_fallback:
        counts = cand[cand["keep"]].groupby("cell_type").size()
        underfilled = [t for t in type_means.index if counts.get(t, 0) < top_n]
        if underfilled:
            relax = (cand["cell_type"].isin(underfilled) & detected & belongs
                     & (cand["tau"] >= tau_floor))
            cand.loc[relax, "keep"] = True
            log(f"  tau-fallback applied to under-filled types: {underfilled} "
                f"(floor={tau_floor})")

    kept = cand[cand["keep"]].copy()
    kept = kept.sort_values(["cell_type", "spec_score"], ascending=[True, False])
    kept["rank"] = kept.groupby("cell_type").cumcount() + 1
    trimmed = kept[kept["rank"] <= top_n].copy()

    # Diagnostics + outputs.
    trimmed_out = trimmed[["cell_type", "gene", "tau", "mean_expr_in_type",
                           "argmax_type", "frac_expressing", "rank"]].rename(
        columns={"cell_type": "source"})
    trimmed_out.to_csv(tab_dir / "trimmed_markers.csv", index=False)

    summary_rows = []
    for t in type_means.index:
        n_moh = int((membership["source"] == t).sum())
        sub = cand[cand["cell_type"] == t]
        n_det = int((sub["frac_expressing"] >= min_detect).sum())
        n_belong = int(((sub["argmax_type"] == t) & (sub["frac_expressing"] >= min_detect)).sum())
        n_kept = int((trimmed["cell_type"] == t).sum())
        flag = ("severe_low_panel" if n_kept < 5 else
                "low_panel_quality" if n_kept < top_n else "ok")
        summary_rows.append({"source": t, "n_mohammadi_in": n_moh, "n_detected": n_det,
                             "n_belong": n_belong, "n_kept": n_kept, "panel_quality_flag": flag})
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(tab_dir / "trimmed_markers_summary.csv", index=False)
    log("  per-type panel sizes (n_kept):")
    for _, r in summary.iterrows():
        log(f"    {r['source']:20s} kept={r['n_kept']:3d}  ({r['panel_quality_flag']})")

    # Overlap before vs after (Jaccard over all type pairs). "Before" = each
    # type's Mohammadi set clipped to its top-N by mean expr (comparable basis).
    before_sets = {}
    for t in type_means.index:
        tg = membership.loc[membership["source"] == t, "target"]
        tg = tg[tg.isin(type_means.columns)]
        order = type_means.loc[t, list(tg)].sort_values(ascending=False)
        before_sets[t] = set(order.index[:top_n])
    after_sets = {t: set(trimmed.loc[trimmed["cell_type"] == t, "gene"]) for t in type_means.index}

    def jac(a, b):
        u = len(a | b)
        return (len(a & b) / u) if u else 0.0

    types = list(type_means.index)
    ov_rows = []
    for i in range(len(types)):
        for j in range(i + 1, len(types)):
            a, b = types[i], types[j]
            jb = jac(before_sets[a], before_sets[b])
            ja = jac(after_sets[a], after_sets[b])
            ov_rows.append({"type_a": a, "type_b": b, "jaccard_before": jb,
                            "jaccard_after": ja, "delta": ja - jb})
    overlap = pd.DataFrame(ov_rows).sort_values("jaccard_before", ascending=False)
    overlap.to_csv(tab_dir / "marker_overlap_before_after.csv", index=False)
    log(f"  mean pairwise Jaccard: before={overlap['jaccard_before'].mean():.3f}  "
        f"after={overlap['jaccard_after'].mean():.3f}")

    # Return net with weight = tau (used by ULM); AUCell net derived as weight=1.
    net = trimmed[["cell_type", "gene", "tau"]].rename(
        columns={"cell_type": "source", "gene": "target", "tau": "weight"})
    return net


# ---------- Step 2: scoring ----------

def build_score_adata(adata: ad.AnnData, trimmed_net: pd.DataFrame) -> ad.AnnData:
    """Build a compact AnnData on a bounded gene universe for scoring.

    decoupler's dc.mt.* cast the *entire* input matrix to float up front (it does
    not chunk the read), so scoring directly on the full adata.raw (2.4M x 36,601,
    ~6.6e9 nnz -> ~53 GiB of int64 indices alone) OOMs even at 300G. AUCell ranks
    marker genes within each cell, and the per-cell argmax across the 17 signatures
    is invariant to the gene universe as long as it is shared across signatures.
    So we score on the HVG set (the conventional ranking space, ~680M nnz) UNION
    every trimmed marker gene, which bounds memory ~10x while preserving the
    assignment semantics. log-normalized values come from adata.raw.
    """
    raw = adata.raw
    raw_names = pd.Index(raw.var_names)
    hvg = (adata.var_names[adata.var["highly_variable"]].tolist()
           if "highly_variable" in adata.var.columns else list(adata.var_names))
    markers = trimmed_net["target"].unique().tolist()
    universe = [g for g in pd.Index(hvg).union(pd.Index(markers)) if g in set(raw_names)]
    idx = raw_names.get_indexer(universe)
    idx = idx[idx >= 0]
    sub = raw.X[:, idx]
    if sp.issparse(sub):
        sub = sub.tocsr().astype(np.float32)
    else:
        sub = np.asarray(sub, dtype=np.float32)
    score_ad = ad.AnnData(X=sub)
    score_ad.obs_names = list(adata.obs_names)
    score_ad.var_names = list(raw_names[idx])
    n_markers_in = sum(1 for g in markers if g in set(score_ad.var_names))
    log(f"  score universe: {score_ad.n_vars} genes "
        f"(HVG {len(hvg)} ∪ markers {len(markers)}); "
        f"markers present: {n_markers_in}/{len(markers)}")
    return score_ad


def _run_dc(method, score_ad: ad.AnnData, net: pd.DataFrame, tmin: int,
            bsize: int, expect_key: str) -> None:
    """Call a dc.mt.* method and assert it wrote the expected obsm key, with a
    clear diagnostic if it didn't (a bare KeyError otherwise hides which keys
    decoupler actually produced). empty=False so decoupler cannot silently drop
    all-zero observations and desync the obsm length from the AnnData."""
    mname = getattr(method, "name", None) or getattr(method, "__name__", repr(method))
    before = set(score_ad.obsm.keys())
    method(score_ad, net=net, tmin=tmin, raw=False, empty=False, bsize=bsize, verbose=True)
    new = set(score_ad.obsm.keys()) - before
    if expect_key not in score_ad.obsm:
        raise RuntimeError(
            f"decoupler '{mname}' did not write obsm['{expect_key}']. "
            f"New obsm keys: {sorted(new)}; all obsm keys: {sorted(score_ad.obsm.keys())}; "
            f"score_ad: {score_ad.n_obs} obs x {score_ad.n_vars} vars; "
            f"net sources={net['source'].nunique()} entries={len(net)}")


def score(adata: ad.AnnData, trimmed_net: pd.DataFrame, tmin: int, bsize: int) -> None:
    log("\n=== Step 2: scoring (AUCell primary, ULM weighted) ===")
    score_ad = build_score_adata(adata, trimmed_net)
    aucell_net = trimmed_net.assign(weight=1.0)
    log(f"  AUCell on {trimmed_net['source'].nunique()} sources, "
        f"{len(aucell_net)} marker entries (raw=False on score universe, single batch)")
    _run_dc(dc.mt.aucell, score_ad, aucell_net, tmin, bsize, "score_aucell")
    adata.obsm["score_aucell"] = score_ad.obsm["score_aucell"]
    log(f"  -> obsm['score_aucell'] {adata.obsm['score_aucell'].shape}")
    log("  ULM (weight=tau) ...")
    try:
        _run_dc(dc.mt.ulm, score_ad, trimmed_net, tmin, bsize, "score_ulm")
        adata.obsm["score_ulm"] = score_ad.obsm["score_ulm"]
        if "padj_ulm" in score_ad.obsm:
            adata.obsm["padj_ulm"] = score_ad.obsm["padj_ulm"]
        log(f"  -> obsm['score_ulm'] {adata.obsm['score_ulm'].shape}")
    except Exception as e:
        # ULM has been observed to fail at full scale with
        # "ValueError: `ps` must be within [0, 1]" inside decoupler's BH FDR
        # routine — a numerical / float-precision edge case in 2.1.5. AUCell
        # is the primary score for assignment; ULM is a cross-check only, so
        # we surface the error and continue rather than block the whole run.
        log(f"  WARN: ULM scoring failed ({type(e).__name__}: {e}); "
            f"continuing with AUCell only.")
    del score_ad
    gc.collect()


# ---------- Step 3: assignment ----------

def assign_labels(adata: ad.AnnData, cluster_key: str, score_key: str,
                  gap_min: float | None, majority_min: float) -> dict:
    log("\n=== Step 3: assignment (per-cell argmax -> cluster majority) ===")
    scores = adata.obsm[score_key]
    if not isinstance(scores, pd.DataFrame):
        raise RuntimeError(f"obsm['{score_key}'] is not a DataFrame")
    cols = list(scores.columns)
    arr = scores.to_numpy(dtype=float)
    order = np.argsort(-arr, axis=1)
    top1_idx = order[:, 0]
    top2_idx = order[:, 1]
    top1 = arr[np.arange(arr.shape[0]), top1_idx]
    top2 = arr[np.arange(arr.shape[0]), top2_idx]
    percell = np.array(cols)[top1_idx]
    percell_top2 = np.array(cols)[top2_idx]
    gap = top1 - top2

    adata.obs["cell_type_percell"] = pd.Categorical(percell, categories=cols)
    adata.obs["cell_type_percell_top2"] = pd.Categorical(percell_top2, categories=cols)
    adata.obs["cell_type_percell_gap"] = gap

    if gap_min is None:
        # per-type 5th-percentile of gap, then take the median across types as a
        # single global threshold (documented heuristic).
        gdf = pd.DataFrame({"t": percell, "gap": gap})
        gap_min = float(gdf.groupby("t")["gap"].quantile(0.05).median())
    log(f"  gap_min (low-confidence threshold) = {gap_min:.4f}")

    # Cluster majority vote.
    clusters = adata.obs[cluster_key].astype(str).values
    pc = pd.Series(percell, index=adata.obs_names)
    maj = {}
    maj_frac = {}
    for cl in pd.unique(clusters):
        vc = pc[clusters == cl].value_counts()
        maj[cl] = vc.idxmax()
        maj_frac[cl] = float(vc.iloc[0] / vc.sum())
    cell_type = pd.Series(clusters, index=adata.obs_names).map(maj)
    adata.obs["cell_type"] = pd.Categorical(
        cell_type.values, categories=[c for c in CELLTYPE_ORDER if c in set(maj.values())])
    adata.obs["cluster_majority_frac"] = pd.Series(clusters, index=adata.obs_names).map(maj_frac).values

    low_cell = gap < gap_min
    low_cluster = adata.obs["cluster_majority_frac"].to_numpy() < majority_min
    adata.obs["low_confidence"] = low_cell | low_cluster

    n_low = int(adata.obs["low_confidence"].sum())
    log(f"  low-confidence cells: {n_low:,} / {adata.n_obs:,} "
        f"({100*n_low/adata.n_obs:.1f}%)")
    log("  cluster -> majority cell_type:")
    for cl in sorted(maj, key=lambda c: (maj[c], c)):
        log(f"    cl {cl:>3s}: {maj[cl]:20s} (maj_frac={maj_frac[cl]:.2f})")

    return {"gap_min": gap_min, "cluster_majority": maj, "cluster_majority_frac": maj_frac}


def synth_rankings(adata: ad.AnnData, cluster_key: str, score_key: str,
                   tab_dir: Path) -> None:
    """Emit cluster_annotation_rankings.csv in donor_qc's long schema
    (group, reference, names, statistic, meanchange, pvals, pvals_adj)."""
    log("\n  synthesizing cluster_annotation_rankings.csv (donor_qc schema)")
    scores = adata.obsm[score_key]
    cols = list(scores.columns)
    arr = scores.to_numpy(dtype=float)
    clusters = adata.obs[cluster_key].astype(str).values
    rows = []
    for cl in pd.unique(clusters):
        in_mask = clusters == cl
        out_mask = ~in_mask
        for j, name in enumerate(cols):
            a = arr[in_mask, j]
            b = arr[out_mask, j]
            t, p = stats.ttest_ind(a, b, equal_var=False, nan_policy="omit")
            rows.append({"group": cl, "reference": "rest", "names": name,
                         "statistic": float(t), "meanchange": float(np.nanmean(a) - np.nanmean(b)),
                         "pvals": float(p)})
    rk = pd.DataFrame(rows)
    # BH adjust within each group.
    rk["pvals_adj"] = np.nan
    for cl, sub in rk.groupby("group"):
        p = sub["pvals"].to_numpy()
        ok = ~np.isnan(p)
        adj = np.full_like(p, np.nan, dtype=float)
        if ok.any():
            pv = p[ok]
            order = np.argsort(pv)
            ranked = pv[order]
            m = len(pv)
            bh = ranked * m / (np.arange(1, m + 1))
            bh = np.minimum.accumulate(bh[::-1])[::-1]
            out = np.empty(m); out[order] = np.clip(bh, 0, 1)
            adj[ok] = out
        rk.loc[sub.index, "pvals_adj"] = adj
    rk = rk.sort_values(["group", "meanchange"], ascending=[True, False])
    rk.to_csv(tab_dir / "cluster_annotation_rankings.csv", index=False)
    rk.groupby("group").head(3).to_csv(tab_dir / "cluster_annotation_top3.csv", index=False)
    log(f"    wrote rankings ({len(rk)} rows) + top3")


def alias_ora_estimate(adata: ad.AnnData, score_key: str) -> None:
    """Alias the primary score into obsm['ora_estimate'] with /->| sanitized
    columns so donor_qc Section C (read_ora_estimate) works unchanged."""
    df = adata.obsm[score_key].copy()
    df.columns = [c.replace("/", "|") for c in df.columns]
    adata.obsm["ora_estimate"] = df


def sanitize_obsm_columns(adata: ad.AnnData) -> None:
    """anndata's h5ad writer treats '/' in a DataFrame column name as an HDF5
    group path separator, so a column like 'Ex-L4/5' collides with 'Ex-L4'.
    Replace '/' with '|' in every obsm DataFrame's columns before writing
    (matches the convention donor_qc.read_ora_estimate already maps back)."""
    for key, val in list(adata.obsm.items()):
        if isinstance(val, pd.DataFrame) and any("/" in str(c) for c in val.columns):
            val = val.copy()
            val.columns = [str(c).replace("/", "|") for c in val.columns]
            adata.obsm[key] = val


# ---------- Step 4: validation ----------

def validation(adata: ad.AnnData, tab_dir: Path) -> None:
    log("\n=== Step 4: validation (A/B vs v1) ===")
    if "cell_type_ora_v1" in adata.obs:
        conf = pd.crosstab(adata.obs["cell_type_ora_v1"].astype(str),
                           adata.obs["cell_type"].astype(str))
        conf.to_csv(tab_dir / "confusion_old_vs_new.csv")
        rows = []
        for old in conf.index:
            r = conf.loc[old]
            total = int(r.sum())
            same = int(r.get(old, 0))
            dest = r.drop(labels=[old], errors="ignore")
            top_dest = dest.idxmax() if len(dest) and dest.max() > 0 else "-"
            rows.append({"old_cell_type": old, "n_cells": total, "n_changed": total - same,
                         "frac_changed": (total - same) / total if total else 0.0,
                         "top_destination": top_dest})
        pd.DataFrame(rows).to_csv(tab_dir / "label_change_summary.csv", index=False)
        log(f"  wrote confusion_old_vs_new.csv + label_change_summary.csv")


# ---------- main ----------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Stage 3b specificity annotation")
    p.add_argument("--input-h5ad", type=Path, default=None)
    p.add_argument("--output-h5ad", type=Path, default=None)
    p.add_argument("--markers-rds", type=Path, default=None)
    p.add_argument("--cluster-key", default="leiden_res0_5")
    p.add_argument("--top-n", type=int, default=25)
    p.add_argument("--tau-min", type=float, default=0.6)
    p.add_argument("--min-detect", type=float, default=0.10)
    p.add_argument("--allow-tau-fallback", action="store_true")
    p.add_argument("--tau-floor", type=float, default=0.4)
    p.add_argument("--gap-min", type=float, default=None)
    p.add_argument("--majority-min", type=float, default=0.5)
    p.add_argument("--tmin", type=int, default=5)
    p.add_argument("--bsize", type=int, default=250000)
    p.add_argument("--primary-score", default="score_aucell",
                   choices=["score_aucell", "score_ulm"])
    p.add_argument("--subsample", type=int, default=0, help="Smoke: subsample N cells.")
    p.add_argument("--trim-only", action="store_true", help="Run Step 1 only and exit.")
    p.add_argument("--seed", type=int, default=0)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    paths = default_paths()
    integrated_dir = resolve_integrated_dir(args.input_h5ad.parent if args.input_h5ad
                                            else paths["integrated_cfg"])
    in_h5ad = args.input_h5ad or (integrated_dir / "tsai_annotated.h5ad")
    out_h5ad = args.output_h5ad or (integrated_dir / "tsai_annotated_v2.h5ad")
    markers_rds = args.markers_rds or Path(paths["markers_rds"])
    stage_dir = integrated_dir / "stage3b"
    tab_dir = stage_dir / "tables"
    fig_dir = stage_dir / "figures"
    tab_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    log(f"Input H5AD : {in_h5ad}")
    log(f"Output H5AD: {out_h5ad}")
    log(f"Markers    : {markers_rds}")
    log(f"Cluster key: {args.cluster_key}  top_n={args.top_n}  tau_min={args.tau_min}")
    if not in_h5ad.exists():
        raise SystemExit(f"ERROR: input not found: {in_h5ad}")
    if out_h5ad.resolve() == in_h5ad.resolve():
        raise SystemExit("ERROR: refusing to overwrite the v1 annotated object.")

    log("Loading adata ...")
    adata = ad.read_h5ad(in_h5ad)
    log(f"  n_obs={adata.n_obs:,}  n_vars={adata.n_vars:,}  raw_vars={adata.raw.n_vars:,}")
    if args.subsample and args.subsample < adata.n_obs:
        log(f"  [smoke] subsampling to {args.subsample}")
        sc.pp.subsample(adata, n_obs=args.subsample, random_state=args.seed)

    markers_df = load_markers(markers_rds)
    trimmed_net = trim_markers(adata, markers_df, args.cluster_key, args.top_n,
                               args.tau_min, args.min_detect, args.allow_tau_fallback,
                               args.tau_floor, tab_dir)
    if args.trim_only:
        log("\n[trim-only] done.")
        return

    score(adata, trimmed_net, args.tmin, args.bsize)

    # Preserve old labels before overwriting.
    if "cell_type" in adata.obs:
        adata.obs["cell_type_ora_v1"] = adata.obs["cell_type"].astype(str).values

    info = assign_labels(adata, args.cluster_key, args.primary_score,
                         args.gap_min, args.majority_min)
    synth_rankings(adata, args.cluster_key, args.primary_score, tab_dir)
    alias_ora_estimate(adata, args.primary_score)
    validation(adata, tab_dir)

    meta = {
        "input_h5ad": str(in_h5ad), "output_h5ad": str(out_h5ad),
        "markers_rds": str(markers_rds), "cluster_key": args.cluster_key,
        "top_n": args.top_n, "tau_min": args.tau_min, "min_detect": args.min_detect,
        "primary_score": args.primary_score, "gap_min": info["gap_min"],
        "n_cells": int(adata.n_obs), "decoupler": dc.__version__,
        "timestamp": datetime.utcnow().isoformat() + "Z",
    }
    (tab_dir / "stage3b_run_metadata.json").write_text(json.dumps(meta, indent=2))

    # anndata >= 0.11 refuses to serialize nullable/Arrow string arrays (e.g.
    # the raw var _index inherited from the v1 object) unless opted in.
    try:
        ad.settings.allow_write_nullable_strings = True
    except Exception:
        pass

    # '/' in obsm DataFrame columns breaks the h5ad writer (treated as HDF5
    # path separator); sanitize to '|' (donor_qc maps it back).
    sanitize_obsm_columns(adata)

    log(f"\nWriting {out_h5ad} ...")
    adata.write_h5ad(out_h5ad)
    log("Done.")
    gc.collect()


if __name__ == "__main__":
    main()
