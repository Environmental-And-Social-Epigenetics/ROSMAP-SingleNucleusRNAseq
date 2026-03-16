#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gc
import os
from pathlib import Path

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import anndata as ad
import decoupler as dc
import harmonypy as hm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from rpy2.robjects import r
sc.settings.verbosity = 0
sc.set_figure_params(figsize=(6, 6), frameon=False)


def default_paths() -> dict[str, Path]:
    script_path = Path(__file__).resolve()
    repo_root = script_path.parents[3]
    workspace_root = repo_root.parent
    return {
        "repo_root": repo_root,
        "workspace_root": workspace_root,
        "input_dir": workspace_root / "Tsai_Data" / "Processing_Outputs" / "02_Doublet_Removed",
        "output_dir": workspace_root / "Tsai_Data" / "Processing_Outputs" / "03_Integrated",
        "metadata_csv": repo_root
        / "Preprocessing"
        / "Tsai"
        / "02_Cellranger_Counts"
        / "Tracking"
        / "patient_metadata.csv",
        "markers_rds": script_path.parent / "Resources" / "Brain_Human_PFC_Markers_Mohammadi2020.rds",
        "derived_batches_csv": script_path.parent / "Resources" / "derived_batches.csv",
    }


def parse_args() -> argparse.Namespace:
    paths = default_paths()
    parser = argparse.ArgumentParser(
        description="Stage 3 Tsai integration, clustering, and ORA-based annotation."
    )
    parser.add_argument("--input-dir", type=Path, default=paths["input_dir"])
    parser.add_argument("--output-dir", type=Path, default=paths["output_dir"])
    parser.add_argument("--metadata-csv", type=Path, default=paths["metadata_csv"])
    parser.add_argument("--markers-rds", type=Path, default=paths["markers_rds"])
    parser.add_argument(
        "--sample-ids",
        type=str,
        default="",
        help="Comma-separated sample IDs to integrate. Defaults to all available singlet files.",
    )
    parser.add_argument(
        "--annotation-cluster-key",
        type=str,
        default="leiden_res0_5",
        help="Cluster key used for ORA-based cell type annotation.",
    )
    parser.add_argument("--n-pcs", type=int, default=30)
    parser.add_argument("--n-neighbors", type=int, default=30)
    parser.add_argument("--n-hvgs", type=int, default=3000, help="Number of highly variable genes.")
    parser.add_argument("--neighbor-metric", type=str, default="cosine", help="Distance metric for neighbor graph.")
    parser.add_argument("--umap-min-dist", type=float, default=0.15, help="UMAP minimum distance.")
    parser.add_argument(
        "--harmony-batch-key",
        type=str,
        default="derived_batch",
        help="obs column for Harmony batch correction. "
        "Options: 'derived_batch' (flowcell-based, ~41 groups), "
        "'projid' (per-patient, 480 groups — old default), "
        "or any other obs column.",
    )
    parser.add_argument(
        "--harmony-theta",
        type=float,
        default=2.0,
        help="Harmony diversity penalty. Lower = less aggressive correction. Default: 2.0.",
    )
    parser.add_argument(
        "--derived-batches-csv",
        type=Path,
        default=paths["derived_batches_csv"],
        help="CSV mapping projid -> derived_batch (from derive_batches.py).",
    )
    parser.add_argument(
        "--skip-harmony",
        action="store_true",
        help="Skip Harmony entirely. Use PCA embedding directly for "
        "neighbors/clustering (useful as a no-correction baseline).",
    )
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def discover_samples(input_dir: Path, metadata_csv: Path) -> list[str]:
    available = sorted(
        path.name.replace("_singlets.h5ad", "")
        for path in input_dir.glob("*_singlets.h5ad")
        if path.is_file()
    )
    if not metadata_csv.exists():
        return available

    metadata = pd.read_csv(metadata_csv, dtype={"projid": str})
    ordered = [sample for sample in metadata["projid"].dropna().astype(str) if sample in set(available)]
    remaining = sorted(set(available) - set(ordered))
    return ordered + remaining


def parse_requested_samples(raw_sample_ids: str, available_samples: list[str]) -> list[str]:
    if not raw_sample_ids.strip():
        return available_samples

    requested = [sample.strip() for sample in raw_sample_ids.split(",") if sample.strip()]
    missing = [sample for sample in requested if sample not in set(available_samples)]
    if missing:
        missing_str = ", ".join(missing)
        raise ValueError(f"Requested sample IDs are not available as singlet inputs: {missing_str}")
    return requested


def load_metadata(metadata_csv: Path) -> pd.DataFrame:
    metadata = pd.read_csv(metadata_csv, dtype={"projid": str, "batch": str})
    metadata = metadata.drop_duplicates(subset="projid", keep="first").set_index("projid")
    return metadata


def ensure_qc_metrics(adata: ad.AnnData) -> ad.AnnData:
    required_metrics = {
        "pct_counts_mt",
        "log1p_total_counts",
        "log1p_n_genes_by_counts",
        "pct_counts_in_top_20_genes",
    }
    if required_metrics.issubset(set(adata.obs.columns)):
        return adata

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
    return adata


def load_markers(markers_rds: Path) -> pd.DataFrame:
    marker_list = r["readRDS"](str(markers_rds))
    records: list[dict[str, object]] = []
    marker_names = list(marker_list.names)

    for cell_type, genes in zip(marker_names, list(marker_list)):
        for gene in list(genes):
            records.append(
                {
                    "source": str(cell_type),
                    "target": str(gene),
                    "weight": 1.0,
                }
            )

    markers_df = pd.DataFrame(records).drop_duplicates()
    if markers_df.empty:
        raise ValueError(f"No marker genes were loaded from {markers_rds}")
    return markers_df


def attach_metadata_columns(adata: ad.AnnData, sample_id: str, metadata: pd.DataFrame) -> ad.AnnData:
    row = metadata.loc[sample_id]
    adata.obs["sample_id"] = sample_id
    adata.obs["projid"] = sample_id
    adata.obs["batch"] = str(row["batch"])

    for column, value in row.items():
        if pd.isna(value):
            continue
        adata.obs[column] = str(value)

    return adata


def load_adatas(
    sample_ids: list[str], args: argparse.Namespace, metadata: pd.DataFrame
) -> tuple[list[ad.AnnData], list[str]]:
    """Return (adatas, loaded_sample_ids) — skips samples missing from metadata."""
    adatas: list[ad.AnnData] = []
    loaded_ids: list[str] = []
    skipped: list[str] = []
    for sample_id in sample_ids:
        input_path = args.input_dir / f"{sample_id}_singlets.h5ad"
        if sample_id not in metadata.index:
            print(f"[WARN] {sample_id}: missing from metadata — skipping")
            skipped.append(sample_id)
            continue
        adata = ad.read_h5ad(input_path)
        adata.layers.clear()  # Drop redundant counts layer to halve memory
        if hasattr(adata.X, "astype"):
            adata.X = adata.X.astype(np.float32)
        adata.obs_names_make_unique()
        adata.var_names_make_unique()
        adata = attach_metadata_columns(adata, sample_id, metadata)
        adata = ensure_qc_metrics(adata)
        adatas.append(adata)
        loaded_ids.append(sample_id)
        print(f"[load] {sample_id}: {adata.n_obs} cells x {adata.n_vars} genes")
    if skipped:
        print(f"[WARN] Skipped {len(skipped)} samples missing from metadata: {', '.join(skipped)}")
    return adatas, loaded_ids


def run_annotation(
    adata: ad.AnnData,
    markers_df: pd.DataFrame,
    cluster_key: str,
    output_dir: Path,
) -> ad.AnnData:
    if "ora_estimate" in adata.obsm:
        del adata.obsm["ora_estimate"]

    dc.run_ora(adata, markers_df, source="source", target="target", use_raw=True)
    if "ora_estimate" not in adata.obsm:
        raise RuntimeError(
            "dc.run_ora() produced no 'ora_estimate' in obsm. "
            f"Check that marker genes overlap with adata.raw.var_names "
            f"({adata.raw.n_vars if adata.raw is not None else 'NO .raw'} genes)."
        )
    acts = dc.get_acts(adata, obsm_key="ora_estimate")

    acts_values = np.asarray(acts.X)
    finite_mask = np.isfinite(acts_values)
    if finite_mask.any():
        max_finite = np.nanmax(acts_values[finite_mask])
        acts_values[~finite_mask] = max_finite
        acts.X = acts_values

    ranked = dc.rank_sources_groups(
        acts,
        groupby=cluster_key,
        reference="rest",
        method="t-test_overestim_var",
    )
    ranked.to_csv(output_dir / "cluster_annotation_rankings.csv", index=False)

    top_three = ranked.groupby("group").head(3)
    top_three.to_csv(output_dir / "cluster_annotation_top3.csv", index=False)

    annotation_dict = ranked.groupby("group").head(1).set_index("group")["names"].to_dict()
    adata.obs["cell_type"] = (
        adata.obs[cluster_key].astype(str).map(annotation_dict).fillna("Unassigned")
    )
    adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")
    adata.uns["cell_type_annotation"] = {
        "cluster_key": cluster_key,
        "markers_ranking_csv": str(output_dir / "cluster_annotation_rankings.csv"),
    }
    return adata


def save_figures(
    adata: ad.AnnData,
    output_dir: Path,
    cluster_key: str,
    markers_df: pd.DataFrame,
) -> None:
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = str(figures_dir)

    # --- 1. PCA elbow plot ---
    if "pca" in adata.uns:
        variance_ratio = adata.uns["pca"]["variance_ratio"]
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.plot(range(1, len(variance_ratio) + 1), variance_ratio, "o-", markersize=3)
        ax.set_xlabel("Principal Component")
        ax.set_ylabel("Variance Ratio")
        ax.set_title("PCA Elbow Plot")
        ax.axvline(30, color="red", linestyle="--", alpha=0.5, label="n_pcs=30")
        ax.legend()
        fig.savefig(figures_dir / "pca_elbow.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

    # --- 2. Pre/post Harmony comparison ---
    # Color by the batch key that was used for correction (stored in uns).
    harmony_params = adata.uns.get("harmony_params", {})
    batch_color_key = harmony_params.get("batch_key", "batch")
    if batch_color_key == "SKIPPED" or batch_color_key not in adata.obs.columns:
        batch_color_key = "batch" if "batch" in adata.obs.columns else None

    if "X_umap_pca" in adata.obsm and batch_color_key is not None:
        batch_values = adata.obs[batch_color_key].astype(str)
        unique_batches = sorted(batch_values.unique())
        cmap = plt.cm.get_cmap("tab20", len(unique_batches))
        batch_colors = {b: cmap(i) for i, b in enumerate(unique_batches)}
        colors = [batch_colors[b] for b in batch_values]

        correction_label = harmony_params.get("batch_key", "Harmony")
        if correction_label == "SKIPPED":
            correction_label = "No Correction (PCA)"
        else:
            correction_label = f"After Harmony ({correction_label})"

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        ax1.scatter(
            adata.obsm["X_umap_pca"][:, 0], adata.obsm["X_umap_pca"][:, 1],
            c=colors, s=1, alpha=0.3, rasterized=True,
        )
        ax1.set_title("Before Harmony (PCA)")
        ax1.set_xlabel("UMAP1")
        ax1.set_ylabel("UMAP2")
        ax1.set_aspect("equal")

        ax2.scatter(
            adata.obsm["X_umap"][:, 0], adata.obsm["X_umap"][:, 1],
            c=colors, s=1, alpha=0.3, rasterized=True,
        )
        ax2.set_title(correction_label)
        ax2.set_xlabel("UMAP1")
        ax2.set_ylabel("UMAP2")
        ax2.set_aspect("equal")

        fig.suptitle(f"Batch Integration Comparison (colored by {batch_color_key})", fontsize=13)
        fig.tight_layout()
        fig.savefig(figures_dir / "harmony_comparison.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

    # --- 3. QC metrics on UMAP ---
    qc_cols = [c for c in ["total_counts", "n_genes_by_counts", "pct_counts_mt"] if c in adata.obs.columns]
    if qc_cols:
        sc.pl.umap(adata, color=qc_cols, show=False, save="_qc_metrics.png")

    # --- 4. Multi-resolution clustering ---
    res_keys = [k for k in ["leiden_res0_2", "leiden_res0_5", "leiden_res1"] if k in adata.obs.columns]
    if res_keys:
        sc.pl.umap(adata, color=res_keys, show=False, save="_multi_resolution.png")

    # --- 5 & 6. Existing integration + cell type UMAPs ---
    umap_color_keys = [k for k in [batch_color_key, cluster_key] if k and k in adata.obs.columns]
    if umap_color_keys:
        sc.pl.umap(adata, color=umap_color_keys, show=False, save="_tsai_integration.png")
    if "cell_type" in adata.obs.columns:
        sc.pl.umap(adata, color=["cell_type"], legend_loc="on data", show=False, save="_tsai_celltypes.png")

    # --- 7. Cluster-batch composition ---
    composition_key = batch_color_key or "batch"
    if composition_key in adata.obs.columns and cluster_key in adata.obs.columns:
        ct = pd.crosstab(adata.obs[cluster_key], adata.obs[composition_key], normalize="index")
        fig, ax = plt.subplots(figsize=(max(8, len(ct) * 0.5), 6))
        ct.plot(kind="bar", stacked=True, ax=ax, width=0.85)
        ax.set_xlabel("Cluster")
        ax.set_ylabel("Proportion")
        ax.set_title(f"Batch Composition per Cluster ({composition_key})")
        ax.legend(title=composition_key, bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=6)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
        fig.savefig(figures_dir / "cluster_batch_composition.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

    # --- 8. ORA enrichment heatmap ---
    if "ora_estimate" in adata.obsm:
        ora_scores = pd.DataFrame(
            adata.obsm["ora_estimate"],
            index=adata.obs_names,
            columns=adata.obsm["ora_estimate"].columns
            if hasattr(adata.obsm["ora_estimate"], "columns")
            else [f"CT{i}" for i in range(adata.obsm["ora_estimate"].shape[1])],
        )
        ora_scores[cluster_key] = adata.obs[cluster_key].values
        mean_scores = ora_scores.groupby(cluster_key).mean()
        fig, ax = plt.subplots(figsize=(max(8, len(mean_scores.columns) * 0.6), max(6, len(mean_scores) * 0.4)))
        im = ax.imshow(mean_scores.values, aspect="auto", cmap="RdBu_r")
        ax.set_xticks(range(len(mean_scores.columns)))
        ax.set_xticklabels(mean_scores.columns, rotation=45, ha="right", fontsize=8)
        ax.set_yticks(range(len(mean_scores.index)))
        ax.set_yticklabels(mean_scores.index, fontsize=8)
        ax.set_xlabel("Cell Type (marker set)")
        ax.set_ylabel("Cluster")
        ax.set_title("ORA Enrichment Scores")
        fig.colorbar(im, ax=ax, shrink=0.8, label="Mean ORA Score")
        fig.savefig(figures_dir / "ora_enrichment_heatmap.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

    # --- 9. Marker dot plot ---
    if adata.raw is not None and not markers_df.empty and "cell_type" in adata.obs.columns:
        raw_genes = set(adata.raw.var_names)
        top_markers = markers_df.groupby("source").head(5)
        top_markers = top_markers[top_markers["target"].isin(raw_genes)]
        if not top_markers.empty:
            marker_genes = top_markers["target"].unique().tolist()
            sc.pl.dotplot(
                adata, var_names=marker_genes, groupby="cell_type",
                use_raw=True, show=False, save="_markers.png",
            )

    # --- 10. Cell type proportions ---
    if "cell_type" in adata.obs.columns:
        counts = adata.obs["cell_type"].value_counts().sort_values(ascending=True)
        fig, ax = plt.subplots(figsize=(8, max(4, len(counts) * 0.35)))
        ax.barh(range(len(counts)), counts.values, color="#4292c6", edgecolor="none")
        ax.set_yticks(range(len(counts)))
        ax.set_yticklabels(counts.index, fontsize=9)
        ax.set_xlabel("Number of Cells")
        ax.set_title("Cell Type Proportions")
        for i, v in enumerate(counts.values):
            ax.text(v + counts.max() * 0.01, i, f"{v:,}", va="center", fontsize=8)
        fig.savefig(figures_dir / "cell_type_proportions.png", dpi=150, bbox_inches="tight")
        plt.close(fig)


def main() -> None:
    args = parse_args()
    args.input_dir = args.input_dir.resolve()
    args.output_dir = args.output_dir.resolve()
    args.metadata_csv = args.metadata_csv.resolve()
    args.markers_rds = args.markers_rds.resolve()
    args.derived_batches_csv = args.derived_batches_csv.resolve()

    integrated_path = args.output_dir / "tsai_integrated.h5ad"
    annotated_path = args.output_dir / "tsai_annotated.h5ad"
    if annotated_path.exists() and not args.overwrite:
        print(f"[skip] {annotated_path} already exists")
        return

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # ── Checkpoint: resume from integrated object if annotation failed ──
    if integrated_path.exists() and not annotated_path.exists() and not args.overwrite:
        print(f"[resume] Loading existing integrated object from {integrated_path}")
        combined = ad.read_h5ad(integrated_path)
    else:
        metadata = load_metadata(args.metadata_csv)
        available_samples = discover_samples(args.input_dir, args.metadata_csv)
        selected_samples = parse_requested_samples(args.sample_ids, available_samples)
        if not selected_samples:
            raise SystemExit("No singlet inputs were found to integrate.")

        adatas, loaded_sample_ids = load_adatas(selected_samples, args, metadata)
        if not adatas:
            raise SystemExit("No samples could be loaded (all missing from metadata?).")
        combined = ad.concat(
            adatas,
            join="inner",
            merge="same",
            label="concat_sample_id",
            keys=loaded_sample_ids,
            index_unique="-",
        )
        del adatas
        gc.collect()

        combined.obs["sample_id"] = combined.obs["sample_id"].astype(str)
        combined.obs["batch"] = combined.obs["batch"].astype(str).astype("category")

        # Merge derived batch assignments if using flowcell-based batches.
        if args.harmony_batch_key == "derived_batch":
            derived_csv = args.derived_batches_csv.resolve()
            if not derived_csv.exists():
                raise SystemExit(
                    f"Derived batches CSV not found: {derived_csv}\n"
                    "Run derive_batches.py first, or use --harmony-batch-key projid"
                )
            derived = pd.read_csv(derived_csv, dtype=str).set_index("projid")
            combined.obs["derived_batch"] = (
                combined.obs["projid"]
                .map(derived["derived_batch"])
                .fillna("UNKNOWN")
                .astype("category")
            )
            n_unknown = (combined.obs["derived_batch"] == "UNKNOWN").sum()
            if n_unknown > 0:
                print(f"[WARN] {n_unknown} cells have no derived_batch assignment")
            n_groups = combined.obs["derived_batch"].nunique()
            print(f"[batch] Using derived_batch ({n_groups} groups) for Harmony")
        elif not args.skip_harmony:
            print(f"[batch] Using '{args.harmony_batch_key}' for Harmony")

        print(f"[concat] {combined.n_obs} cells x {combined.n_vars} genes")

        combined.layers["counts"] = combined.X.copy()
        sc.pp.normalize_total(combined)
        sc.pp.log1p(combined)
        combined.X = combined.X.astype(np.float32)

        # Save normalized full-gene data to disk for later ORA annotation,
        # instead of keeping it in memory as .raw through PCA/Harmony.
        raw_path = args.output_dir / "_raw_normalized.h5ad"
        combined.write_h5ad(raw_path)

        sc.pp.highly_variable_genes(combined, flavor="seurat_v3", n_top_genes=args.n_hvgs, layer="counts")
        combined = combined[:, combined.var["highly_variable"]].copy()
        sc.tl.pca(combined, svd_solver="arpack", n_comps=args.n_pcs, use_highly_variable=True)
        if "counts" in combined.layers:
            del combined.layers["counts"]

        # Compute pre-Harmony UMAP for batch-effect comparison visualization.
        sc.pp.neighbors(combined, use_rep="X_pca", n_neighbors=args.n_neighbors,
                         n_pcs=args.n_pcs, metric=args.neighbor_metric, key_added="pca_neighbors")
        sc.tl.umap(combined, neighbors_key="pca_neighbors", min_dist=args.umap_min_dist, random_state=0)
        combined.obsm["X_umap_pca"] = combined.obsm["X_umap"].copy()
        del combined.obsp["pca_neighbors_distances"], combined.obsp["pca_neighbors_connectivities"]
        combined.uns.pop("pca_neighbors", None)

        if not args.skip_harmony:
            batch_key = args.harmony_batch_key
            if batch_key not in combined.obs.columns:
                raise SystemExit(
                    f"Harmony batch key '{batch_key}' not found in obs columns. "
                    f"Available: {', '.join(combined.obs.columns[:20])}"
                )
            print(f"[harmony] batch_key={batch_key}, theta={args.harmony_theta}")
            harmony_result = hm.run_harmony(
                combined.obsm["X_pca"],
                combined.obs,
                batch_key,
                theta=args.harmony_theta,
            )
            combined.obsm["X_harmony"] = harmony_result.Z_corr.T
            neighbor_rep = "X_harmony"
        else:
            print("[harmony] SKIPPED — using PCA embedding directly")
            neighbor_rep = "X_pca"

        # Store provenance for downstream reference.
        combined.uns["harmony_params"] = {
            "batch_key": "SKIPPED" if args.skip_harmony else args.harmony_batch_key,
            "theta": None if args.skip_harmony else args.harmony_theta,
            "neighbor_rep": neighbor_rep,
        }

        sc.pp.neighbors(
            combined,
            use_rep=neighbor_rep,
            n_neighbors=args.n_neighbors,
            n_pcs=args.n_pcs,
            metric=args.neighbor_metric,
        )
        sc.tl.leiden(combined, key_added="leiden_res0_2", resolution=0.2)
        sc.tl.leiden(combined, key_added="leiden_res0_5", resolution=0.5)
        sc.tl.leiden(combined, key_added="leiden_res1", resolution=1.0)
        sc.tl.umap(combined, min_dist=args.umap_min_dist, random_state=0)

        # Reload full-gene normalized data as .raw (needed for ORA annotation
        # and preserved in the integrated checkpoint for resume support).
        raw_adata = ad.read_h5ad(raw_path)
        combined.raw = raw_adata
        del raw_adata
        gc.collect()

        combined.write_h5ad(integrated_path)
        print(f"[write] integrated object -> {integrated_path}")
        raw_path.unlink(missing_ok=True)

    markers_df = load_markers(args.markers_rds)
    combined = run_annotation(
        combined,
        markers_df=markers_df,
        cluster_key=args.annotation_cluster_key,
        output_dir=args.output_dir,
    )
    save_figures(combined, args.output_dir, args.annotation_cluster_key, markers_df)

    # Sanitize obsm column names: '/' is HDF5's path separator and breaks write.
    for key in list(combined.obsm.keys()):
        elem = combined.obsm[key]
        if hasattr(elem, "columns") and elem.columns.str.contains("/").any():
            elem.columns = elem.columns.str.replace("/", "|", regex=False)
            combined.obsm[key] = elem

    combined.write_h5ad(annotated_path)
    print(f"[write] annotated object -> {annotated_path}")


if __name__ == "__main__":
    main()
