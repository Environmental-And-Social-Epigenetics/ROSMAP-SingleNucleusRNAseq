#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from pathlib import Path

import anndata as ad
import decoupler as dc
import harmonypy as hm
import numpy as np
import pandas as pd
import scanpy as sc
from rpy2.robjects import r


os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
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
    parser.add_argument("--n-top-genes", type=int, default=3000)
    parser.add_argument("--n-pcs", type=int, default=50)
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--min-dist", type=float, default=0.15)
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


def save_figures(adata: ad.AnnData, output_dir: Path, cluster_key: str) -> None:
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = str(figures_dir)

    sc.pl.umap(adata, color=["batch", cluster_key], show=False, save="_tsai_integration.png")
    sc.pl.umap(adata, color=["cell_type"], legend_loc="on data", show=False, save="_tsai_celltypes.png")


def main() -> None:
    args = parse_args()
    args.input_dir = args.input_dir.resolve()
    args.output_dir = args.output_dir.resolve()
    args.metadata_csv = args.metadata_csv.resolve()
    args.markers_rds = args.markers_rds.resolve()

    integrated_path = args.output_dir / "tsai_integrated.h5ad"
    annotated_path = args.output_dir / "tsai_annotated.h5ad"
    if annotated_path.exists() and not args.overwrite:
        print(f"[skip] {annotated_path} already exists")
        return

    args.output_dir.mkdir(parents=True, exist_ok=True)

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
    combined.obs["sample_id"] = combined.obs["sample_id"].astype(str)
    combined.obs["batch"] = combined.obs["batch"].astype(str).astype("category")
    combined.layers["counts"] = combined.X.copy()

    print(f"[concat] {combined.n_obs} cells x {combined.n_vars} genes")

    sc.pp.normalize_total(combined, target_sum=1e4)
    sc.pp.log1p(combined)
    combined.raw = combined.copy()

    sc.pp.highly_variable_genes(
        combined,
        flavor="seurat_v3",
        layer="counts",
        n_top_genes=args.n_top_genes,
    )
    combined = combined[:, combined.var["highly_variable"]].copy()
    sc.pp.scale(combined, max_value=10)
    sc.tl.pca(combined, svd_solver="arpack", n_comps=args.n_pcs)

    harmony_result = hm.run_harmony(combined.obsm["X_pca"], combined.obs, "batch")
    combined.obsm["X_harmony"] = harmony_result.Z_corr.T

    sc.pp.neighbors(
        combined,
        use_rep="X_harmony",
        n_neighbors=args.n_neighbors,
    )
    sc.tl.leiden(combined, key_added="leiden_res0_2", resolution=0.2)
    sc.tl.leiden(combined, key_added="leiden_res0_5", resolution=0.5)
    sc.tl.leiden(combined, key_added="leiden_res1", resolution=1.0)
    sc.tl.umap(combined, min_dist=args.min_dist)

    combined.write_h5ad(integrated_path)
    print(f"[write] integrated object -> {integrated_path}")

    markers_df = load_markers(args.markers_rds)
    combined = run_annotation(
        combined,
        markers_df=markers_df,
        cluster_key=args.annotation_cluster_key,
        output_dir=args.output_dir,
    )
    save_figures(combined, args.output_dir, args.annotation_cluster_key)
    combined.write_h5ad(annotated_path)
    print(f"[write] annotated object -> {annotated_path}")


if __name__ == "__main__":
    main()
