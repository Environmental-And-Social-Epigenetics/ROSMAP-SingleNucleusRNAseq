#!/usr/bin/env python3
"""
Prepare per-cell-type raw-count H5AD files for Tsai Resilient pseudobulk DEG.

This script reads the integrated annotated object for cell-type labels and the
doublet-removed per-sample H5AD files for raw counts, then writes one H5AD per
cell type plus broad excitatory/inhibitory groups into ANALYSIS_OUTPUT_ROOT.
"""

from __future__ import annotations

import argparse
import gc
import os
import sys
from pathlib import Path

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
sys.stdout.reconfigure(line_buffering=True)

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

BATCH_SIZE = 50
REQUIRED_PHENO_COLUMNS = {
    "projid",
    "cogdx",
    "braaksc",
    "ceradsc",
    "msex",
    "age_death",
    "pmi",
    "niareagansc",
}


def repo_root() -> Path:
    return Path(__file__).resolve().parents[4]


def workspace_root() -> Path:
    return repo_root().parent


def env_path(name: str, fallback: Path) -> Path:
    return Path(os.environ.get(name, str(fallback))).expanduser()


def clean_ct_name(ct: str) -> str:
    return ct.replace("/", "_").replace(" ", "_").replace("(", "").replace(")", "")


def broad_group(ct_clean: str) -> str:
    if ct_clean.startswith("Ex-"):
        return "Exc"
    if ct_clean.startswith("In-"):
        return "Inh"
    return ct_clean


def normalize_integration(integration: str) -> str:
    aliases = {
        "batch": "derived_batch",
        "derived_batch": "derived_batch",
        "projid": "projid",
    }
    if integration not in aliases:
        choices = ", ".join(sorted(aliases))
        raise SystemExit(f"Unsupported integration '{integration}'. Expected one of: {choices}")
    return aliases[integration]


def default_paths(integration: str) -> tuple[Path, Path, Path]:
    integration = normalize_integration(integration)
    annotated = (
        env_path("TSAI_INTEGRATED_PROJID", workspace_root() / "Tsai_Data" / "Processing_Outputs" / "03_Integrated_projid")
        if integration == "projid"
        else env_path("TSAI_INTEGRATED", workspace_root() / "Tsai_Data" / "Processing_Outputs" / "03_Integrated")
    ) / "tsai_annotated.h5ad"
    doublet_dir = env_path(
        "TSAI_DOUBLET_REMOVED",
        workspace_root() / "Tsai_Data" / "Processing_Outputs" / "02_Doublet_Removed",
    )
    output_dir = (
        env_path("ANALYSIS_OUTPUT_ROOT", workspace_root() / "Analysis_Outputs")
        / "Resilient"
        / "DEG"
        / "Tsai"
        / f"celltype_splits_{integration}"
    )
    return annotated, doublet_dir, output_dir


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare per-cell-type Tsai Resilient H5AD files.")
    parser.add_argument("--integration", default="derived_batch", choices=["batch", "derived_batch", "projid"])
    parser.add_argument("--annotated-h5ad", type=Path, default=None)
    parser.add_argument("--doublet-dir", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument(
        "--celltype-filter",
        type=str,
        default="",
        help="Comma-separated cleaned cell type names to build (e.g. broad_Exc,Ast).",
    )
    parser.add_argument("--sample-limit", type=int, default=0, help="Limit the number of doublet-removed inputs.")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def read_obs_h5py(h5ad_path: Path) -> pd.DataFrame:
    with h5py.File(h5ad_path, "r") as f:
        obs = f["obs"]
        raw_idx = obs["_index"][:]
        barcodes = pd.array(raw_idx.astype("U"), dtype="string")

        def read_cat(name: str) -> pd.Categorical:
            grp = obs[name]
            codes = grp["codes"][:]
            raw_cats = grp["categories"][:]
            cats = [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in raw_cats]
            return pd.Categorical.from_codes(codes, categories=cats)

        required_obs = {"cell_type", "sample_id"}
        missing = sorted(required_obs - set(obs.keys()))
        if missing:
            raise SystemExit(f"Annotated H5AD is missing required obs fields: {', '.join(missing)}")

        cell_types = read_cat("cell_type")
        sample_ids = read_cat("sample_id")

    return pd.DataFrame(
        {
            "anno_barcode": barcodes,
            "cell_type": cell_types,
            "sample_id": sample_ids,
        }
    )


def validate_resilient_phenotypes() -> None:
    pheno_csv = env_path("RESILIENT_PHENOTYPE_CSV", repo_root() / "Data" / "Phenotypes" / "ROSMAP_clinical.csv")
    if not pheno_csv.exists():
        raise SystemExit(f"Resilient phenotype CSV not found: {pheno_csv}")
    df = pd.read_csv(pheno_csv, nrows=5)
    missing = sorted(REQUIRED_PHENO_COLUMNS - set(df.columns))
    if missing:
        raise SystemExit(f"Resilient phenotype CSV is missing required columns: {', '.join(missing)}")


def batched_concat(chunk_files: list[Path], batch_size: int = BATCH_SIZE) -> ad.AnnData:
    combined: ad.AnnData | None = None
    for i in range(0, len(chunk_files), batch_size):
        batch = chunk_files[i:i + batch_size]
        adatas = [ad.read_h5ad(cf) for cf in batch]
        batch_result = ad.concat(adatas, join="inner")
        del adatas
        gc.collect()
        if combined is None:
            combined = batch_result
        else:
            combined = ad.concat([combined, batch_result], join="inner")
            del batch_result
            gc.collect()
        print(f"    Batch {i // batch_size + 1}: {combined.n_obs} cells accumulated")
    assert combined is not None
    return combined


def selected_filter_values(raw_filter: str) -> set[str]:
    return {item.strip() for item in raw_filter.split(",") if item.strip()}


def should_write(name: str, filter_values: set[str]) -> bool:
    return not filter_values or name in filter_values


def main() -> None:
    args = parse_args()
    validate_resilient_phenotypes()
    args.integration = normalize_integration(args.integration)

    annotated_default, doublet_default, output_default = default_paths(args.integration)
    anno_path = (args.annotated_h5ad or annotated_default).resolve()
    dr_dir = (args.doublet_dir or doublet_default).resolve()
    out_dir = (args.output_dir or output_default).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir = out_dir / "_chunks"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    filter_values = selected_filter_values(args.celltype_filter)

    all_cell_types: set[str] = set()

    existing_chunks = sorted(tmp_dir.glob("*.h5ad"))
    if existing_chunks and not args.overwrite:
        print(f"Found {len(existing_chunks)} existing chunk files, skipping Pass 1.")
        for cf in existing_chunks:
            all_cell_types.add(cf.name.rsplit("__", 1)[0])
    else:
        print(f"Reading obs from {anno_path} via h5py...")
        obs_df = read_obs_h5py(anno_path)
        print(f"  {len(obs_df)} cells with annotations")

        obs_df["broad_type"] = obs_df["cell_type"].apply(lambda ct: broad_group(clean_ct_name(str(ct))))
        obs_df["dr_barcode"] = obs_df["anno_barcode"].str.rsplit("-", n=1).str[0]
        obs_df = obs_df.set_index("dr_barcode")
        obs_by_sample = {sid: grp for sid, grp in obs_df.groupby("sample_id")}
        sample_ids = sorted(obs_by_sample.keys())
        if args.sample_limit > 0:
            sample_ids = sample_ids[: args.sample_limit]
        print(f"  {len(sample_ids)} samples to process")

        del obs_df
        gc.collect()

        print("\n--- Pass 1: Splitting samples into per-celltype chunks ---")
        for i, sample_id in enumerate(sample_ids):
            dr_path = dr_dir / f"{sample_id}_singlets.h5ad"
            if not dr_path.exists():
                print(f"  [skip] {sample_id}: missing {dr_path}")
                continue

            adata_dr = ad.read_h5ad(dr_path)
            anno_grp = obs_by_sample.get(sample_id)
            if anno_grp is None:
                del adata_dr
                continue

            matched_mask = adata_dr.obs_names.isin(anno_grp.index)
            n_match = int(matched_mask.sum())
            if n_match == 0:
                del adata_dr
                continue

            if (i + 1) % 25 == 0 or i == 0:
                print(f"  [{i + 1}/{len(sample_ids)}] {sample_id}: {n_match}/{adata_dr.n_obs} annotated cells")

            matched_bcs = adata_dr.obs_names[matched_mask]
            adata_sub = adata_dr[matched_bcs].copy()
            del adata_dr

            anno_matched = anno_grp.loc[matched_bcs.tolist()]
            adata_sub.obs["cell_type"] = anno_matched["cell_type"].astype(str).values
            adata_sub.obs["broad_type"] = anno_matched["broad_type"].astype(str).values
            adata_sub.obs["projid"] = sample_id
            adata_sub.obs["sample_id"] = sample_id
            adata_sub.obs_names = pd.Index(anno_matched["anno_barcode"].astype(str).values)

            if sp.issparse(adata_sub.X):
                adata_sub.X = adata_sub.X.astype(np.int32)
            else:
                adata_sub.X = np.asarray(adata_sub.X, dtype=np.int32)

            for ct in adata_sub.obs["cell_type"].unique():
                ct_clean = clean_ct_name(str(ct))
                all_cell_types.add(ct_clean)
                if not should_write(ct_clean, filter_values):
                    continue
                mask = np.array(adata_sub.obs["cell_type"] == ct)
                sub = adata_sub[mask].copy()
                sub.write_h5ad(tmp_dir / f"{ct_clean}__{sample_id}.h5ad")
                del sub

            del adata_sub
            gc.collect()

    print("\n--- Pass 2: Concatenating per-celltype files ---")
    for ct_clean in sorted(all_cell_types):
        if not should_write(ct_clean, filter_values):
            continue
        out_path = out_dir / f"{ct_clean}.h5ad"
        if out_path.exists() and not args.overwrite:
            print(f"  {ct_clean}: already exists, skipping")
            continue
        chunk_files = sorted(tmp_dir.glob(f"{ct_clean}__*.h5ad"))
        if not chunk_files:
            continue
        print(f"  {ct_clean}: concatenating {len(chunk_files)} chunks...")
        combined = batched_concat(chunk_files)
        combined.obs_names_make_unique()
        combined.write_h5ad(out_path)
        print(f"  {ct_clean}: {combined.n_obs} cells x {combined.n_vars} genes")
        del combined
        gc.collect()

    broad_map: dict[str, list[str]] = {}
    for ct_clean in all_cell_types:
        broad_map.setdefault(broad_group(ct_clean), []).append(ct_clean)

    for broad_name, subtypes in broad_map.items():
        broad_file = f"broad_{broad_name}"
        if not should_write(broad_file, filter_values):
            continue
        if len(subtypes) == 1 and subtypes[0] == broad_name:
            continue
        out_path = out_dir / f"{broad_file}.h5ad"
        if out_path.exists() and not args.overwrite:
            print(f"  {broad_file}: already exists, skipping")
            continue
        chunk_files: list[Path] = []
        for subtype in subtypes:
            chunk_files.extend(sorted(tmp_dir.glob(f"{subtype}__*.h5ad")))
        if not chunk_files:
            continue
        print(f"  {broad_file}: concatenating {len(chunk_files)} chunks...")
        combined = batched_concat(chunk_files)
        combined.obs_names_make_unique()
        combined.write_h5ad(out_path)
        print(f"  {broad_file}: {combined.n_obs} cells x {combined.n_vars} genes")
        del combined
        gc.collect()

    var_names_path = out_dir / "var_names.csv"
    if (not var_names_path.exists()) or args.overwrite:
        dr_files = sorted(dr_dir.glob("*_singlets.h5ad"))
        if dr_files:
            with h5py.File(dr_files[0], "r") as f:
                raw_vn = f["var"]["_index"][:]
            var_names = [x.decode("utf-8") if isinstance(x, bytes) else str(x) for x in raw_vn]
            pd.DataFrame({"gene": var_names}).to_csv(var_names_path, index=False)
            print(f"Saved {var_names_path}")

    print(f"\nDone. Output directory: {out_dir}")


if __name__ == "__main__":
    main()
