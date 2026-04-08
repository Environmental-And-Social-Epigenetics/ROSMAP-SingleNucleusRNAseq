#!/usr/bin/env python3
"""
Prepare per-cell-type raw-count H5AD files for DeJager ACE pseudobulk DEG.
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
    "tot_adverse_exp",
    "early_hh_ses",
    "msex",
    "age_death",
    "pmi",
    "niareagansc",
}
INTEGRATION_TO_ENV = {
    "library_id": "DEJAGER_INTEGRATED",
    "patient_id": "DEJAGER_INTEGRATED_PATIENT_ID",
    "pool_batch": "DEJAGER_INTEGRATED_POOL_BATCH",
    "derived_batch": "DEJAGER_INTEGRATED_DERIVED_BATCH",
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
    if ct_clean.startswith("Ex-") or ct_clean.startswith("Ex_"):
        return "Exc"
    if ct_clean.startswith("In-") or ct_clean.startswith("In_"):
        return "Inh"
    return ct_clean


def default_paths(integration: str) -> tuple[Path, Path, Path]:
    integrated_dir = env_path(
        INTEGRATION_TO_ENV[integration],
        workspace_root() / "DeJager_Data" / "Processing_Outputs" / f"03_Integrated_{integration}",
    )
    if integration == "library_id":
        integrated_dir = env_path(
            "DEJAGER_INTEGRATED",
            workspace_root() / "DeJager_Data" / "Processing_Outputs" / "03_Integrated",
        )
    annotated = integrated_dir / "dejager_annotated.h5ad"
    doublet_dir = env_path(
        "DEJAGER_DOUBLET_REMOVED",
        workspace_root() / "DeJager_Data" / "Processing_Outputs" / "02_Doublet_Removed",
    )
    output_dir = (
        env_path("ANALYSIS_OUTPUT_ROOT", workspace_root() / "Analysis_Outputs")
        / "ACE"
        / "DEG"
        / "DeJager"
        / f"celltype_splits_{integration}"
    )
    return annotated, doublet_dir, output_dir


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare per-cell-type DeJager ACE H5AD files.")
    parser.add_argument("--integration", default="library_id", choices=sorted(INTEGRATION_TO_ENV))
    parser.add_argument("--annotated-h5ad", type=Path, default=None)
    parser.add_argument("--doublet-dir", type=Path, default=None)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--celltype-filter", type=str, default="")
    parser.add_argument("--sample-limit", type=int, default=0)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def validate_ace_phenotypes() -> None:
    pheno_csv = env_path("ACE_SCORES_CSV", repo_root() / "Data" / "Phenotypes" / "TSAI_DEJAGER_all_patients_wACEscores.csv")
    if not pheno_csv.exists():
        raise SystemExit(f"ACE phenotype CSV not found: {pheno_csv}")
    df = pd.read_csv(pheno_csv, nrows=5)
    missing = sorted(REQUIRED_PHENO_COLUMNS - set(df.columns))
    if missing:
        raise SystemExit(f"ACE phenotype CSV is missing required columns: {', '.join(missing)}")


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

        required_obs = {"cell_type", "library_id", "patient_id"}
        missing = sorted(required_obs - set(obs.keys()))
        if missing:
            raise SystemExit(f"Annotated H5AD is missing required obs fields: {', '.join(missing)}")

        return pd.DataFrame(
            {
                "anno_barcode": barcodes,
                "cell_type": read_cat("cell_type"),
                "library_id": read_cat("library_id"),
                "patient_id": read_cat("patient_id"),
            }
        )


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
    assert combined is not None
    return combined


def selected_filter_values(raw_filter: str) -> set[str]:
    return {item.strip() for item in raw_filter.split(",") if item.strip()}


def should_write(name: str, filter_values: set[str]) -> bool:
    return not filter_values or name in filter_values


def main() -> None:
    args = parse_args()
    validate_ace_phenotypes()

    annotated_default, doublet_default, output_default = default_paths(args.integration)
    anno_path = (args.annotated_h5ad or annotated_default).resolve()
    dr_dir = (args.doublet_dir or doublet_default).resolve()
    out_dir = (args.output_dir or output_default).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir = out_dir / "_chunks"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    filter_values = selected_filter_values(args.celltype_filter)

    obs_df = read_obs_h5py(anno_path)
    obs_df["broad_type"] = obs_df["cell_type"].apply(lambda ct: broad_group(clean_ct_name(str(ct))))
    obs_df["dr_barcode"] = obs_df["anno_barcode"].str.rsplit("-", n=1).str[0]
    obs_df = obs_df.set_index("dr_barcode")
    obs_by_library = {lib: grp for lib, grp in obs_df.groupby("library_id")}
    library_ids = sorted(obs_by_library)
    if args.sample_limit > 0:
        library_ids = library_ids[: args.sample_limit]

    all_cell_types: set[str] = set()
    print(f"Processing {len(library_ids)} DeJager libraries")

    for library_id in library_ids:
        input_path = dr_dir / f"{library_id}_singlets.h5ad"
        if not input_path.exists():
            print(f"[skip] {library_id}: missing {input_path}")
            continue

        adata_dr = ad.read_h5ad(input_path)
        anno_grp = obs_by_library.get(library_id)
        if anno_grp is None:
            del adata_dr
            continue

        matched_mask = adata_dr.obs_names.isin(anno_grp.index)
        if int(matched_mask.sum()) == 0:
            del adata_dr
            continue

        matched_bcs = adata_dr.obs_names[matched_mask]
        adata_sub = adata_dr[matched_bcs].copy()
        del adata_dr

        anno_matched = anno_grp.loc[matched_bcs.tolist()]
        adata_sub.obs["cell_type"] = anno_matched["cell_type"].astype(str).values
        adata_sub.obs["broad_type"] = anno_matched["broad_type"].astype(str).values
        adata_sub.obs["projid"] = anno_matched["patient_id"].astype(str).values
        adata_sub.obs["patient_id"] = anno_matched["patient_id"].astype(str).values
        adata_sub.obs["library_id"] = library_id
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
            sub = adata_sub[np.array(adata_sub.obs["cell_type"] == ct)].copy()
            sub.write_h5ad(tmp_dir / f"{ct_clean}__{library_id}.h5ad")
            del sub

        del adata_sub
        gc.collect()

    for ct_clean in sorted(all_cell_types):
        if not should_write(ct_clean, filter_values):
            continue
        out_path = out_dir / f"{ct_clean}.h5ad"
        if out_path.exists() and not args.overwrite:
            continue
        chunk_files = sorted(tmp_dir.glob(f"{ct_clean}__*.h5ad"))
        if not chunk_files:
            continue
        combined = batched_concat(chunk_files)
        combined.obs_names_make_unique()
        combined.write_h5ad(out_path)
        del combined
        gc.collect()

    broad_map: dict[str, list[str]] = {}
    for ct_clean in all_cell_types:
        broad_map.setdefault(broad_group(ct_clean), []).append(ct_clean)

    for broad_name, subtypes in broad_map.items():
        broad_file = f"broad_{broad_name}"
        if not should_write(broad_file, filter_values):
            continue
        chunk_files: list[Path] = []
        for subtype in subtypes:
            chunk_files.extend(sorted(tmp_dir.glob(f"{subtype}__*.h5ad")))
        if not chunk_files:
            continue
        out_path = out_dir / f"{broad_file}.h5ad"
        if out_path.exists() and not args.overwrite:
            continue
        combined = batched_concat(chunk_files)
        combined.obs_names_make_unique()
        combined.write_h5ad(out_path)
        del combined
        gc.collect()

    print(f"Done. Output directory: {out_dir}")


if __name__ == "__main__":
    main()
