#!/usr/bin/env python3
"""
Pseudobulk broad_Exc.h5ad in Python — too large for R's dgCMatrix (>2^31 nnz).

Reads the broad_Exc.h5ad in chunks via h5py, sums raw counts per patient_id,
then writes a count matrix CSV (genes × patients) plus a metadata CSV that
the small R companion script can feed directly into DESeq2.

Usage:
    python pseudobulk_Exc_python.py --integration library_id
    python pseudobulk_Exc_python.py --integration patient_id --phenotype msex
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
sys.stdout.reconfigure(line_buffering=True)

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--integration", required=True,
                   choices=["library_id", "patient_id", "pool_batch", "derived_batch", "sequencing_date"])
    p.add_argument("--input-dir", type=Path, required=True,
                   help="celltype_splits_<integration>/")
    p.add_argument("--output-dir", type=Path, required=True,
                   help="results_<integration>/Exc_pseudobulk/")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    h5ad_path = args.input_dir / "broad_Exc.h5ad"
    if not h5ad_path.exists():
        h5ad_path = args.input_dir / "Exc.h5ad"
    if not h5ad_path.exists():
        raise SystemExit(f"No broad_Exc.h5ad or Exc.h5ad in {args.input_dir}")

    print(f"Loading {h5ad_path} (backed=r) ...")
    adata = ad.read_h5ad(h5ad_path, backed="r")
    print(f"  cells={adata.n_obs}  genes={adata.n_vars}")

    # Make sure projid (==patient_id) is present
    if "projid" not in adata.obs.columns:
        if "patient_id" in adata.obs.columns:
            adata.obs["projid"] = adata.obs["patient_id"].astype(str)
        else:
            raise SystemExit("Neither projid nor patient_id in obs")

    projids = adata.obs["projid"].astype(str).values
    unique_pid = sorted(np.unique(projids))
    pid_to_idx = {p: i for i, p in enumerate(unique_pid)}
    print(f"  unique patients: {len(unique_pid)}")

    # Aggregate counts per patient in chunks of cells.
    # adata.X is a backed sparse matrix; iterate row-blocks.
    n_genes = adata.n_vars
    n_pat = len(unique_pid)
    pseudo = np.zeros((n_genes, n_pat), dtype=np.float64)

    chunk = 100_000
    n = adata.n_obs
    for start in range(0, n, chunk):
        end = min(start + chunk, n)
        X = adata.X[start:end]
        if sp.issparse(X):
            X = X.toarray()
        pid_block = projids[start:end]
        # Map each cell to its patient column
        pid_idx = np.array([pid_to_idx[p] for p in pid_block])
        # Aggregate by adding each cell's gene vector to its patient column
        np.add.at(pseudo.T, pid_idx, X)
        if (start // chunk) % 10 == 0:
            print(f"  pseudobulked {end}/{n} cells")

    pseudo = pseudo.astype(np.int64)
    var_names = adata.var_names.astype(str).tolist()
    counts_df = pd.DataFrame(pseudo, index=var_names, columns=unique_pid)
    counts_path = args.output_dir / "pseudobulk_counts_Exc.csv"
    counts_df.to_csv(counts_path)
    print(f"[write] {counts_path}  ({counts_df.shape[0]} genes x {counts_df.shape[1]} patients)")

    # Metadata: one row per patient, with first-cell values for covariates.
    obs_df = adata.obs.reset_index(drop=False).rename(columns={"index": "barcode"})
    meta = obs_df.groupby("projid", as_index=False).first()
    meta["projid"] = meta["projid"].astype(str)
    meta = meta.set_index("projid").loc[unique_pid].reset_index()
    meta_path = args.output_dir / "pseudobulk_meta_Exc.csv"
    meta.to_csv(meta_path, index=False)
    print(f"[write] {meta_path}")

    print("Done.")


if __name__ == "__main__":
    main()
