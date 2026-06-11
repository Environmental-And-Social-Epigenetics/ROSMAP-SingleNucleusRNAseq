#!/usr/bin/env python3
"""
Pseudobulk broad_Exc in Python (too large for R's dgCMatrix: >2^31 nonzeros).

Reads broad_Exc.h5ad in chunks via h5py, sums counts per patient_id × sex,
and saves pseudobulk count matrices + metadata as CSV files that the R
DESeq2 script can read directly.

Usage:
    python pseudobulk_broad_Exc.py --integration projid
    python pseudobulk_broad_Exc.py --integration batch
"""
import argparse
import os
import sys

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
sys.stdout.reconfigure(line_buffering=True)

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp


def read_obs_from_h5ad(h5ad_path):
    """Read obs DataFrame from h5ad without loading X."""
    with h5py.File(h5ad_path, "r") as f:
        obs = f["obs"]
        barcodes = obs["_index"][:].astype("U")

        cols = {}
        for key in obs.keys():
            if key == "_index":
                continue
            item = obs[key]
            if isinstance(item, h5py.Group) and "codes" in item:
                codes = item["codes"][:]
                cats = item["categories"][:].astype("U")
                vals = cats[codes]
                vals[codes < 0] = ""
                cols[key] = vals
            elif isinstance(item, h5py.Dataset):
                data = item[:]
                if data.dtype.kind in ("S", "O"):
                    data = data.astype("U")
                cols[key] = data

    return pd.DataFrame(cols, index=barcodes)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--integration", default="batch", choices=["batch", "projid"])
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    repo_dir = "/orcd/data/lhtsai/001/mabdel03/ROSMAP_Code/Transcriptomics"

    h5ad_path = os.path.join(script_dir, f"celltype_splits_{args.integration}", "broad_Exc.h5ad")
    pheno_path = os.path.join(repo_dir, "Data", "Phenotypes", "dataset_652_basic_12-23-2021.csv")

    if not os.path.exists(h5ad_path):
        print(f"ERROR: {h5ad_path} does not exist")
        sys.exit(1)

    # Read var names
    with h5py.File(h5ad_path, "r") as f:
        var_names = f["var"]["_index"][:].astype("U")
    n_genes = len(var_names)
    print(f"Genes: {n_genes:,}")

    # Read obs metadata (indexed by barcode)
    print("Reading obs metadata...")
    obs_df = read_obs_from_h5ad(h5ad_path)
    n_total = len(obs_df)
    print(f"Total cells: {n_total:,}")

    # Add positional index for h5py row access
    obs_df["_h5_row"] = np.arange(n_total)

    # Load phenotype data
    pheno = pd.read_csv(pheno_path)
    pheno = pheno.rename(columns={"projid": "patient_id"})
    pheno_ace = pheno[pheno["tot_adverse_exp"].notna()].copy()
    pheno_ace["patient_id"] = pheno_ace["patient_id"].astype(str)
    pheno_ace["ace_aggregate"] = (
        (pheno_ace["tot_adverse_exp"] - pheno_ace["tot_adverse_exp"].mean()) / pheno_ace["tot_adverse_exp"].std()
        + (pheno_ace["early_hh_ses"] - pheno_ace["early_hh_ses"].mean()) / pheno_ace["early_hh_ses"].std()
    )
    print(f"Patients with ACE data: {len(pheno_ace)}")

    # Filter cells to those with ACE phenotype data
    obs_df["patient_id"] = obs_df["projid"].astype(str)
    ace_patients = set(pheno_ace["patient_id"])
    obs_df = obs_df[obs_df["patient_id"].isin(ace_patients)]
    print(f"Cells with ACE phenotype: {len(obs_df):,}")

    # Get sex from phenotype data
    sex_lookup = pheno_ace.drop_duplicates("patient_id").set_index("patient_id")["msex"].to_dict()
    obs_df["msex_pheno"] = obs_df["patient_id"].map(sex_lookup)

    # For each sex, pseudobulk and save
    for sex_code in [0, 1]:
        sex_label = "Fem" if sex_code == 0 else "Male"
        sex_cells = obs_df[obs_df["msex_pheno"] == sex_code]

        if len(sex_cells) < 10:
            print(f"\n  SKIP {sex_label}: too few cells ({len(sex_cells)})")
            continue

        unique_patients = sorted(sex_cells["patient_id"].unique())
        n_patients = len(unique_patients)
        print(f"\n  Sex: {sex_label} - cells: {len(sex_cells):,}, patients: {n_patients}")

        if n_patients < 5:
            print(f"  SKIP {sex_label}: too few patients ({n_patients})")
            continue

        # Build mapping: for each cell, its h5ad row index and patient index
        patient_to_idx = {p: i for i, p in enumerate(unique_patients)}
        h5_rows = sex_cells["_h5_row"].values
        patient_indices = np.array([patient_to_idx[p] for p in sex_cells["patient_id"].values])

        # Sort by h5_row for sequential disk access
        sort_order = np.argsort(h5_rows)
        h5_rows = h5_rows[sort_order]
        patient_indices = patient_indices[sort_order]

        # Pseudobulk: sum counts per patient by reading sparse rows from h5ad
        pseudo_counts = np.zeros((n_genes, n_patients), dtype=np.int64)

        print(f"  Reading and pseudobulking {len(h5_rows):,} cells from h5ad...")
        with h5py.File(h5ad_path, "r") as f:
            x_data = f["X"]["data"]
            x_indices = f["X"]["indices"]
            x_indptr = f["X"]["indptr"][:]

            # Process in batches for progress reporting
            batch_size = 50000
            for i in range(0, len(h5_rows), batch_size):
                if i % (batch_size * 5) == 0:
                    print(f"    {i:,} / {len(h5_rows):,} cells processed")

                batch_rows = h5_rows[i:i + batch_size]
                batch_patients = patient_indices[i:i + batch_size]

                for j, (row, p_idx) in enumerate(zip(batch_rows, batch_patients)):
                    ptr_s = x_indptr[row]
                    ptr_e = x_indptr[row + 1]
                    if ptr_e <= ptr_s:
                        continue
                    col_idx = x_indices[ptr_s:ptr_e]
                    vals = x_data[ptr_s:ptr_e]
                    pseudo_counts[col_idx, p_idx] += vals

        print(f"  Pseudobulk matrix: {n_genes} genes x {n_patients} patients")

        # Save count matrix (genes x patients)
        for pheno_name in ["tot_adverse_exp", "early_hh_ses", "ace_aggregate"]:
            out_dir = os.path.join(script_dir, f"results_{args.integration}", pheno_name)
            os.makedirs(out_dir, exist_ok=True)

        out_base = os.path.join(script_dir, f"results_{args.integration}")
        counts_df = pd.DataFrame(pseudo_counts, index=var_names, columns=unique_patients)
        counts_path = os.path.join(out_base, f"pseudobulk_counts_Exc_{sex_label}.csv")
        counts_df.to_csv(counts_path)
        print(f"  Saved: {counts_path}")

        # Save patient metadata
        patient_meta = pheno_ace[pheno_ace["patient_id"].isin(unique_patients)].copy()
        patient_meta = patient_meta.drop_duplicates(subset="patient_id")
        patient_meta = patient_meta.set_index("patient_id").loc[unique_patients].reset_index()
        meta_path = os.path.join(out_base, f"pseudobulk_meta_Exc_{sex_label}.csv")
        patient_meta.to_csv(meta_path, index=False)
        print(f"  Saved: {meta_path}")

    print("\nDone! Run run_deseq_broad_Exc.R to complete the analysis.")


if __name__ == "__main__":
    main()
