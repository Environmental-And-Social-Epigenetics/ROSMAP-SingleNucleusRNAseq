#!/usr/bin/env python3
"""
Per-arm SCENIC regulon-activity association (no GRN re-run).

The expensive SCENIC steps (GRNBoost -> cisTarget -> AUCell) are phenotype-
independent and already produced an `auc_matrix.csv` per cell type (micropool x
regulon, indexed by "<projid>_pool<N>"). This script loads that cached AUCell
matrix and runs an OLS of each regulon's activity against the ACE phenotype with
the covariate set of a specified male AD-model arm:

    AUCell ~ tot_adverse_exp + age_death + pmi [+ arm AD covars] [+ AD_binary:pheno]

Males-only (the cached AUCell already comes from a Male_<CT> SCENIC run).

Cohort-agnostic: the cached AUCell matrix is indexed by "<projid>_pool<N>" for
BOTH the Tsai and DeJager cohorts (scenic_analysis.py micropools by `patient_id`
which equals the ROSMAP projid in both cohorts), so the phenotype merge on
`projid` is identical regardless of cohort. A ``--cohort {tsai,dejager}`` flag
is accepted for symmetry with scenic_analysis.py and for log provenance; it does
not change the association method.

Usage:
    python scenic_associate.py \
        --cohort tsai \
        --auc-matrix .../Male_Mic/auc_matrix.csv \
        --pheno-csv  $ACE_SCORES_CSV \
        --arm MaleContAD --phenotype tot_adverse_exp \
        --cell-type Mic \
        --output-dir .../results_derived_batch_MaleContAD/tot_adverse_exp/Mic
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "_shared"))
import arm_covariates as armcov  # noqa: E402

MIN_OBS = 5


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Per-arm SCENIC regulon association over cached AUCell.")
    p.add_argument("--cohort", default="tsai", choices=["tsai", "dejager"],
                   help="Cohort label for provenance/logging (default: tsai). The "
                        "cached AUCell is projid-keyed for both cohorts, so this does "
                        "not change the association method.")
    p.add_argument("--auc-matrix", required=True, type=Path,
                   help="Cached AUCell matrix CSV (micropool x regulon, index <projid>_pool<N>).")
    p.add_argument("--pheno-csv", required=True, type=Path)
    p.add_argument("--arm", required=True, help="Male AD-model arm (e.g. MaleContAD).")
    p.add_argument("--phenotype", default="tot_adverse_exp")
    p.add_argument("--cell-type", required=True)
    p.add_argument("--output-dir", required=True, type=Path)
    return p.parse_args()


def parse_projid(index: pd.Index) -> pd.Series:
    """Micropool id '<projid>_pool<N>' -> projid (str)."""
    return pd.Series(index, index=index).str.replace(r"_pool\d+$", "", regex=True).astype(str)


def main() -> None:
    args = parse_args()
    spec = armcov.get_spec(args.arm)
    phenotype = args.phenotype
    args.output_dir.mkdir(parents=True, exist_ok=True)

    auc = pd.read_csv(args.auc_matrix, index_col=0)
    auc.index = auc.index.astype(str)
    projid = parse_projid(auc.index)

    pheno = pd.read_csv(args.pheno_csv)
    pheno["projid"] = pheno["projid"].astype(str)
    pheno = pheno.drop_duplicates(subset=["projid"]).set_index("projid")
    if spec["needs_ad_binary"]:
        armcov.add_ad_binary(pheno)

    # covariate columns needed for this arm
    cov_cols = ["age_death", "pmi", phenotype]
    for c in spec["ad_covars"]:
        cov_cols.append(c)
    if spec["interaction"] and "AD_binary" not in cov_cols:
        cov_cols.append("AD_binary")
    cov_cols = [c for c in dict.fromkeys(cov_cols)]  # dedup, keep order
    missing = [c for c in cov_cols if c not in pheno.columns]
    if missing:
        raise SystemExit(f"ERROR: phenotype CSV missing columns for arm {args.arm}: {missing}")

    # per-micropool covariate frame, aligned to AUCell rows via projid
    meta = pheno.reindex(projid.values)[cov_cols].reset_index(drop=True)
    meta.index = auc.index

    # z-score continuous covariates for stable regression (match scenic_analysis)
    for col in ["age_death", "pmi"]:
        v = pd.to_numeric(meta[col], errors="coerce")
        sd = v.std()
        meta[col] = (v - v.mean()) / sd if sd and sd > 0 else 0.0

    formula = armcov.ols_formula(args.arm, phenotype, response="AUCell")
    print(f"=== SCENIC associate: {args.cell_type} / arm={args.arm} / cohort={args.cohort} ===")
    print(f"  AUCell: {auc.shape[0]} micropools x {auc.shape[1]} regulons")
    print(f"  Formula: {formula}")

    results = []
    for regulon in auc.columns:
        data = meta.copy()
        data["AUCell"] = pd.to_numeric(auc[regulon].values, errors="coerce")
        data[phenotype] = pd.to_numeric(data[phenotype], errors="coerce")
        data = data.dropna()
        if len(data) < MIN_OBS:
            continue
        try:
            fit = smf.ols(formula, data=data).fit()
            # ACE main-effect term (exact phenotype name)
            if phenotype not in fit.params.index:
                continue
            results.append({
                "regulon": regulon,
                "coef": fit.params[phenotype],
                "stderr": fit.bse[phenotype],
                "pvalue": fit.pvalues[phenotype],
                "n_obs": int(fit.nobs),
            })
        except Exception as exc:  # noqa: BLE001
            print(f"  OLS failed for {regulon}: {exc}")
            continue

    if not results:
        out = pd.DataFrame(columns=["regulon", "coef", "stderr", "pvalue", "padj", "n_obs"])
    else:
        out = pd.DataFrame(results)
        _, padj, _, _ = multipletests(out["pvalue"], method="fdr_bh")
        out["padj"] = padj
        out = out.sort_values("pvalue")

    out_path = args.output_dir / "regression_results.csv"
    out.to_csv(out_path, index=False)
    n_sig = int((out["padj"] < 0.05).sum()) if "padj" in out.columns and len(out) else 0
    print(f"  Wrote {out_path}  (regulons: {len(out)}, padj<0.05: {n_sig})")


if __name__ == "__main__":
    main()
