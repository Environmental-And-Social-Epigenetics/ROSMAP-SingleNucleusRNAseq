#!/usr/bin/env python3
"""
ACE Transcription Factor & Pathway Activity Analysis

Infers TF activity (DoRothEA, CollecTRI) and pathway activity (PROGENy)
per patient using decoupler, then tests association with ACE phenotypes
via OLS regression with covariates.

Designed to process pre-split per-cell-type h5ad files produced by the
DEG preprocessing step.  Each file is pseudobulked per patient, then
TF/pathway activities are estimated with decoupler's MLM method and
tested for association with the specified ACE phenotype, stratified by
sex.

Usage:
    python tf_activity_analysis.py \
        --integration derived_batch \
        --phenotype tot_adverse_exp \
        --input-dir ${ACE_OUTPUT_ROOT}/DEG/Tsai/celltype_splits_derived_batch \
        --pheno-csv ${ACE_SCORES_CSV} \
        --output-dir ${ACE_OUTPUT_ROOT}/TFActivity/Tsai/results_derived_batch/tot_adverse_exp \
        --run-progeny
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path
from typing import Optional

# Force line-buffered stdout so SLURM logs stream in real time
sys.stdout.reconfigure(line_buffering=True)

# Prevent HDF5 file-locking issues on networked filesystems
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

import anndata as ad
import decoupler as dc
import numpy as np
import pandas as pd
import scanpy as sc
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

# Shared ACE phenotype loader (handles pooled-cohort duplicate projids) and the
# per-arm AD-covariate spec (single source of truth for the male AD-model arms).
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "_shared"))
from load_ace_phenotype import load_for_projids  # noqa: E402
import arm_covariates as armcov  # noqa: E402

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    stream=sys.stdout,
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SEX_MAP = {1: "Male", 0: "Female"}
MIN_PATIENTS = 5
PHENOTYPE_COLS = [
    "projid",
    "tot_adverse_exp",
    "early_hh_ses",
    "ace_aggregate",
    "msex",
    "age_death",
    "pmi",
    "niareagansc",
    "amylsqrt",
    "tangsqrt",
]
COVARIATE_FORMULA = "{phenotype} + age_death + pmi + niareagansc"

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


def discover_h5ad_files(input_dir: Path) -> list[Path]:
    """Return sorted list of .h5ad files in *input_dir*."""
    files = sorted(input_dir.glob("*.h5ad"))
    if not files:
        raise FileNotFoundError(f"No .h5ad files found in {input_dir}")
    return files


def load_phenotype(pheno_csv: Path) -> pd.DataFrame:
    """Load and validate the ACE phenotype CSV."""
    pheno = pd.read_csv(pheno_csv)
    pheno["projid"] = pheno["projid"].astype(str)
    # Compute ace_aggregate before schema check (it's a derived column).
    if "ace_aggregate" not in pheno.columns:
        if "tot_adverse_exp" in pheno.columns and "early_hh_ses" in pheno.columns:
            z_adv = (pheno["tot_adverse_exp"] - pheno["tot_adverse_exp"].mean()) / pheno["tot_adverse_exp"].std()
            z_ses = (pheno["early_hh_ses"] - pheno["early_hh_ses"].mean()) / pheno["early_hh_ses"].std()
            pheno["ace_aggregate"] = z_adv - z_ses
    missing = set(PHENOTYPE_COLS) - set(pheno.columns)
    if missing:
        raise ValueError(f"Phenotype CSV missing columns: {missing}")
    return pheno


def pseudobulk_per_patient(
    adata: ad.AnnData,
    patient_col: str = "patient_id",
) -> ad.AnnData:
    """Aggregate raw counts per patient (sum), returning (n_patients, n_genes)."""
    # Ensure dense for groupby summation
    from scipy import sparse

    X = adata.X
    if sparse.issparse(X):
        X = X.toarray()

    df = pd.DataFrame(X, index=adata.obs_names, columns=adata.var_names)
    df[patient_col] = adata.obs[patient_col].values

    summed = df.groupby(patient_col).sum()
    pb = ad.AnnData(
        X=summed.values.astype(np.float32),
        obs=pd.DataFrame(index=summed.index),
        var=pd.DataFrame(index=summed.columns),
    )
    pb.obs.index.name = patient_col
    return pb


def normalize_pseudobulk(pb: ad.AnnData) -> ad.AnnData:
    """CPM + log1p normalization of pseudobulk AnnData."""
    sc.pp.normalize_total(pb, target_sum=1e6)
    sc.pp.log1p(pb)
    return pb


def run_mlm(pb: ad.AnnData, net: pd.DataFrame) -> pd.DataFrame:
    """Run decoupler MLM and return the activity estimate matrix."""
    # Pseudobulk AnnData has no .raw — force decoupler to use .X.
    dc.run_mlm(
        mat=pb,
        net=net,
        source="source",
        target="target",
        weight="weight",
        use_raw=False,
        verbose=True,
    )
    # MLM stores results in obsm
    act = pb.obsm["mlm_estimate"]
    if isinstance(act, pd.DataFrame):
        return act
    return pd.DataFrame(act, index=pb.obs_names, columns=pb.var_names)


def test_associations(
    activity_df: pd.DataFrame,
    pheno_df: pd.DataFrame,
    phenotype: str,
    patient_col: str = "patient_id",
    arm: Optional[str] = None,
) -> pd.DataFrame:
    """OLS regression of each activity score against *phenotype* with covariates.

    If *arm* is given, the covariate set matches that male AD-model arm
    (e.g. MaleContAD -> + amylsqrt + tangsqrt; MaleAceByAD -> + AD_binary +
    AD_binary:phenotype). Otherwise the baseline formula
    (phenotype + age_death + pmi + niareagansc) is used.
    """
    results = []
    if arm is not None:
        spec = armcov.get_spec(arm)
        ad_covars = list(spec["ad_covars"])
        formula = armcov.ols_formula(arm, phenotype, response="activity")
        base_model_cols = ["activity", phenotype, "age_death", "pmi"]
        # AD_binary is derived from niareagansc; require niareagansc non-NA for those arms
        req_cols = base_model_cols + [c for c in ad_covars if c != "AD_binary"]
        if spec["needs_ad_binary"]:
            req_cols.append("niareagansc")
    else:
        ad_covars = ["niareagansc"]
        formula = f"activity ~ {COVARIATE_FORMULA.format(phenotype=phenotype)}"
        req_cols = ["activity", phenotype, "age_death", "pmi", "niareagansc"]

    for name in activity_df.columns:
        merged = pheno_df.copy()
        merged["activity"] = activity_df[name].values

        # Derive AD_binary when the arm needs it (1 if niareagansc in {1,2}).
        if arm is not None and armcov.get_spec(arm)["needs_ad_binary"]:
            armcov.add_ad_binary(merged)

        # Drop rows with missing values in model columns
        model_cols = req_cols
        sub = merged.dropna(subset=model_cols)

        if len(sub) < MIN_PATIENTS:
            continue

        try:
            fit = smf.ols(formula, data=sub).fit()
            results.append(
                {
                    "name": name,
                    "coef": fit.params[phenotype],
                    "stderr": fit.bse[phenotype],
                    "tstat": fit.tvalues[phenotype],
                    "pvalue": fit.pvalues[phenotype],
                    "n_obs": int(fit.nobs),
                }
            )
        except Exception as exc:
            log.warning("OLS failed for %s: %s", name, exc)
            continue

    if not results:
        return pd.DataFrame(
            columns=["name", "coef", "stderr", "tstat", "pvalue", "padj", "n_obs"]
        )

    df = pd.DataFrame(results)
    # FDR correction
    _, padj, _, _ = multipletests(df["pvalue"].values, method="fdr_bh")
    df["padj"] = padj
    return df.sort_values("pvalue")


def count_targets(net: pd.DataFrame) -> dict[str, int]:
    """Count number of targets per source in a regulon/pathway network."""
    return net.groupby("source")["target"].nunique().to_dict()


# ---------------------------------------------------------------------------
# Per-cell-type processing
# ---------------------------------------------------------------------------


def process_celltype(
    h5ad_path: Path,
    pheno: pd.DataFrame,
    phenotype: str,
    output_dir: Path,
    run_progeny: bool = False,
    arm: Optional[str] = None,
    sex_filter: Optional[str] = None,
) -> list[pd.DataFrame]:
    """Process one cell type h5ad: pseudobulk, infer activity, test associations.

    *arm* selects the per-arm covariate set for the OLS (None = baseline).
    *sex_filter* restricts to a single sex label ("Male"/"Female"); None = both.
    """
    celltype = h5ad_path.stem
    log.info("Loading %s", h5ad_path)
    adata = ad.read_h5ad(h5ad_path)
    log.info("  Shape: %s", adata.shape)

    # Determine patient column
    if "patient_id" in adata.obs.columns:
        patient_col = "patient_id"
    elif "projid" in adata.obs.columns:
        patient_col = "projid"
    else:
        raise ValueError(f"No patient_id or projid column in {h5ad_path}")

    adata.obs[patient_col] = adata.obs[patient_col].astype(str)

    # Merge phenotype info onto obs for sex stratification
    obs_with_pheno = adata.obs.copy()
    obs_with_pheno = obs_with_pheno.merge(
        pheno[PHENOTYPE_COLS],
        how="left",
        left_on=patient_col,
        right_on="projid",
    )

    all_results = []

    for sex_code, sex_label in SEX_MAP.items():
        if sex_filter is not None and sex_label != sex_filter:
            continue
        log.info("  Processing %s / %s (arm=%s)", celltype, sex_label, arm or "baseline")

        # Filter to this sex
        sex_mask = obs_with_pheno["msex"] == sex_code
        patient_ids_sex = obs_with_pheno.loc[sex_mask, patient_col].unique()

        if len(patient_ids_sex) < MIN_PATIENTS:
            log.warning(
                "  Skipping %s/%s: only %d patients (need %d)",
                celltype,
                sex_label,
                len(patient_ids_sex),
                MIN_PATIENTS,
            )
            continue

        # Subset adata to this sex
        cell_mask = adata.obs[patient_col].isin(patient_ids_sex)
        adata_sex = adata[cell_mask].copy()
        log.info("    Cells: %d, Patients: %d", adata_sex.n_obs, len(patient_ids_sex))

        # Pseudobulk
        pb = pseudobulk_per_patient(adata_sex, patient_col=patient_col)
        pb = normalize_pseudobulk(pb)
        log.info("    Pseudobulk shape: %s", pb.shape)

        # Build phenotype df aligned to pseudobulk patients.
        # Filter to cohort projids and dedup on projid (pooled CSV has
        # duplicate projids from Tsai+DeJager; keep first).
        pb_pheno = pheno[pheno["projid"].isin(pb.obs_names)].drop_duplicates(subset=["projid"])
        pb_pheno = pb_pheno.set_index("projid").reindex(pb.obs_names)
        extra = [c for c in PHENOTYPE_COLS if c != "projid" and c in pb_pheno.columns]
        pb_pheno = pb_pheno[extra]

        # --- DoRothEA ---
        log.info("    Running DoRothEA (levels A, B, C)...")
        try:
            dorothea = dc.get_dorothea(organism="human", levels=["A", "B", "C"])
            dorothea_act = run_mlm(pb.copy(), dorothea)
            dorothea_res = test_associations(
                dorothea_act, pb_pheno.reset_index(), phenotype, patient_col="projid", arm=arm
            )
            n_targets = count_targets(dorothea)
            dorothea_res["n_targets"] = dorothea_res["name"].map(n_targets).fillna(0).astype(int)
            dorothea_res["source_db"] = "DoRothEA"
            dorothea_res["cell_type"] = celltype
            dorothea_res["sex"] = sex_label
            dorothea_res["phenotype"] = phenotype

            out_path = output_dir / f"tf_dorothea_{celltype}_{sex_label}.csv"
            dorothea_res.to_csv(out_path, index=False)
            log.info("    DoRothEA: %d TFs tested, saved to %s", len(dorothea_res), out_path)
            all_results.append(dorothea_res)
        except Exception as exc:
            log.error("    DoRothEA failed for %s/%s: %s", celltype, sex_label, exc)

        # --- CollecTRI ---
        log.info("    Running CollecTRI...")
        try:
            collectri = dc.get_collectri(organism="human")
            collectri_act = run_mlm(pb.copy(), collectri)
            collectri_res = test_associations(
                collectri_act, pb_pheno.reset_index(), phenotype, patient_col="projid", arm=arm
            )
            n_targets_c = count_targets(collectri)
            collectri_res["n_targets"] = collectri_res["name"].map(n_targets_c).fillna(0).astype(int)
            collectri_res["source_db"] = "CollecTRI"
            collectri_res["cell_type"] = celltype
            collectri_res["sex"] = sex_label
            collectri_res["phenotype"] = phenotype

            out_path = output_dir / f"tf_collectri_{celltype}_{sex_label}.csv"
            collectri_res.to_csv(out_path, index=False)
            log.info("    CollecTRI: %d TFs tested, saved to %s", len(collectri_res), out_path)
            all_results.append(collectri_res)
        except Exception as exc:
            log.error("    CollecTRI failed for %s/%s: %s", celltype, sex_label, exc)

        # --- PROGENy ---
        if run_progeny:
            log.info("    Running PROGENy (top 500 genes)...")
            try:
                progeny = dc.get_progeny(organism="human", top=500)
                progeny_act = run_mlm(pb.copy(), progeny)
                progeny_res = test_associations(
                    progeny_act, pb_pheno.reset_index(), phenotype, patient_col="projid", arm=arm
                )
                n_targets_p = count_targets(progeny)
                progeny_res["n_targets"] = progeny_res["name"].map(n_targets_p).fillna(0).astype(int)
                progeny_res["source_db"] = "PROGENy"
                progeny_res["cell_type"] = celltype
                progeny_res["sex"] = sex_label
                progeny_res["phenotype"] = phenotype

                out_path = output_dir / f"pathway_progeny_{celltype}_{sex_label}.csv"
                progeny_res.to_csv(out_path, index=False)
                log.info("    PROGENy: %d pathways tested, saved to %s", len(progeny_res), out_path)
                all_results.append(progeny_res)
            except Exception as exc:
                log.error("    PROGENy failed for %s/%s: %s", celltype, sex_label, exc)

    return all_results


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="ACE TF & pathway activity analysis using decoupler."
    )
    parser.add_argument(
        "--integration",
        type=str,
        default="derived_batch",
        help="Integration method label (default: derived_batch)",
    )
    parser.add_argument(
        "--phenotype",
        type=str,
        required=True,
        help="ACE phenotype column to test (e.g. tot_adverse_exp)",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Directory containing per-cell-type h5ad files",
    )
    parser.add_argument(
        "--pheno-csv",
        type=Path,
        required=True,
        help="Path to ACE phenotype CSV",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory for output CSVs",
    )
    parser.add_argument(
        "--run-progeny",
        action="store_true",
        help="Also run PROGENy pathway activity inference",
    )
    parser.add_argument(
        "--arm",
        type=str,
        default=None,
        help="Male AD-model arm (e.g. MaleContAD). Sets the OLS covariate set to "
             "match that arm; implies males-only. Default: baseline formula.",
    )
    parser.add_argument(
        "--sex-filter",
        type=str,
        default=None,
        choices=["Male", "Female"],
        help="Restrict to a single sex. Male AD-model arms use Male.",
    )
    parser.add_argument(
        "--celltypes",
        type=str,
        default=None,
        help="Comma-separated cell-type list to restrict processing (by h5ad stem).",
    )
    parser.add_argument(
        "--smoke",
        action="store_true",
        help="Smoke-test mode: process first h5ad only, skip CollecTRI and PROGENy",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    log.info("=" * 60)
    log.info("ACE TF Activity Analysis")
    log.info("=" * 60)
    log.info("  Integration: %s", args.integration)
    log.info("  Phenotype:   %s", args.phenotype)
    log.info("  Input dir:   %s", args.input_dir)
    log.info("  Pheno CSV:   %s", args.pheno_csv)
    log.info("  Output dir:  %s", args.output_dir)
    log.info("  PROGENy:     %s", args.run_progeny)
    log.info("  Smoke mode:  %s", args.smoke)
    log.info("=" * 60)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # A male AD-model arm implies males-only unless an explicit sex filter is given.
    sex_filter = args.sex_filter
    if args.arm is not None and sex_filter is None:
        sex_filter = "Male"
    if args.arm is not None:
        log.info("  Arm:         %s  (sex filter: %s)", args.arm, sex_filter)

    # Discover input files
    h5ad_files = discover_h5ad_files(args.input_dir)
    log.info("Found %d cell-type h5ad files", len(h5ad_files))

    # Restrict to requested cell types (match by h5ad stem). The pooled inhibitory
    # label "Inh" corresponds to the file "broad_Inh".
    if args.celltypes:
        want = set()
        for c in (x.strip() for x in args.celltypes.split(",")):
            want.add(c)
            want.add("broad_Inh" if c == "Inh" else c)
        h5ad_files = [f for f in h5ad_files if f.stem in want]
        log.info("Restricted to %d cell-type files: %s",
                 len(h5ad_files), ", ".join(f.stem for f in h5ad_files))
        if not h5ad_files:
            raise SystemExit(f"ERROR: --celltypes filter removed all files: {args.celltypes}")

    if args.smoke:
        h5ad_files = h5ad_files[:1]
        log.info("Smoke mode: processing only %s", h5ad_files[0].name)

    # Load phenotype
    pheno = load_phenotype(args.pheno_csv)
    log.info("Loaded phenotype data: %d patients", len(pheno))

    # In smoke mode, disable CollecTRI and PROGENy
    run_progeny = args.run_progeny and not args.smoke

    # Process each cell type
    all_results: list[pd.DataFrame] = []
    for h5ad_path in h5ad_files:
        if args.smoke:
            # In smoke mode, only run DoRothEA
            celltype = h5ad_path.stem
            log.info("Loading %s", h5ad_path)
            adata = ad.read_h5ad(h5ad_path)

            patient_col = "patient_id" if "patient_id" in adata.obs.columns else "projid"
            adata.obs[patient_col] = adata.obs[patient_col].astype(str)

            obs_with_pheno = adata.obs.copy()
            obs_with_pheno = obs_with_pheno.merge(
                pheno[PHENOTYPE_COLS],
                how="left",
                left_on=patient_col,
                right_on="projid",
            )

            for sex_code, sex_label in SEX_MAP.items():
                if sex_filter is not None and sex_label != sex_filter:
                    continue
                sex_mask = obs_with_pheno["msex"] == sex_code
                patient_ids_sex = obs_with_pheno.loc[sex_mask, patient_col].unique()
                if len(patient_ids_sex) < MIN_PATIENTS:
                    continue

                adata_sex = adata[adata.obs[patient_col].isin(patient_ids_sex)].copy()
                pb = pseudobulk_per_patient(adata_sex, patient_col=patient_col)
                pb = normalize_pseudobulk(pb)

                pb_pheno = load_for_projids(args.pheno_csv, pb.obs_names)
                extra = [c for c in PHENOTYPE_COLS if c != "projid" and c in pb_pheno.columns]
                pb_pheno = pb_pheno[extra]

                try:
                    dorothea = dc.get_dorothea(organism="human", levels=["A", "B", "C"])
                    dorothea_act = run_mlm(pb.copy(), dorothea)
                    dorothea_res = test_associations(
                        dorothea_act, pb_pheno.reset_index(), args.phenotype,
                        patient_col="projid", arm=args.arm
                    )
                    n_targets = count_targets(dorothea)
                    dorothea_res["n_targets"] = dorothea_res["name"].map(n_targets).fillna(0).astype(int)
                    dorothea_res["source_db"] = "DoRothEA"
                    dorothea_res["cell_type"] = celltype
                    dorothea_res["sex"] = sex_label
                    dorothea_res["phenotype"] = args.phenotype

                    out_path = args.output_dir / f"smoke_tf_dorothea_{celltype}_{sex_label}.csv"
                    dorothea_res.to_csv(out_path, index=False)
                    log.info("Smoke: saved %s", out_path)
                    all_results.append(dorothea_res)
                except Exception as exc:
                    log.error("Smoke DoRothEA failed: %s", exc)
        else:
            results = process_celltype(
                h5ad_path,
                pheno,
                args.phenotype,
                args.output_dir,
                run_progeny=run_progeny,
                arm=args.arm,
                sex_filter=sex_filter,
            )
            all_results.extend(results)

    # Combine summaries
    if all_results:
        combined = pd.concat(all_results, ignore_index=True)

        tf_mask = combined["source_db"].isin(["DoRothEA", "CollecTRI"])
        pw_mask = combined["source_db"] == "PROGENy"

        if tf_mask.any():
            tf_summary = combined[tf_mask].copy()
            summary_path = args.output_dir / "tf_summary.csv"
            tf_summary.to_csv(summary_path, index=False)
            log.info("TF summary: %d rows -> %s", len(tf_summary), summary_path)

        if pw_mask.any():
            pw_summary = combined[pw_mask].copy()
            summary_path = args.output_dir / "pathway_summary.csv"
            pw_summary.to_csv(summary_path, index=False)
            log.info("Pathway summary: %d rows -> %s", len(pw_summary), summary_path)

        if args.smoke:
            smoke_path = args.output_dir / f"smoke_summary_{args.phenotype}.csv"
            combined.to_csv(smoke_path, index=False)
            log.info("Smoke summary -> %s", smoke_path)

    log.info("Done.")


if __name__ == "__main__":
    main()
