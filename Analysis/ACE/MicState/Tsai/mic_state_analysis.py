#!/usr/bin/env python3
"""
ACE Microglial State Analysis

Scores microglia for established activation state signatures (homeostatic,
DAM, interferon-response, inflammatory), computes diffusion maps for
continuous state visualization, and tests ACE associations with state scores
via OLS regression.

Usage:
    python mic_state_analysis.py \
        --sex Male \
        --phenotype tot_adverse_exp \
        --input-h5ad /path/to/Mic.h5ad \
        --pheno-csv /path/to/ACE_scores.csv \
        --output-dir /path/to/output \
        [--smoke]
"""

import sys
import os
import argparse
import warnings

sys.stdout.reconfigure(line_buffering=True)
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
warnings.filterwarnings("ignore", category=FutureWarning)

import anndata as ad
import pandas as pd
import numpy as np
import scanpy as sc
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

# Shared ACE phenotype loader (handles dup-projid + computes ace_aggregate)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "_shared"))
from load_ace_phenotype import load_for_projids  # noqa: E402


# =============================================================================
# Microglial state gene signatures
# =============================================================================
STATE_SIGNATURES = {
    "homeostatic": [
        "P2RY12", "CX3CR1", "TMEM119", "CSF1R", "HEXB", "SELPLG",
        "P2RY13", "SIGLEC", "TGFBR1", "FCRLS", "OLFML3", "GPR34",
        "SLC2A5", "SALL1"
    ],
    "DAM_stage1": [
        "TYROBP", "APOE", "B2M", "FTH1", "LYZ",
        "CTSB", "CTSD", "TIMP2", "FTL"
    ],
    "DAM_stage2": [
        "TREM2", "AXL", "CST7", "LPL", "SPP1", "ITGAX",
        "CLEC7A", "LILRB4", "GPNMB", "CD9"
    ],
    "interferon_response": [
        "ISG15", "IFI44L", "IFIT1", "MX1", "OAS1",
        "IFI6", "IFITM3", "IRF7", "STAT1", "IFI44"
    ],
    "inflammatory": [
        "IL1B", "TNF", "CCL2", "CCL3", "CXCL10",
        "CCL4", "IL6", "NFKBIA", "PTGS2", "CXCL8"
    ],
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="ACE Microglial State Analysis"
    )
    parser.add_argument("--sex", required=True, choices=["Male", "Female"])
    parser.add_argument("--phenotype", default="tot_adverse_exp")
    parser.add_argument("--input-h5ad", required=True)
    parser.add_argument("--pheno-csv", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--n-diffmap-components", type=int, default=10)
    parser.add_argument("--n-state-clusters", type=int, default=4)
    parser.add_argument("--smoke", action="store_true")
    return parser.parse_args()


def load_and_filter(args):
    """Load microglia h5ad and filter by sex."""
    print(f"Loading {args.input_h5ad}...")
    adata = ad.read_h5ad(args.input_h5ad)
    print(f"  Loaded: {adata.shape[0]} cells x {adata.shape[1]} genes")

    # Identify patient ID column
    pid_col = None
    for col in ["projid", "patient_id", "sample_id"]:
        if col in adata.obs.columns:
            pid_col = col
            break
    if pid_col is None:
        raise ValueError("No patient ID column found in obs")

    adata.obs[pid_col] = adata.obs[pid_col].astype(str)

    # Use shared loader: dedups pooled-cohort CSV + computes ace_aggregate
    cohort_projids = adata.obs[pid_col].unique()
    pheno_dedup = load_for_projids(args.pheno_csv, cohort_projids)
    for col in [args.phenotype, "msex", "age_death", "pmi", "niareagansc"]:
        if col in pheno_dedup.columns:
            adata.obs[col] = adata.obs[pid_col].map(pheno_dedup[col])

    # Filter by sex
    sex_code = 1 if args.sex == "Male" else 0
    mask = (adata.obs["msex"] == sex_code) & adata.obs[args.phenotype].notna()
    adata = adata[mask].copy()

    n_patients = adata.obs[pid_col].nunique()
    print(f"  After filtering: {adata.shape[0]} cells, {n_patients} patients")

    if n_patients < 10:
        print("WARNING: Fewer than 10 patients. Skipping analysis.")
        return None, None

    # Smoke mode subsample
    if args.smoke:
        print("SMOKE MODE: subsampling to 3000 cells")
        np.random.seed(42)
        n_keep = min(3000, adata.shape[0])
        idx = np.random.choice(adata.shape[0], n_keep, replace=False)
        adata = adata[idx].copy()

    return adata, pid_col


def score_states(adata):
    """Score each cell for microglial state signatures."""
    print("Scoring microglial states...")

    # Normalize for scoring
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Score each state signature
    for state_name, genes in STATE_SIGNATURES.items():
        # Filter to genes present in the data
        present_genes = [g for g in genes if g in adata.var_names]
        if len(present_genes) < 3:
            print(f"  WARNING: Only {len(present_genes)} genes found for "
                  f"{state_name}. Skipping.")
            adata.obs[f"score_{state_name}"] = np.nan
            continue

        sc.tl.score_genes(adata, gene_list=present_genes,
                          score_name=f"score_{state_name}")
        print(f"  {state_name}: {len(present_genes)}/{len(genes)} genes used")

    return adata


def compute_diffusion_map(adata, n_components=10):
    """Compute diffusion map for continuous state visualization."""
    print("Computing diffusion map...")

    # HVG selection and PCA
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    sc.pp.pca(adata, n_comps=30, use_highly_variable=True)

    # Neighbors and diffusion map
    sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30)
    sc.tl.diffmap(adata, n_comps=n_components)

    return adata


def aggregate_to_patient(adata, pid_col, phenotype):
    """Aggregate state scores to patient level."""
    print("Aggregating scores to patient level...")

    score_cols = [c for c in adata.obs.columns if c.startswith("score_")]
    covariate_cols = ["age_death", "pmi", "niareagansc", phenotype]

    # Group by patient and compute mean scores
    patient_df = (
        adata.obs
        .groupby(pid_col)
        .agg({
            **{col: "mean" for col in score_cols if col in adata.obs.columns},
            **{col: "first" for col in covariate_cols if col in adata.obs.columns},
            pid_col: "count"
        })
        .rename(columns={pid_col: "n_cells"})
    )

    # Reset index to get pid_col back as a column
    patient_df = patient_df.reset_index()

    print(f"  {len(patient_df)} patients with mean state scores")
    return patient_df


def test_associations(patient_df, phenotype):
    """Test ACE association with each state score via OLS."""
    print("Testing ACE-state associations...")

    score_cols = [c for c in patient_df.columns if c.startswith("score_")]
    results = []

    for score_col in score_cols:
        state_name = score_col.replace("score_", "")

        # Drop rows with NaN in this score
        df = patient_df[[score_col, phenotype, "age_death", "pmi",
                         "niareagansc"]].dropna()

        if len(df) < 10:
            print(f"  {state_name}: insufficient data ({len(df)} patients)")
            continue

        # Z-score covariates for numerical stability
        for col in ["age_death", "pmi"]:
            if df[col].std() > 0:
                df[col] = (df[col] - df[col].mean()) / df[col].std()

        # Rename for formula compatibility
        df = df.rename(columns={score_col: "score", phenotype: "pheno"})

        try:
            model = smf.ols("score ~ pheno + age_death + pmi + niareagansc",
                            data=df)
            fit = model.fit()

            results.append({
                "state": state_name,
                "coef": fit.params["pheno"],
                "stderr": fit.bse["pheno"],
                "tstat": fit.tvalues["pheno"],
                "pvalue": fit.pvalues["pheno"],
                "n_patients": len(df),
                "r_squared": fit.rsquared,
            })
            print(f"  {state_name}: coef={fit.params['pheno']:.4f}, "
                  f"p={fit.pvalues['pheno']:.4e}")
        except Exception as e:
            print(f"  {state_name}: regression failed: {e}")

    results_df = pd.DataFrame(results)
    if len(results_df) > 0:
        _, pvals_fdr, _, _ = multipletests(
            results_df["pvalue"], alpha=0.05, method="fdr_bh"
        )
        results_df["padj"] = pvals_fdr

    return results_df


def cluster_states(adata, n_clusters=4):
    """Cluster cells by state scores using k-means."""
    from sklearn.cluster import KMeans

    print(f"Clustering cells into {n_clusters} states...")

    score_cols = [c for c in adata.obs.columns
                  if c.startswith("score_") and not adata.obs[c].isna().all()]

    if len(score_cols) < 2:
        print("  WARNING: Not enough state scores for clustering")
        return None

    score_matrix = adata.obs[score_cols].dropna()

    # Standardize
    score_std = (score_matrix - score_matrix.mean()) / score_matrix.std()

    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    labels = kmeans.fit_predict(score_std)

    adata.obs.loc[score_matrix.index, "state_cluster"] = labels.astype(str)

    return adata


def test_state_proportions(adata, pid_col, phenotype):
    """Test whether ACE associates with different state cluster proportions."""
    print("Testing state proportion associations...")

    if "state_cluster" not in adata.obs.columns:
        return None

    # Compute per-patient proportions
    counts = (
        adata.obs
        .groupby([pid_col, "state_cluster"])
        .size()
        .unstack(fill_value=0)
    )
    proportions = counts.div(counts.sum(axis=1), axis=0)

    # Merge with phenotype
    meta_unique = (
        adata.obs
        .groupby(pid_col)[[phenotype, "age_death", "pmi", "niareagansc"]]
        .first()
    )
    prop_pheno = proportions.join(meta_unique)

    results = []
    for cluster in counts.columns:
        col = str(cluster)
        df = prop_pheno[[col, phenotype, "age_death", "pmi",
                         "niareagansc"]].dropna()
        if len(df) < 10:
            continue

        df = df.rename(columns={col: "proportion", phenotype: "pheno"})
        for cov in ["age_death", "pmi"]:
            if df[cov].std() > 0:
                df[cov] = (df[cov] - df[cov].mean()) / df[cov].std()

        try:
            model = smf.ols(
                "proportion ~ pheno + age_death + pmi + niareagansc", data=df
            )
            fit = model.fit()
            results.append({
                "cluster": cluster,
                "coef": fit.params["pheno"],
                "stderr": fit.bse["pheno"],
                "pvalue": fit.pvalues["pheno"],
                "n_patients": len(df),
                "mean_proportion": df["proportion"].mean(),
            })
        except Exception as e:
            print(f"  Cluster {cluster} failed: {e}")

    results_df = pd.DataFrame(results)
    if len(results_df) > 0:
        _, pvals_fdr, _, _ = multipletests(
            results_df["pvalue"], alpha=0.05, method="fdr_bh"
        )
        results_df["padj"] = pvals_fdr

    return results_df


def generate_figures(adata, patient_df, regression_df, proportion_df,
                     output_dir, phenotype, sex):
    """Generate visualization figures."""
    fig_dir = os.path.join(output_dir, "..", "..", "figures")
    os.makedirs(fig_dir, exist_ok=True)

    score_cols = [c for c in patient_df.columns if c.startswith("score_")]

    # 1. State score boxplots by ACE group
    if len(score_cols) > 0:
        print("Generating state score boxplots...")
        median_val = patient_df[phenotype].median()
        patient_df["ace_group"] = np.where(
            patient_df[phenotype] > median_val, "ACE-high", "ACE-low"
        )

        fig, axes = plt.subplots(1, len(score_cols), figsize=(4 * len(score_cols), 5))
        if len(score_cols) == 1:
            axes = [axes]

        for ax, col in zip(axes, score_cols):
            state = col.replace("score_", "")
            data = patient_df[[col, "ace_group"]].dropna()
            sns.boxplot(data=data, x="ace_group", y=col, ax=ax,
                        palette={"ACE-high": "#D73027", "ACE-low": "#4575B4"})
            ax.set_title(state, fontsize=12)
            ax.set_xlabel("")
            ax.set_ylabel("Score")

        plt.suptitle(f"Microglial State Scores - {sex}", fontsize=14, y=1.02)
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f"state_score_boxplots_{sex}.pdf"),
                    bbox_inches="tight", dpi=150)
        plt.close()

    # 2. Regression results barplot
    if regression_df is not None and len(regression_df) > 0:
        print("Generating regression coefficient barplot...")
        fig, ax = plt.subplots(figsize=(8, 5))
        colors = ["#D73027" if p < 0.05 else "#999999"
                  for p in regression_df["padj"]]

        ax.barh(regression_df["state"], regression_df["coef"], color=colors,
                edgecolor="black", linewidth=0.5)
        ax.axvline(0, color="black", linewidth=0.5)
        ax.set_xlabel(f"Coefficient ({phenotype})")
        ax.set_title(f"Microglial State ~ ACE Association ({sex})")

        # Add significance markers
        for i, row in regression_df.iterrows():
            if row["padj"] < 0.05:
                ax.text(row["coef"], i, " *", va="center", fontsize=14,
                        fontweight="bold")

        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f"state_regression_{sex}.pdf"),
                    bbox_inches="tight", dpi=150)
        plt.close()

    # 3. Diffusion map colored by state scores
    if "X_diffmap" in adata.obsm:
        print("Generating diffusion map plots...")
        for col in [c for c in adata.obs.columns if c.startswith("score_")]:
            if adata.obs[col].isna().all():
                continue
            state = col.replace("score_", "")
            fig, ax = plt.subplots(figsize=(7, 6))
            sc_vals = adata.obs[col].values
            scatter = ax.scatter(
                adata.obsm["X_diffmap"][:, 0],
                adata.obsm["X_diffmap"][:, 1],
                c=sc_vals, cmap="RdBu_r", s=1, alpha=0.5,
                vmin=np.nanpercentile(sc_vals, 5),
                vmax=np.nanpercentile(sc_vals, 95)
            )
            plt.colorbar(scatter, label=f"{state} score")
            ax.set_xlabel("Diffusion Component 1")
            ax.set_ylabel("Diffusion Component 2")
            ax.set_title(f"Microglia - {state} ({sex})")
            plt.tight_layout()
            plt.savefig(
                os.path.join(fig_dir, f"diffmap_{state}_{sex}.pdf"),
                bbox_inches="tight", dpi=150
            )
            plt.close()

    # 4. State proportion barplot
    if proportion_df is not None and len(proportion_df) > 0:
        print("Generating state proportion results...")
        fig, ax = plt.subplots(figsize=(6, 4))
        colors = ["#D73027" if p < 0.05 else "#999999"
                  for p in proportion_df["padj"]]
        ax.barh(proportion_df["cluster"].astype(str), proportion_df["coef"],
                color=colors, edgecolor="black", linewidth=0.5)
        ax.axvline(0, color="black", linewidth=0.5)
        ax.set_xlabel(f"Proportion ~ {phenotype} coefficient")
        ax.set_title(f"State Cluster Proportions ~ ACE ({sex})")
        plt.tight_layout()
        plt.savefig(os.path.join(fig_dir, f"state_proportions_{sex}.pdf"),
                    bbox_inches="tight", dpi=150)
        plt.close()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    print("=" * 60)
    print("ACE Microglial State Analysis")
    print("=" * 60)
    print(f"Sex:       {args.sex}")
    print(f"Phenotype: {args.phenotype}")
    print(f"Input:     {args.input_h5ad}")
    print(f"Output:    {args.output_dir}")
    print()

    # Load and filter
    adata, pid_col = load_and_filter(args)
    if adata is None:
        return

    # Score states
    adata = score_states(adata)

    # Save per-cell scores
    score_cols = [c for c in adata.obs.columns if c.startswith("score_")]
    cell_scores = adata.obs[[pid_col] + score_cols].copy()
    cell_scores.to_csv(os.path.join(args.output_dir,
                                     "state_scores_per_cell.csv"))

    # Compute diffusion map
    adata = compute_diffusion_map(adata, args.n_diffmap_components)

    # Save diffusion map coordinates
    diffmap_df = pd.DataFrame(
        adata.obsm["X_diffmap"][:, :3],
        columns=["DC1", "DC2", "DC3"],
        index=adata.obs_names
    )
    diffmap_df.to_csv(os.path.join(args.output_dir,
                                    "diffmap_coordinates.csv"))

    # Aggregate to patient level
    patient_df = aggregate_to_patient(adata, pid_col, args.phenotype)
    patient_df.to_csv(os.path.join(args.output_dir,
                                    "state_scores_per_patient.csv"),
                       index=False)

    # Test associations
    regression_df = test_associations(patient_df, args.phenotype)
    regression_df.to_csv(os.path.join(args.output_dir,
                                       "state_regression_results.csv"),
                          index=False)

    # Cluster states
    adata = cluster_states(adata, args.n_state_clusters)

    # Test state proportions
    proportion_df = test_state_proportions(adata, pid_col, args.phenotype)
    if proportion_df is not None:
        proportion_df.to_csv(os.path.join(args.output_dir,
                                           "state_proportion_test.csv"),
                              index=False)

    # Save cluster assignments
    if "state_cluster" in adata.obs.columns:
        cluster_counts = (
            adata.obs
            .groupby([pid_col, "state_cluster"])
            .size()
            .unstack(fill_value=0)
        )
        cluster_counts.to_csv(os.path.join(args.output_dir,
                                            "state_cluster_proportions.csv"))

    # Generate figures
    generate_figures(adata, patient_df, regression_df, proportion_df,
                     args.output_dir, args.phenotype, args.sex)

    # Summary
    print()
    print("=" * 60)
    print("Microglial State Analysis Complete")
    print("=" * 60)
    if len(regression_df) > 0:
        sig = regression_df[regression_df["padj"] < 0.05]
        print(f"Significant state associations (padj < 0.05): {len(sig)}")
        for _, row in sig.iterrows():
            direction = "higher" if row["coef"] > 0 else "lower"
            print(f"  {row['state']}: {direction} with ACE "
                  f"(coef={row['coef']:.4f}, padj={row['padj']:.4e})")
    print(f"Output: {args.output_dir}")


if __name__ == "__main__":
    main()
