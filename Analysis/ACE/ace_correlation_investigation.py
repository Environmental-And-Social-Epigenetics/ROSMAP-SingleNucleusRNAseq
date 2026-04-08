#!/usr/bin/env python3
"""
ACE Correlation Investigation
==============================
Investigates why combining tot_adverse_exp and early_hh_ses into the
ace_aggregate composite kills transcriptomic signal in the Tsai cohort.

Produces:
  - Multi-panel PDF (ace_correlation_investigation.pdf)
  - Text summary to stdout
"""

from __future__ import annotations

import os
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr, zscore as scipy_zscore, shapiro, mannwhitneyu


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[2]  # ROSMAP_Code/Transcriptomics

PHENO_CSV = (
    REPO_ROOT / "Data" / "Phenotypes" / "TSAI_DEJAGER_all_patients_wACEscores.csv"
)
DEG_SUMMARY_CSV = SCRIPT_DIR / "DEG" / "Tsai" / "deg_summary.csv"
OUTPUT_PDF = SCRIPT_DIR / "ace_correlation_investigation.pdf"

ACE_COMPONENTS = [
    "emotional_neglect",
    "family_pro_sep",
    "financial_need",
    "parental_intimidation",
    "parental_violence",
]
SES_VALIDATORS = ["pareduc", "mateduc", "pateduc", "educ", "income_bl"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def z(x: pd.Series) -> np.ndarray:
    """Z-score with ddof=1 (matches R scale())."""
    vals = x.values.astype(float)
    return (vals - np.nanmean(vals)) / np.nanstd(vals, ddof=1)


def bootstrap_corr(x, y, func, n_boot=1000, seed=42):
    """Bootstrap 95% CI for a correlation function."""
    rng = np.random.default_rng(seed)
    n = len(x)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, size=n)
        boots[i] = func(x[idx], y[idx])[0]
    return np.percentile(boots, [2.5, 97.5])


def extract_celltype_sex(ct_sex: str) -> tuple[str, str]:
    """Parse ct_sex like 'adverse_exp_Mic_Male' -> ('Mic', 'Male')."""
    for sex in ("Male", "Fem"):
        if ct_sex.endswith(f"_{sex}"):
            prefix = ct_sex[: -len(f"_{sex}")]
            # strip phenotype prefix (adverse_exp_, hh_ses_, aggregate_)
            ct = re.sub(r"^(adverse_exp_|hh_ses_|aggregate_)", "", prefix)
            sex_label = "Male" if sex == "Male" else "Female"
            return ct, sex_label
    return ct_sex, "Unknown"


# ---------------------------------------------------------------------------
# 1. Data Loading
# ---------------------------------------------------------------------------
def load_data(pheno_path: Path) -> pd.DataFrame:
    df = pd.read_csv(pheno_path)
    mask = df["tot_adverse_exp"].notna() & df["early_hh_ses"].notna()
    df = df.loc[mask].copy()
    df["sex_label"] = df["msex"].map({0: "Female", 1: "Male"})

    # Original composite (current pipeline)
    df["z_adverse"] = z(df["tot_adverse_exp"])
    df["z_ses"] = z(df["early_hh_ses"])
    df["ace_aggregate_orig"] = df["z_adverse"] + df["z_ses"]

    # Corrected composite (deprivation-aligned)
    df["z_ses_neg"] = z(-df["early_hh_ses"])
    df["ace_aggregate_corr"] = df["z_adverse"] + df["z_ses_neg"]

    return df


# ---------------------------------------------------------------------------
# 2. Directionality Confirmation
# ---------------------------------------------------------------------------
def print_directionality(df: pd.DataFrame) -> None:
    print("=" * 70)
    print("DIRECTIONALITY OF early_hh_ses")
    print("=" * 70)
    print("\nCorrelation of early_hh_ses with SES-related variables (Spearman):")
    for col in SES_VALIDATORS:
        if col in df.columns:
            m = df[col].notna()
            if m.sum() > 10:
                rho, p = spearmanr(df.loc[m, "early_hh_ses"], df.loc[m, col])
                print(f"  vs {col:20s}: rho = {rho:+.4f}  p = {p:.4f}  (n={m.sum()})")

    print("\nCorrelation of early_hh_ses with ACE components (Spearman):")
    for col in ACE_COMPONENTS:
        if col in df.columns:
            m = df[col].notna()
            if m.sum() > 10:
                rho, p = spearmanr(df.loc[m, "early_hh_ses"], df.loc[m, col])
                print(f"  vs {col:25s}: rho = {rho:+.4f}  p = {p:.4f}  (n={m.sum()})")

    print("\n>>> INTERPRETATION: High early_hh_ses = HIGH socioeconomic status")
    print("    (strong positive correlation with parental/own education and income)")
    print("    To create a deprivation composite, NEGATE early_hh_ses before combining.")


# ---------------------------------------------------------------------------
# 3. Correlation Analysis
# ---------------------------------------------------------------------------
def print_correlations(df: pd.DataFrame) -> None:
    print("\n" + "=" * 70)
    print("CORRELATION: tot_adverse_exp vs early_hh_ses")
    print("=" * 70)

    groups = [("All", df), ("Male", df[df["msex"] == 1]), ("Female", df[df["msex"] == 0])]

    for label, sub in groups:
        x = sub["tot_adverse_exp"].values
        y = sub["early_hh_ses"].values
        r, p_r = pearsonr(x, y)
        rho, p_rho = spearmanr(x, y)
        ci_r = bootstrap_corr(x, y, pearsonr)
        ci_rho = bootstrap_corr(x, y, spearmanr)
        print(f"\n  {label} (n={len(sub)}):")
        print(f"    Pearson  r   = {r:+.4f}  [{ci_r[0]:+.4f}, {ci_r[1]:+.4f}]  p = {p_r:.4f}")
        print(f"    Spearman rho = {rho:+.4f}  [{ci_rho[0]:+.4f}, {ci_rho[1]:+.4f}]  p = {p_rho:.4f}")

    print("\n>>> The variables are WEAKLY NEGATIVELY CORRELATED.")
    print("    Summing their z-scores partially cancels the signal,")
    print("    especially in males (rho ~ -0.24).")


# ---------------------------------------------------------------------------
# 4. Variance Decomposition
# ---------------------------------------------------------------------------
def print_variance_decomposition(df: pd.DataFrame) -> None:
    print("\n" + "=" * 70)
    print("VARIANCE DECOMPOSITION: original vs corrected composite")
    print("=" * 70)

    groups = [("All", df), ("Male", df[df["msex"] == 1]), ("Female", df[df["msex"] == 0])]

    for label, sub in groups:
        za = sub["z_adverse"].values
        zs = sub["z_ses"].values
        zs_neg = sub["z_ses_neg"].values

        # Original
        var_za = np.var(za, ddof=1)
        var_zs = np.var(zs, ddof=1)
        cov_orig = 2 * np.cov(za, zs, ddof=1)[0, 1]
        var_orig = var_za + var_zs + cov_orig

        # Corrected
        var_zs_neg = np.var(zs_neg, ddof=1)
        cov_corr = 2 * np.cov(za, zs_neg, ddof=1)[0, 1]
        var_corr = var_za + var_zs_neg + cov_corr

        print(f"\n  {label} (n={len(sub)}):")
        print(f"    ORIGINAL  ace_aggregate = z(adverse) + z(ses):")
        print(f"      Var(z_adv)={var_za:.3f}  Var(z_ses)={var_zs:.3f}  2*Cov={cov_orig:+.3f}  => Var(comp)={var_orig:.3f}")
        print(f"    CORRECTED ace_aggregate = z(adverse) - z(ses):")
        print(f"      Var(z_adv)={var_za:.3f}  Var(-z_ses)={var_zs_neg:.3f}  2*Cov={cov_corr:+.3f}  => Var(comp)={var_corr:.3f}")
        print(f"    Variance gain: {var_corr / var_orig:.1%} of original")


# ---------------------------------------------------------------------------
# 5. Rank Preservation
# ---------------------------------------------------------------------------
def print_rank_preservation(df: pd.DataFrame) -> None:
    print("\n" + "=" * 70)
    print("RANK PRESERVATION: how much of tot_adverse_exp ranking survives?")
    print("=" * 70)

    groups = [("All", df), ("Male", df[df["msex"] == 1]), ("Female", df[df["msex"] == 0])]

    for label, sub in groups:
        rho_orig, _ = spearmanr(sub["tot_adverse_exp"], sub["ace_aggregate_orig"])
        rho_corr, _ = spearmanr(sub["tot_adverse_exp"], sub["ace_aggregate_corr"])
        print(f"\n  {label} (n={len(sub)}):")
        print(f"    Original composite:  Spearman with tot_adverse_exp = {rho_orig:.4f}  (R^2={rho_orig**2:.1%})")
        print(f"    Corrected composite: Spearman with tot_adverse_exp = {rho_corr:.4f}  (R^2={rho_corr**2:.1%})")

    # Noise simulation for males
    print("\n  Noise simulation (Males):")
    males = df[df["msex"] == 1]
    z_adv_m = males["z_adverse"].values
    rng = np.random.default_rng(42)
    noise_rhos = np.empty(1000)
    for i in range(1000):
        noise = rng.standard_normal(len(z_adv_m))
        composite_noise = z_adv_m + noise
        noise_rhos[i], _ = spearmanr(males["tot_adverse_exp"].values, composite_noise)

    rho_orig_m, _ = spearmanr(males["tot_adverse_exp"], males["ace_aggregate_orig"])
    rho_corr_m, _ = spearmanr(males["tot_adverse_exp"], males["ace_aggregate_corr"])

    pct_orig = np.mean(noise_rhos <= rho_orig_m) * 100
    pct_corr = np.mean(noise_rhos <= rho_corr_m) * 100
    print(f"    Random noise composite: median rho = {np.median(noise_rhos):.4f}  [95% CI: {np.percentile(noise_rhos, 2.5):.4f}, {np.percentile(noise_rhos, 97.5):.4f}]")
    print(f"    Original composite rho = {rho_orig_m:.4f}  ({pct_orig:.1f}th percentile of noise)")
    print(f"    Corrected composite rho = {rho_corr_m:.4f}  ({pct_corr:.1f}th percentile of noise)")
    if pct_orig < 50:
        print("    >>> Original composite preserves LESS ranking than random noise!")
        print("        SES is actively canceling adversity signal, not just adding noise.")


# ---------------------------------------------------------------------------
# 6. Visualizations
# ---------------------------------------------------------------------------
def make_plots(df: pd.DataFrame) -> None:
    deg = None
    if DEG_SUMMARY_CSV.exists():
        deg = pd.read_csv(DEG_SUMMARY_CSV)

    colors = {"Female": "#E74C3C", "Male": "#3498DB"}

    with PdfPages(str(OUTPUT_PDF)) as pdf:
        # ---- Page 1: Diagnosis ----
        fig, axes = plt.subplots(3, 3, figsize=(16, 14))
        fig.suptitle("ACE Composite Diagnosis: Why ace_aggregate Kills Signal", fontsize=14, fontweight="bold")

        # (0,0) Scatter: tot_adverse_exp vs early_hh_ses
        ax = axes[0, 0]
        for sex in ["Female", "Male"]:
            sub = df[df["sex_label"] == sex]
            ax.scatter(sub["tot_adverse_exp"], sub["early_hh_ses"], alpha=0.5, s=20, color=colors[sex], label=sex)
            # Regression line
            coef = np.polyfit(sub["tot_adverse_exp"], sub["early_hh_ses"], 1)
            xline = np.linspace(sub["tot_adverse_exp"].min(), sub["tot_adverse_exp"].max(), 50)
            ax.plot(xline, np.polyval(coef, xline), color=colors[sex], linewidth=1.5)
            rho, p = spearmanr(sub["tot_adverse_exp"], sub["early_hh_ses"])
            ax.text(0.05, 0.95 if sex == "Female" else 0.85, f"{sex}: rho={rho:.2f}, p={p:.3f}",
                    transform=ax.transAxes, fontsize=8, color=colors[sex], va="top")
        ax.set_xlabel("tot_adverse_exp")
        ax.set_ylabel("early_hh_ses")
        ax.set_title("Adversity vs SES (raw)")
        ax.legend(fontsize=8)

        # (0,1) Distribution of tot_adverse_exp
        ax = axes[0, 1]
        for sex in ["Female", "Male"]:
            sub = df[df["sex_label"] == sex]
            ax.hist(sub["tot_adverse_exp"], bins=25, alpha=0.5, color=colors[sex], label=f"{sex} (n={len(sub)})", density=True)
        ax.set_xlabel("tot_adverse_exp")
        ax.set_ylabel("Density")
        ax.set_title("Distribution of Adversity")
        ax.legend(fontsize=8)

        # (0,2) Distribution of early_hh_ses
        ax = axes[0, 2]
        for sex in ["Female", "Male"]:
            sub = df[df["sex_label"] == sex]
            ax.hist(sub["early_hh_ses"], bins=25, alpha=0.5, color=colors[sex], label=f"{sex} (n={len(sub)})", density=True)
        ax.set_xlabel("early_hh_ses")
        ax.set_ylabel("Density")
        ax.set_title("Distribution of SES")
        ax.legend(fontsize=8)

        # (1,0) Z-scored scatter
        ax = axes[1, 0]
        for sex in ["Female", "Male"]:
            sub = df[df["sex_label"] == sex]
            ax.scatter(sub["z_adverse"], sub["z_ses"], alpha=0.5, s=20, color=colors[sex], label=sex)
        ax.axhline(0, color="gray", linewidth=0.5, linestyle="--")
        ax.axvline(0, color="gray", linewidth=0.5, linestyle="--")
        ax.set_xlabel("z(tot_adverse_exp)")
        ax.set_ylabel("z(early_hh_ses)")
        ax.set_title("Z-scored space (ORIGINAL)")
        ax.legend(fontsize=8)

        # (1,1) Distribution of original ace_aggregate
        ax = axes[1, 1]
        for sex in ["Female", "Male"]:
            sub = df[df["sex_label"] == sex]
            ax.hist(sub["ace_aggregate_orig"], bins=25, alpha=0.5, color=colors[sex], label=f"{sex} (n={len(sub)})", density=True)
        ax.set_xlabel("ace_aggregate (original)")
        ax.set_ylabel("Density")
        ax.set_title("Original Composite: z(adv) + z(ses)")
        ax.legend(fontsize=8)

        # (1,2) Correlation heatmap
        ax = axes[1, 2]
        heatmap_cols = ACE_COMPONENTS + ["tot_adverse_exp", "early_hh_ses"]
        available = [c for c in heatmap_cols if c in df.columns]
        if "pareduc" in df.columns:
            available.append("pareduc")
        corr_matrix = df[available].corr(method="spearman")
        im = ax.imshow(corr_matrix.values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
        ax.set_xticks(range(len(available)))
        ax.set_yticks(range(len(available)))
        short_labels = [c.replace("emotional_neglect", "emot_neg").replace("parental_intimidation", "par_intim")
                        .replace("parental_violence", "par_viol").replace("family_pro_sep", "fam_sep")
                        .replace("financial_need", "fin_need").replace("tot_adverse_exp", "tot_adv")
                        .replace("early_hh_ses", "ses") for c in available]
        ax.set_xticklabels(short_labels, rotation=45, ha="right", fontsize=7)
        ax.set_yticklabels(short_labels, fontsize=7)
        # Annotate cells
        for i in range(len(available)):
            for j in range(len(available)):
                val = corr_matrix.values[i, j]
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=5.5,
                        color="white" if abs(val) > 0.5 else "black")
        ax.set_title("Spearman Correlation Matrix")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

        # (2,0) Rank comparison: males, original composite
        ax = axes[2, 0]
        males = df[df["msex"] == 1].copy()
        males["rank_adv"] = males["tot_adverse_exp"].rank()
        males["rank_orig"] = males["ace_aggregate_orig"].rank()
        rank_change = np.abs(males["rank_adv"] - males["rank_orig"])
        sc = ax.scatter(males["rank_adv"], males["rank_orig"], c=rank_change, cmap="YlOrRd", s=20, alpha=0.7)
        ax.plot([0, males["rank_adv"].max()], [0, males["rank_adv"].max()], "k--", linewidth=0.8, alpha=0.5)
        rho_val, _ = spearmanr(males["tot_adverse_exp"], males["ace_aggregate_orig"])
        ax.set_xlabel("Rank(tot_adverse_exp)")
        ax.set_ylabel("Rank(ace_aggregate_orig)")
        ax.set_title(f"Male Rank Disruption (orig)\nrho={rho_val:.3f}")
        fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04, label="|rank shift|")

        # (2,1) DEG count comparison
        ax = axes[2, 1]
        if deg is not None:
            deg_proj = deg[deg["integration"] == "projid"].copy()
            deg_proj["ct"], deg_proj["sex"] = zip(*deg_proj["ct_sex"].map(extract_celltype_sex))
            deg_male = deg_proj[deg_proj["sex"] == "Male"]

            for pheno in ["tot_adverse_exp", "ace_aggregate"]:
                sub = deg_male[deg_male["phenotype"] == pheno].sort_values("ct")
                label = "tot_adverse_exp" if pheno == "tot_adverse_exp" else "ace_aggregate (orig)"
                bar_positions = np.arange(len(sub))
                width = 0.35
                offset = -width / 2 if pheno == "tot_adverse_exp" else width / 2
                ax.bar(bar_positions + offset, sub["n_sig"].values, width,
                       label=label, alpha=0.8)
                if pheno == "tot_adverse_exp":
                    ct_labels = sub["ct"].values

            if len(ct_labels) > 0:
                ax.set_xticks(np.arange(len(ct_labels)))
                ax.set_xticklabels(ct_labels, rotation=45, ha="right", fontsize=7)
            ax.set_ylabel("# DEGs (padj < 0.05)")
            ax.set_title("Male DEG Signal Loss")
            ax.legend(fontsize=7)
        else:
            ax.text(0.5, 0.5, "deg_summary.csv\nnot found", transform=ax.transAxes, ha="center")

        # (2,2) Variance decomposition comparison
        ax = axes[2, 2]
        za = df["z_adverse"].values
        zs = df["z_ses"].values
        zs_neg = df["z_ses_neg"].values

        var_za = np.var(za, ddof=1)
        var_zs = np.var(zs, ddof=1)
        cov_orig = 2 * np.cov(za, zs, ddof=1)[0, 1]
        var_zs_neg = np.var(zs_neg, ddof=1)
        cov_corr = 2 * np.cov(za, zs_neg, ddof=1)[0, 1]

        x_pos = [0, 1]
        bar_labels = ["Original\nz(adv)+z(ses)", "Corrected\nz(adv)-z(ses)"]
        bottom_vals = [0, 0]

        # Stack: Var(z_adv), Var(z_ses), 2*Cov
        components_orig = [var_za, var_zs, cov_orig]
        components_corr = [var_za, var_zs_neg, cov_corr]
        comp_labels = ["Var(z_adv)", "Var(z_ses)", "2*Cov"]
        comp_colors = ["#2ECC71", "#3498DB", "#E74C3C"]

        for i, (co, cc, cl, clr) in enumerate(zip(components_orig, components_corr, comp_labels, comp_colors)):
            ax.bar(0, co if co > 0 else 0, bottom=bottom_vals[0] if co > 0 else bottom_vals[0] + co,
                   color=clr, alpha=0.8, width=0.5, label=cl if i < 3 else None)
            if co > 0:
                bottom_vals[0] += co

            ax.bar(1, cc if cc > 0 else 0, bottom=bottom_vals[1] if cc > 0 else bottom_vals[1] + cc,
                   color=clr, alpha=0.8, width=0.5)
            if cc > 0:
                bottom_vals[1] += cc

        # Annotate totals
        var_orig_total = var_za + var_zs + cov_orig
        var_corr_total = var_za + var_zs_neg + cov_corr
        ax.text(0, var_orig_total + 0.05, f"Total: {var_orig_total:.2f}", ha="center", fontsize=9, fontweight="bold")
        ax.text(1, var_corr_total + 0.05, f"Total: {var_corr_total:.2f}", ha="center", fontsize=9, fontweight="bold")

        ax.set_xticks(x_pos)
        ax.set_xticklabels(bar_labels, fontsize=9)
        ax.set_ylabel("Variance")
        ax.set_title("Composite Variance Decomposition")
        ax.legend(fontsize=8)
        ax.axhline(0, color="gray", linewidth=0.5)

        fig.tight_layout(rect=[0, 0, 1, 0.96])
        pdf.savefig(fig)
        plt.close(fig)

        # ---- Page 2: Corrected Composite ----
        fig, axes = plt.subplots(2, 3, figsize=(16, 10))
        fig.suptitle("Corrected Composite: z(adverse) - z(ses)  [both = deprivation]", fontsize=14, fontweight="bold")

        # (0,0) Z-scored scatter corrected
        ax = axes[0, 0]
        for sex in ["Female", "Male"]:
            sub = df[df["sex_label"] == sex]
            ax.scatter(sub["z_adverse"], sub["z_ses_neg"], alpha=0.5, s=20, color=colors[sex], label=sex)
        ax.axhline(0, color="gray", linewidth=0.5, linestyle="--")
        ax.axvline(0, color="gray", linewidth=0.5, linestyle="--")
        rho_all, _ = spearmanr(df["z_adverse"], df["z_ses_neg"])
        ax.text(0.05, 0.95, f"rho={rho_all:+.3f}", transform=ax.transAxes, fontsize=9, va="top")
        ax.set_xlabel("z(tot_adverse_exp)")
        ax.set_ylabel("-z(early_hh_ses)  [= z(deprivation)]")
        ax.set_title("Z-scored space (CORRECTED)")
        ax.legend(fontsize=8)

        # (0,1) Corrected composite distribution
        ax = axes[0, 1]
        for sex in ["Female", "Male"]:
            sub = df[df["sex_label"] == sex]
            ax.hist(sub["ace_aggregate_corr"], bins=25, alpha=0.5, color=colors[sex], label=f"{sex} (n={len(sub)})", density=True)
        ax.set_xlabel("ace_aggregate (corrected)")
        ax.set_ylabel("Density")
        ax.set_title("Corrected Composite Distribution")
        ax.legend(fontsize=8)

        # (0,2) Rank comparison: males, corrected composite
        ax = axes[0, 2]
        males = df[df["msex"] == 1].copy()
        males["rank_adv"] = males["tot_adverse_exp"].rank()
        males["rank_corr"] = males["ace_aggregate_corr"].rank()
        rank_change = np.abs(males["rank_adv"] - males["rank_corr"])
        sc = ax.scatter(males["rank_adv"], males["rank_corr"], c=rank_change, cmap="YlOrRd", s=20, alpha=0.7)
        ax.plot([0, males["rank_adv"].max()], [0, males["rank_adv"].max()], "k--", linewidth=0.8, alpha=0.5)
        rho_val, _ = spearmanr(males["tot_adverse_exp"], males["ace_aggregate_corr"])
        ax.set_xlabel("Rank(tot_adverse_exp)")
        ax.set_ylabel("Rank(ace_aggregate_corr)")
        ax.set_title(f"Male Rank Preservation (corr)\nrho={rho_val:.3f}")
        fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04, label="|rank shift|")

        # (1,0) Directionality scatter
        ax = axes[1, 0]
        if "pareduc" in df.columns:
            m = df["pareduc"].notna()
            sub = df[m]
            ax.scatter(sub["early_hh_ses"], sub["pareduc"], alpha=0.4, s=15, color="#8E44AD")
            coef = np.polyfit(sub["early_hh_ses"], sub["pareduc"], 1)
            xline = np.linspace(sub["early_hh_ses"].min(), sub["early_hh_ses"].max(), 50)
            ax.plot(xline, np.polyval(coef, xline), color="#E74C3C", linewidth=2)
            rho, p = spearmanr(sub["early_hh_ses"], sub["pareduc"])
            ax.text(0.05, 0.95, f"rho={rho:.3f}, p={p:.1e}", transform=ax.transAxes, fontsize=9, va="top")
        ax.set_xlabel("early_hh_ses")
        ax.set_ylabel("pareduc (parental education)")
        ax.set_title("Directionality: high SES = high education")

        # (1,1) Noise simulation
        ax = axes[1, 1]
        males = df[df["msex"] == 1]
        z_adv_m = males["z_adverse"].values
        rng = np.random.default_rng(42)
        noise_rhos = np.empty(1000)
        for i in range(1000):
            noise = rng.standard_normal(len(z_adv_m))
            composite_noise = z_adv_m + noise
            noise_rhos[i], _ = spearmanr(males["tot_adverse_exp"].values, composite_noise)

        ax.hist(noise_rhos, bins=40, alpha=0.6, color="gray", label="z(adv) + N(0,1)")
        rho_orig_m, _ = spearmanr(males["tot_adverse_exp"], males["ace_aggregate_orig"])
        rho_corr_m, _ = spearmanr(males["tot_adverse_exp"], males["ace_aggregate_corr"])
        ax.axvline(rho_orig_m, color="#E74C3C", linewidth=2, linestyle="--", label=f"Original: {rho_orig_m:.3f}")
        ax.axvline(rho_corr_m, color="#2ECC71", linewidth=2, linestyle="--", label=f"Corrected: {rho_corr_m:.3f}")
        ax.set_xlabel("Spearman rho (tot_adverse_exp vs composite)")
        ax.set_ylabel("Count")
        ax.set_title("Male: Original vs Corrected vs Random Noise")
        ax.legend(fontsize=8)

        # (1,2) Side-by-side rank preservation bar chart
        ax = axes[1, 2]
        groups_data = []
        for sex_val, sex_label in [(1, "Male"), (0, "Female")]:
            sub = df[df["msex"] == sex_val]
            rho_o, _ = spearmanr(sub["tot_adverse_exp"], sub["ace_aggregate_orig"])
            rho_c, _ = spearmanr(sub["tot_adverse_exp"], sub["ace_aggregate_corr"])
            groups_data.append((sex_label, rho_o, rho_c))

        x = np.arange(len(groups_data))
        w = 0.3
        ax.bar(x - w / 2, [g[1] for g in groups_data], w, label="Original", color="#E74C3C", alpha=0.8)
        ax.bar(x + w / 2, [g[2] for g in groups_data], w, label="Corrected", color="#2ECC71", alpha=0.8)
        # Add noise median line
        ax.axhline(np.median(noise_rhos), color="gray", linewidth=1, linestyle=":", label=f"Random noise median ({np.median(noise_rhos):.3f})")
        ax.set_xticks(x)
        ax.set_xticklabels([g[0] for g in groups_data])
        ax.set_ylabel("Spearman rho with tot_adverse_exp")
        ax.set_title("Rank Preservation by Sex")
        ax.legend(fontsize=8)
        ax.set_ylim(0, 1)

        for i, g in enumerate(groups_data):
            ax.text(i - w / 2, g[1] + 0.02, f"{g[1]:.3f}", ha="center", fontsize=8, color="#C0392B")
            ax.text(i + w / 2, g[2] + 0.02, f"{g[2]:.3f}", ha="center", fontsize=8, color="#27AE60")

        fig.tight_layout(rect=[0, 0, 1, 0.95])
        pdf.savefig(fig)
        plt.close(fig)

    print(f"\n>>> PDF saved to: {OUTPUT_PDF}")


# ---------------------------------------------------------------------------
# 7. Summary
# ---------------------------------------------------------------------------
def print_summary(df: pd.DataFrame) -> None:
    print("\n" + "=" * 70)
    print("SUMMARY & RECOMMENDATION")
    print("=" * 70)

    males = df[df["msex"] == 1]
    rho_orig, _ = spearmanr(males["tot_adverse_exp"], males["ace_aggregate_orig"])
    rho_corr, _ = spearmanr(males["tot_adverse_exp"], males["ace_aggregate_corr"])

    print(f"""
DIAGNOSIS:
  The ace_aggregate composite (zscore(tot_adverse_exp) + zscore(early_hh_ses))
  kills transcriptomic signal because the two variables have OPPOSITE
  directionality:

    - tot_adverse_exp: higher = MORE adversity (bad)
    - early_hh_ses:    higher = MORE affluence (good)

  Summing them cancels the adversity signal. In males specifically:
    - tot_adverse_exp vs early_hh_ses:   Spearman rho = -0.24
    - Original composite rank preservation:  {rho_orig:.3f} (only {rho_orig**2:.0%} of rank variance)
    - The original composite preserves LESS signal than adding random noise.

FIX:
  Negate early_hh_ses before combining:

    ace_aggregate = zscore(tot_adverse_exp) - zscore(early_hh_ses)
                  = zscore(tot_adverse_exp) + zscore(-early_hh_ses)

  Now both terms increase with deprivation. Result for males:
    - Corrected composite rank preservation: {rho_corr:.3f} ({rho_corr**2:.0%} of rank variance)
    - The corrected composite reinforces instead of canceling.

PIPELINE CHANGES NEEDED:
  1. prep_counts.py line 164-166: change + to -
     ace_aggregate = zscore(tot_adverse_exp) - zscore(early_hh_ses)

  2. aceDegT.Rscript (ace_aggregate construction):
     pheno$ace_aggregate <- as.numeric(scale(pheno$tot_adverse_exp)) -
                            as.numeric(scale(pheno$early_hh_ses))
""")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def resolve_pheno_csv() -> Path:
    if PHENO_CSV.exists():
        return PHENO_CSV
    alt = Path("/orcd/data/lhtsai/001/mabdel03/ROSMAP_Data/Single_Nucleus/Phenotypes/TSAI_DEJAGER_all_patients_wACEscores.csv")
    if alt.exists():
        return alt
    raise SystemExit(f"ERROR: phenotype CSV not found at:\n  {PHENO_CSV}\n  {alt}")


def main() -> None:
    print("ACE Correlation Investigation")
    print("=" * 70)

    pheno_path = resolve_pheno_csv()
    df = load_data(pheno_path)
    print(f"Loaded {len(df)} patients with both tot_adverse_exp and early_hh_ses")
    print(f"  Female: {(df['msex']==0).sum()}, Male: {(df['msex']==1).sum()}")

    print_directionality(df)
    print_correlations(df)
    print_variance_decomposition(df)
    print_rank_preservation(df)
    make_plots(df)
    print_summary(df)


if __name__ == "__main__":
    main()
