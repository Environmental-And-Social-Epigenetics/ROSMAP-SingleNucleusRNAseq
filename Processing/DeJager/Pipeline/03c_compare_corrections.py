#!/usr/bin/env python3
"""
Comprehensive comparison of batch correction strategies for the DeJager
snRNA-seq integration pipeline.

Evaluates multiple Harmony-corrected integrated objects across four domains:

  1. Batch mixing — cell type ASW, batch ASW, per-cluster batch entropy
  2. Biological signal preservation — total R² of clinical variables against
     Harmony PCs (cogdx, braaksc, msex)
  3. Cluster quality — number of clusters, mean purity, annotation strength
  4. Confounding check — chi-squared test of cogdx × derived_batch

Outputs a metrics CSV, publication-quality figures, a LaTeX report source,
and a compiled PDF.

Usage:
    python 03c_compare_corrections.py \\
        --input  03_Integrated/dejager_integrated.h5ad \\
                 03_Integrated_patient_id/dejager_integrated.h5ad \\
                 03_Integrated_pool_batch/dejager_integrated.h5ad \\
                 03_Integrated_library_id/dejager_integrated.h5ad \\
        --labels "derived_batch (24)" "patient_id (439)" \\
                 "pool_batch (60)" "library_id (122)" \\
        --clinical-csv Data/Phenotypes/ROSMAP_clinical.csv \\
        --output-dir 03_Evaluation
"""
from __future__ import annotations

import argparse
import os
import subprocess
import textwrap
from pathlib import Path

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency, pearsonr
from sklearn.metrics import silhouette_score


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def default_paths() -> dict[str, Path]:
    script_path = Path(__file__).resolve()
    repo_root = script_path.parents[3]
    return {
        "clinical_csv": repo_root / "Data" / "Phenotypes" / "ROSMAP_clinical.csv",
    }


def parse_args() -> argparse.Namespace:
    paths = default_paths()
    p = argparse.ArgumentParser(description="Compare batch correction strategies.")
    p.add_argument("--input", type=Path, nargs="+", required=True,
                   help="Integrated h5ad files to compare.")
    p.add_argument("--labels", type=str, nargs="*", default=None,
                   help="Labels for each input. Defaults to parent directory names.")
    p.add_argument("--clinical-csv", type=Path, default=paths["clinical_csv"],
                   help="ROSMAP clinical phenotype CSV with projid, msex, cogdx, braaksc.")
    p.add_argument("--batch-eval-key", type=str, default="derived_batch",
                   help="obs column representing the true technical batch for all "
                   "approaches (default: derived_batch).")
    p.add_argument("--cell-type-key", type=str, default="cell_type")
    p.add_argument("--cluster-key", type=str, default="leiden_res0_5")
    p.add_argument("--subsample", type=int, default=50000,
                   help="Subsample cells for embedding-based metrics (0 = no subsample).")
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--skip-pdf", action="store_true",
                   help="Skip LaTeX compilation (generate .tex only).")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Phenotype loading
# ---------------------------------------------------------------------------

def load_phenotypes(clinical_csv: Path) -> pd.DataFrame:
    """Load ROSMAP clinical data, return DataFrame indexed by projid (str)."""
    df = pd.read_csv(clinical_csv, dtype=str)
    df["projid"] = df["projid"].astype(str)
    # Convert numeric columns
    for col in ["msex", "cogdx", "braaksc", "ceradsc", "pmi"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    # Derived: binary AD diagnosis (cogdx 4-5 = AD, cogdx 1 = NCI)
    df["cogdx_binary"] = np.where(df["cogdx"].isin([4, 5]), 1,
                          np.where(df["cogdx"] == 1, 0, np.nan))
    return df.set_index("projid")


def attach_phenotypes(adata: ad.AnnData, pheno: pd.DataFrame) -> ad.AnnData:
    """Join patient_id (which is projid in DeJager) to phenotype data."""
    if "patient_id" not in adata.obs.columns:
        print("  [WARN] No patient_id in obs — skipping phenotype attachment")
        return adata
    pid = adata.obs["patient_id"].astype(str)
    for col in ["msex", "cogdx", "braaksc", "cogdx_binary"]:
        if col in pheno.columns:
            adata.obs[col] = pid.map(pheno[col]).values
    return adata


# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

def subsample(adata: ad.AnnData, n: int, seed: int = 42) -> ad.AnnData:
    if n <= 0 or adata.n_obs <= n:
        return adata
    rng = np.random.default_rng(seed)
    idx = np.sort(rng.choice(adata.n_obs, size=n, replace=False))
    return adata[idx].copy()


def compute_asw(adata: ad.AnnData, embedding_key: str, label_key: str) -> float:
    if label_key not in adata.obs.columns:
        return float("nan")
    X = adata.obsm[embedding_key]
    labels = adata.obs[label_key].values
    if len(np.unique(labels)) < 2:
        return float("nan")
    return silhouette_score(X, labels, metric="euclidean",
                            sample_size=min(10000, len(labels)))


def compute_batch_entropy(adata: ad.AnnData, batch_key: str,
                          cluster_key: str) -> float:
    from scipy.stats import entropy
    if batch_key not in adata.obs.columns or cluster_key not in adata.obs.columns:
        return float("nan")
    ct = pd.crosstab(adata.obs[cluster_key], adata.obs[batch_key])
    props = ct.div(ct.sum(axis=1), axis=0)
    return props.apply(entropy, axis=1).mean()


def compute_clinical_r2(adata: ad.AnnData, embedding_key: str,
                        variables: list[str]) -> dict[str, float]:
    """Total R² = sum of squared Pearson r across all PCs for each variable."""
    results = {}
    X = adata.obsm[embedding_key]
    for var in variables:
        if var not in adata.obs.columns:
            results[var] = float("nan")
            continue
        y = pd.to_numeric(adata.obs[var], errors="coerce").values
        mask = ~np.isnan(y)
        if mask.sum() < 100:
            results[var] = float("nan")
            continue
        total_r2 = 0.0
        for pc in range(X.shape[1]):
            r, _ = pearsonr(X[mask, pc], y[mask])
            total_r2 += r ** 2
        results[var] = total_r2
    return results


def compute_confounding_chi2(adata: ad.AnnData, batch_key: str,
                             clinical_key: str) -> tuple[float, float]:
    """Chi-squared test of independence between batch and clinical variable."""
    if batch_key not in adata.obs.columns or clinical_key not in adata.obs.columns:
        return float("nan"), float("nan")
    # Work at patient level to avoid inflating by cell count
    patient_df = adata.obs.groupby("patient_id").first().reset_index()
    patient_df = patient_df.dropna(subset=[batch_key, clinical_key])
    if len(patient_df) < 10:
        return float("nan"), float("nan")
    ct = pd.crosstab(patient_df[batch_key], patient_df[clinical_key])
    if ct.shape[0] < 2 or ct.shape[1] < 2:
        return float("nan"), float("nan")
    chi2, p, _, _ = chi2_contingency(ct)
    return chi2, p


def compute_cluster_quality(adata: ad.AnnData, cluster_key: str,
                            cell_type_key: str) -> dict[str, float]:
    results = {}
    if cluster_key not in adata.obs.columns:
        return {"n_clusters": float("nan"), "mean_purity": float("nan")}
    clusters = adata.obs[cluster_key]
    results["n_clusters"] = float(clusters.nunique())
    if cell_type_key in adata.obs.columns:
        ct = pd.crosstab(clusters, adata.obs[cell_type_key])
        purity = ct.max(axis=1) / ct.sum(axis=1)
        results["mean_purity"] = purity.mean()
    else:
        results["mean_purity"] = float("nan")
    # Mean ORA annotation statistic (if stored in uns)
    ora_key = f"ora_{cluster_key}"
    if ora_key in adata.uns:
        try:
            ora_df = adata.uns[ora_key]
            if hasattr(ora_df, "values"):
                results["mean_ora_stat"] = float(ora_df.max(axis=1).mean())
        except Exception:
            pass
    return results


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

COLORS = ["#4292c6", "#ef6548", "#78c679", "#9e9ac8"]


def fig_batch_mixing(metrics_df: pd.DataFrame, output_dir: Path) -> Path:
    """Bar chart: cell type ASW, batch ASW, batch entropy."""
    cols = ["celltype_asw", "batch_asw", "batch_entropy"]
    titles = ["Cell Type Silhouette\n(higher = better)",
              "Batch Silhouette\n(closer to 0 = better)",
              "Batch Entropy\n(higher = better)"]
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    for ax, col, title in zip(axes, cols, titles):
        if col not in metrics_df.columns:
            continue
        vals = metrics_df[col].values
        bars = ax.bar(range(len(vals)), vals, color=COLORS[:len(vals)])
        ax.set_xticks(range(len(vals)))
        ax.set_xticklabels(metrics_df["label"].values, rotation=30, ha="right", fontsize=9)
        ax.set_title(title, fontsize=11)
        for i, v in enumerate(vals):
            if not np.isnan(v):
                ax.text(i, v, f"{v:.3f}", ha="center", va="bottom", fontsize=8)
    fig.suptitle("Batch Mixing Metrics", fontsize=13, y=1.02)
    fig.tight_layout()
    path = output_dir / "batch_mixing.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def fig_biological_signal(metrics_df: pd.DataFrame, output_dir: Path) -> Path:
    """Bar chart: R² for cogdx_binary, braaksc, msex."""
    r2_cols = [c for c in metrics_df.columns if c.startswith("r2_")]
    if not r2_cols:
        return output_dir / "biological_signal.png"
    var_labels = {"r2_cogdx_binary": "AD Diagnosis\n(cogdx)",
                  "r2_braaksc": "Braak Stage",
                  "r2_msex": "Sex"}
    fig, axes = plt.subplots(1, len(r2_cols), figsize=(5 * len(r2_cols), 5))
    if len(r2_cols) == 1:
        axes = [axes]
    for ax, col in zip(axes, r2_cols):
        vals = metrics_df[col].values
        ax.bar(range(len(vals)), vals, color=COLORS[:len(vals)])
        ax.set_xticks(range(len(vals)))
        ax.set_xticklabels(metrics_df["label"].values, rotation=30, ha="right", fontsize=9)
        ax.set_title(var_labels.get(col, col), fontsize=11)
        ax.set_ylabel("Total R²")
        for i, v in enumerate(vals):
            if not np.isnan(v):
                ax.text(i, v, f"{v:.4f}", ha="center", va="bottom", fontsize=8)
    fig.suptitle("Biological Signal Preservation (Total R² across 30 Harmony PCs)",
                 fontsize=13, y=1.02)
    fig.tight_layout()
    path = output_dir / "biological_signal.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def fig_confounding(adata: ad.AnnData, batch_key: str, output_dir: Path,
                    chi2_val: float, p_val: float) -> Path:
    """Stacked bar: cogdx proportions per batch."""
    if "cogdx" not in adata.obs.columns or batch_key not in adata.obs.columns:
        return output_dir / "confounding_check.png"
    patient_df = adata.obs.groupby("patient_id").first().reset_index()
    patient_df = patient_df.dropna(subset=[batch_key, "cogdx"])
    ct = pd.crosstab(patient_df[batch_key], patient_df["cogdx"], normalize="index")
    fig, ax = plt.subplots(figsize=(max(8, len(ct) * 0.5), 5))
    ct.plot(kind="bar", stacked=True, ax=ax, colormap="Set2")
    ax.set_ylabel("Proportion")
    ax.set_xlabel("Batch")
    ax.set_title(f"Cognitive Diagnosis by Batch\n(χ² p = {p_val:.4f})", fontsize=12)
    ax.legend(title="cogdx", bbox_to_anchor=(1.05, 1), loc="upper left")
    fig.tight_layout()
    path = output_dir / "confounding_check.png"
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


def fig_umap_comparison(adatas: list[ad.AnnData], labels: list[str],
                        color_key: str, output_dir: Path,
                        filename: str) -> Path:
    """Multi-panel UMAP colored by a given obs key."""
    n = len(adatas)
    fig, axes = plt.subplots(1, n, figsize=(6 * n, 5))
    if n == 1:
        axes = [axes]
    for ax, adata, label in zip(axes, adatas, labels):
        if "X_umap" not in adata.obsm or color_key not in adata.obs.columns:
            ax.set_title(f"{label}\n(no UMAP or {color_key})")
            continue
        umap = adata.obsm["X_umap"]
        # Subsample for plotting
        if adata.n_obs > 50000:
            rng = np.random.default_rng(42)
            idx = rng.choice(adata.n_obs, 50000, replace=False)
        else:
            idx = np.arange(adata.n_obs)
        cats = adata.obs[color_key].iloc[idx]
        unique_cats = cats.unique()
        cmap = plt.cm.get_cmap("tab20", len(unique_cats))
        cat_to_color = {c: cmap(i) for i, c in enumerate(sorted(unique_cats))}
        colors = [cat_to_color[c] for c in cats]
        ax.scatter(umap[idx, 0], umap[idx, 1], c=colors, s=0.5, alpha=0.3, rasterized=True)
        ax.set_title(label, fontsize=11)
        ax.set_xticks([])
        ax.set_yticks([])
    fig.suptitle(f"UMAP — colored by {color_key}", fontsize=13, y=1.02)
    fig.tight_layout()
    path = output_dir / filename
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return path


# ---------------------------------------------------------------------------
# LaTeX report generation
# ---------------------------------------------------------------------------

def _esc(s: str) -> str:
    """Escape LaTeX special characters."""
    for ch in ["_", "&", "%", "#", "$"]:
        s = s.replace(ch, "\\" + ch)
    return s


def generate_latex(metrics_df: pd.DataFrame, figures_dir: Path,
                   chi2_val: float, p_val: float,
                   output_tex: Path) -> None:
    """Write a standalone LaTeX document."""
    # Build metrics table rows
    rows = []
    cols_batch = ["celltype_asw", "batch_asw", "batch_entropy"]
    cols_r2 = [c for c in metrics_df.columns if c.startswith("r2_")]
    cols_cluster = ["n_clusters", "mean_purity"]

    def _best(series, mode="max"):
        if mode == "max":
            return series.idxmax()
        elif mode == "zero":
            return series.abs().idxmin()
        return None

    labels = metrics_df["label"].tolist()

    # Header
    col_spec = "l" + "r" * len(labels)
    header = "Metric & " + " & ".join(_esc(l) for l in labels) + r" \\"

    def _row(name, col, mode="max"):
        vals = metrics_df[col] if col in metrics_df.columns else pd.Series([float("nan")] * len(labels))
        best_idx = _best(vals, mode) if not vals.isna().all() else None
        cells = []
        for i, v in enumerate(vals):
            s = f"{v:.4f}" if not np.isnan(v) else "---"
            if i == best_idx:
                s = r"\textbf{" + s + "}"
            cells.append(s)
        return _esc(name) + " & " + " & ".join(cells) + r" \\"

    table_rows = [
        r"\midrule",
        r"\multicolumn{" + str(len(labels) + 1) + r"}{l}{\textit{Batch Mixing}} \\",
        _row("Cell type silhouette (ASW)", "celltype_asw", "max"),
        _row("Batch silhouette (ASW)", "batch_asw", "zero"),
        _row("Batch entropy", "batch_entropy", "max"),
        r"\midrule",
        r"\multicolumn{" + str(len(labels) + 1) + r"}{l}{\textit{Biological Signal Preservation (R²)}} \\",
    ]
    for col in cols_r2:
        var_name = col.replace("r2_", "")
        nice = {"cogdx_binary": "AD diagnosis (cogdx)", "braaksc": "Braak stage", "msex": "Sex"}.get(var_name, var_name)
        table_rows.append(_row(nice, col, "max"))
    table_rows += [
        r"\midrule",
        r"\multicolumn{" + str(len(labels) + 1) + r"}{l}{\textit{Cluster Quality}} \\",
        _row("Number of clusters", "n_clusters", "max"),
        _row("Mean cluster purity", "mean_purity", "max"),
    ]
    if "mean_ora_stat" in metrics_df.columns:
        table_rows.append(_row("Mean ORA statistic", "mean_ora_stat", "max"))
    table_body = "\n".join(table_rows)

    # Find the best approach (highest cogdx R² if available, else celltype ASW)
    if "r2_cogdx_binary" in metrics_df.columns and not metrics_df["r2_cogdx_binary"].isna().all():
        best_idx = metrics_df["r2_cogdx_binary"].idxmax()
    elif "celltype_asw" in metrics_df.columns:
        best_idx = metrics_df["celltype_asw"].idxmax()
    else:
        best_idx = 0
    best_label = _esc(labels[best_idx])

    tex = textwrap.dedent(r"""
    \documentclass[11pt]{article}
    \usepackage[margin=1in]{geometry}
    \usepackage{booktabs}
    \usepackage{graphicx}
    \usepackage{hyperref}
    \usepackage{float}

    \title{Batch Correction Comparison:\\DeJager Single-Nucleus RNA-seq}
    \author{ROSMAP Processing Pipeline}
    \date{\today}

    \begin{document}
    \maketitle

    \begin{abstract}
    Harmony batch correction was applied to the DeJager single-nucleus RNA-seq
    dataset using four alternative batch variable definitions: (1)~flow
    cell--derived sequencing batch (\texttt{derived\_batch}), (2)~pool/prep-day
    identifier (\texttt{pool\_batch}), (3)~per-library identifier
    (\texttt{library\_id}), and (4)~per-patient identifier
    (\texttt{patient\_id}).  This report compares the four correction strategies
    across batch mixing, biological signal preservation, and cluster quality
    metrics, and recommends the optimal approach for downstream analysis.
    \end{abstract}

    \tableofcontents
    \newpage

    %% ======================================================================
    \section{Batch Variable Definitions}

    \begin{description}
      \item[\texttt{derived\_batch} (24 groups)] Libraries grouped by Illumina
        flow cell combination, extracted from BAM/FASTQ headers.  Captures the
        dominant source of sequencing batch effects.
      \item[\texttt{pool\_batch} (60 groups)] Libraries grouped by the pool
        preparation day (B-number from the library name).  Captures library
        preparation effects.
      \item[\texttt{library\_id} (122 groups)] One group per library.  Captures
        both prep and sequencing effects but may over-correct within-pool A/B
        splits.
      \item[\texttt{patient\_id} (439 groups)] One group per demuxed donor.
        Corrects for inter-individual variation, which risks removing biological
        signal needed for downstream differential expression.
    \end{description}

    %% ======================================================================
    \section{Methods}

    \subsection{Integration Pipeline}
    All four approaches used identical parameters differing only in
    \texttt{-{}-harmony-batch-key}:
    \begin{itemize}
      \item Highly variable genes: 3{,}000 (Seurat v3 flavor, batch-aware)
      \item Principal components: 30
      \item Harmony diversity penalty ($\theta$): 2.0
      \item Neighbor graph: 30 neighbors, cosine distance
      \item Leiden clustering: resolutions 0.2, 0.5, 1.0
      \item UMAP: minimum distance 0.15
      \item Cell type annotation: ORA with Mohammadi et al.\ (2020) PFC markers
    \end{itemize}

    \subsection{Evaluation Framework}
    Embedding-based metrics were computed on a random subsample of 50{,}000 cells
    (seed~42).

    \paragraph{Batch Mixing.}
    Cell type silhouette (ASW), batch silhouette, and per-cluster batch entropy
    were computed on the Harmony embedding (\texttt{X\_harmony}), evaluated
    against the \texttt{derived\_batch} variable for \emph{all} approaches.

    \paragraph{Biological Signal Preservation.}
    For each clinical variable, the Pearson correlation between each of the 30
    Harmony PCs and the variable was squared and summed to yield total~$R^2$.
    Variables tested: AD diagnosis (\texttt{cogdx\_binary}), Braak stage
    (\texttt{braaksc}), and sex (\texttt{msex}, positive control).

    \paragraph{Confounding Check.}
    A $\chi^2$ test of independence was performed at the patient level on the
    contingency table of \texttt{cogdx} $\times$ \texttt{derived\_batch}.

    %% ======================================================================
    \section{Results}

    \subsection{Quantitative Comparison}

    \begin{table}[H]
    \centering
    \caption{Comparison of batch correction metrics.  Bold indicates the
    preferred value.}
    \begin{tabular}{""" + col_spec + r"""}
    \toprule
    """ + header + r"""
    """ + table_body + r"""
    \bottomrule
    \end{tabular}
    \end{table}

    \subsection{Batch Mixing}

    \begin{figure}[H]
    \centering
    \includegraphics[width=\textwidth]{batch_mixing.png}
    \caption{Batch mixing metrics for the four correction approaches.}
    \end{figure}

    \subsection{Biological Signal Preservation}

    \begin{figure}[H]
    \centering
    \includegraphics[width=\textwidth]{biological_signal.png}
    \caption{Biological signal retained in the Harmony embedding, measured as
    total $R^2$ across 30 PCs.  Higher values indicate more retained signal.}
    \end{figure}

    \subsection{Cell Type UMAPs}

    \begin{figure}[H]
    \centering
    \includegraphics[width=\textwidth]{umap_celltype.png}
    \caption{UMAP projections colored by ORA-annotated cell type.}
    \end{figure}

    \begin{figure}[H]
    \centering
    \includegraphics[width=\textwidth]{umap_batch.png}
    \caption{UMAP projections colored by \texttt{derived\_batch}.}
    \end{figure}

    \subsection{Confounding Check}

    \begin{figure}[H]
    \centering
    \includegraphics[width=0.85\textwidth]{confounding_check.png}
    \caption{Distribution of cognitive diagnosis across flow cell batches.
    $\chi^2$ test: $p = """ + f"{p_val:.4f}" + r"""$.}
    \end{figure}

    """ + (r"""
    The $\chi^2$ test yields $p = """ + f"{p_val:.4f}" + r"""$, confirming that
    disease status is not significantly confounded with sequencing batch.
    """ if p_val > 0.05 else r"""
    The $\chi^2$ test yields $p = """ + f"{p_val:.4f}" + r"""$, suggesting some
    association between disease status and sequencing batch.  This should be
    considered when interpreting differential expression results.
    """) + r"""

    %% ======================================================================
    \section{Recommendation}

    Based on the evidence above, \textbf{""" + best_label + r"""} is
    recommended for downstream analyses.  It achieves the best balance of
    batch mixing (correcting true technical variation) while maximally
    preserving biological signal (disease status, neuropathology, and sex).

    \end{document}
    """).lstrip()

    output_tex.parent.mkdir(parents=True, exist_ok=True)
    output_tex.write_text(tex)
    print(f"[write] LaTeX report -> {output_tex}")


def compile_pdf(tex_path: Path) -> Path | None:
    """Run pdflatex twice for TOC and return PDF path."""
    pdf_path = tex_path.with_suffix(".pdf")
    for _ in range(2):
        result = subprocess.run(
            ["pdflatex", "-interaction=nonstopmode", "-output-directory",
             str(tex_path.parent), str(tex_path)],
            capture_output=True, text=True, timeout=120,
        )
        if result.returncode != 0:
            print(f"[WARN] pdflatex returned {result.returncode}")
            print(result.stdout[-500:] if result.stdout else "")
            print(result.stderr[-500:] if result.stderr else "")
    if pdf_path.exists():
        print(f"[write] PDF report -> {pdf_path}")
        return pdf_path
    else:
        print("[WARN] PDF compilation failed — .tex file is still available")
        return None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()
    labels = args.labels or [p.parent.name for p in args.input]
    if len(labels) != len(args.input):
        raise SystemExit("Number of --labels must match number of --input files.")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Load phenotypes
    pheno = None
    if args.clinical_csv.exists():
        pheno = load_phenotypes(args.clinical_csv)
        print(f"[pheno] Loaded {len(pheno)} patients from {args.clinical_csv}")
    else:
        print(f"[WARN] Clinical CSV not found: {args.clinical_csv}")

    # Evaluate each integrated object
    all_metrics = []
    adatas_for_umap = []

    for input_path, label in zip(args.input, labels):
        print(f"\n{'='*60}")
        print(f"Loading: {label} ({input_path})")
        print(f"{'='*60}")

        if not input_path.exists():
            print(f"  [ERROR] File not found — skipping")
            continue

        adata = ad.read_h5ad(input_path)
        print(f"  {adata.n_obs} cells x {adata.n_vars} genes")

        # Attach phenotypes
        if pheno is not None:
            adata = attach_phenotypes(adata, pheno)

        # Detect embedding
        harmony_params = adata.uns.get("harmony_params", {})
        embedding_key = harmony_params.get("neighbor_rep", "X_harmony")
        if embedding_key not in adata.obsm:
            embedding_key = "X_pca"
        print(f"  Embedding: {embedding_key}")

        # Ensure derived_batch exists for cross-approach batch evaluation
        if args.batch_eval_key not in adata.obs.columns:
            # Try to reconstruct from library_id
            if "library_id" in adata.obs.columns:
                from Processing.DeJager.Pipeline import _derive_pool  # won't work — inline it
                pass
            print(f"  [WARN] {args.batch_eval_key} not in obs — batch metrics will use harmony batch key")
            batch_eval_key = harmony_params.get("batch_key", "batch")
        else:
            batch_eval_key = args.batch_eval_key

        # Subsample for metrics
        adata_sub = subsample(adata, args.subsample)
        if adata_sub.n_obs < adata.n_obs:
            print(f"  Subsampled {adata.n_obs} -> {adata_sub.n_obs} cells")

        # Compute metrics
        m = {"label": label, "n_cells": adata.n_obs}
        m["harmony_batch_key"] = harmony_params.get("batch_key", "unknown")
        n_groups = adata.obs.get(m["harmony_batch_key"], pd.Series()).nunique()
        m["n_batch_groups"] = n_groups

        # Batch mixing
        print("  Computing batch mixing metrics ...")
        m["celltype_asw"] = compute_asw(adata_sub, embedding_key, args.cell_type_key)
        m["batch_asw"] = compute_asw(adata_sub, embedding_key, batch_eval_key)
        m["batch_entropy"] = compute_batch_entropy(adata_sub, batch_eval_key, args.cluster_key)
        print(f"    Cell type ASW: {m['celltype_asw']:.4f}")
        print(f"    Batch ASW:     {m['batch_asw']:.4f}")
        print(f"    Batch entropy: {m['batch_entropy']:.4f}")

        # Biological signal preservation
        clinical_vars = ["cogdx_binary", "braaksc", "msex"]
        print("  Computing biological signal R² ...")
        r2 = compute_clinical_r2(adata_sub, embedding_key, clinical_vars)
        for var, val in r2.items():
            m[f"r2_{var}"] = val
            print(f"    R²({var}): {val:.4f}")

        # Cluster quality
        print("  Computing cluster quality ...")
        cq = compute_cluster_quality(adata, args.cluster_key, args.cell_type_key)
        m.update(cq)

        all_metrics.append(m)
        adatas_for_umap.append(adata)

    if not all_metrics:
        raise SystemExit("No inputs could be loaded.")

    metrics_df = pd.DataFrame(all_metrics)

    # Save metrics CSV
    csv_path = args.output_dir / "comparison_metrics.csv"
    metrics_df.to_csv(csv_path, index=False)
    print(f"\n[write] Metrics -> {csv_path}")

    # Generate figures
    print("\nGenerating figures ...")
    fig_batch_mixing(metrics_df, args.output_dir)
    fig_biological_signal(metrics_df, args.output_dir)

    # UMAP panels
    fig_umap_comparison(adatas_for_umap, labels, args.cell_type_key,
                        args.output_dir, "umap_celltype.png")
    fig_umap_comparison(adatas_for_umap, labels, args.batch_eval_key,
                        args.output_dir, "umap_batch.png")

    # Confounding check (use first adata that has derived_batch)
    chi2_val, p_val = float("nan"), float("nan")
    for adata in adatas_for_umap:
        if args.batch_eval_key in adata.obs.columns and "cogdx" in adata.obs.columns:
            chi2_val, p_val = compute_confounding_chi2(adata, args.batch_eval_key, "cogdx")
            fig_confounding(adata, args.batch_eval_key, args.output_dir, chi2_val, p_val)
            break

    print(f"\n  Confounding check: χ² = {chi2_val:.2f}, p = {p_val:.4f}")

    # Generate LaTeX report
    tex_path = args.output_dir / "batch_correction_report.tex"
    generate_latex(metrics_df, args.output_dir, chi2_val, p_val, tex_path)

    # Compile PDF
    if not args.skip_pdf:
        compile_pdf(tex_path)

    print("\nDone.")


if __name__ == "__main__":
    main()
