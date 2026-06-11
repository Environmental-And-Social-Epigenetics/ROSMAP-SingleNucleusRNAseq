#!/usr/bin/env python3
"""
ACE SCENIC Gene Regulatory Network Analysis

Infers gene regulatory networks and transcription factor regulons using
pySCENIC (GRNBoost2 + cisTarget + AUCell), then tests regulon activity
association with ACE phenotypes via OLS regression.

This script processes ONE cell type x ONE sex combination per invocation
to enable SLURM parallelization across cell types.

Usage:
    python scenic_analysis.py \
        --cell-type Mic --sex Male \
        --input-h5ad /path/to/Mic.h5ad \
        --pheno-csv /path/to/ace_scores.csv \
        --output-dir /path/to/output/Male_Mic \
        --ranking-dir /path/to/scenic_databases \
        --tf-list /path/to/hg.txt
"""

import argparse
import logging
import os
import pickle
import sys
import warnings

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.stats import zscore
from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf

from arboreto.algo import grnboost2
from arboreto.utils import load_tf_names
from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.aucell import aucell
from pyscenic.prune import df2regulons, prune2df
from pyscenic.utils import modules_from_adjacencies

# Shared ACE phenotype loader (handles pooled-cohort duplicate projids)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "_shared"))
from load_ace_phenotype import load_for_projids  # noqa: E402

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

# Force line-buffered stdout so SLURM logs update in real time
sys.stdout = open(sys.stdout.fileno(), mode="w", buffering=1)

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

warnings.filterwarnings("ignore", category=FutureWarning)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

SEX_MAP = {"Male": 1, "Female": 0}


# ---------------------------------------------------------------------------
# Micropool aggregation
# ---------------------------------------------------------------------------

def fixed_micropool(adata, patient_col="patient_id", pool_size=50):
    """
    Aggregate cells into fixed-size micropools per patient.

    For each patient, groups of ``pool_size`` cells are summed to produce
    pseudo-bulk profiles.  Remainder cells (< pool_size) are discarded.
    Numeric metadata columns are averaged; categorical columns use mode.

    Parameters
    ----------
    adata : AnnData
        Must contain raw counts (uses .raw if available, otherwise .X).
    patient_col : str
        Column in ``adata.obs`` identifying patients.
    pool_size : int
        Number of cells per micropool.

    Returns
    -------
    AnnData
        Micropooled AnnData with summed raw counts in .X.
    """
    pooled_X = []
    pooled_obs = []
    obs_cols = [c for c in adata.obs.columns if c != patient_col]

    # Use raw counts if available
    if adata.raw is not None:
        adata = ad.AnnData(
            X=adata.raw.to_adata().X,
            obs=adata.obs.copy(),
            var=adata.raw.to_adata().var.copy(),
        )

    for patient in adata.obs[patient_col].unique():
        idx = np.where(adata.obs[patient_col] == patient)[0]
        n_cells = len(idx)
        n_pools = n_cells // pool_size

        for i in range(n_pools):
            pool_idx = idx[i * pool_size : (i + 1) * pool_size]
            pooled_counts = adata[pool_idx].X.sum(axis=0)

            if hasattr(pooled_counts, "toarray"):
                pooled_counts = pooled_counts.toarray().ravel()
            else:
                pooled_counts = np.array(pooled_counts).ravel()

            pooled_X.append(pooled_counts)

            pooled_info = {
                "patient_id": patient,
                "pool_id": f"{patient}_pool{i}",
            }
            obs_subset = adata.obs.iloc[pool_idx]
            for col in obs_cols:
                if pd.api.types.is_numeric_dtype(obs_subset[col]):
                    pooled_info[col] = obs_subset[col].mean()
                else:
                    pooled_info[col] = obs_subset[col].mode().iloc[0]

            pooled_obs.append(pooled_info)

    if len(pooled_X) == 0:
        log.warning("No micropools created (all patients had < %d cells)", pool_size)
        return None

    pooled_X = np.vstack(pooled_X)
    pooled_obs = pd.DataFrame(pooled_obs)

    adata_pooled = sc.AnnData(X=pooled_X)
    adata_pooled.var_names = adata.var_names
    adata_pooled.obs = pooled_obs.set_index("pool_id")

    return adata_pooled


# ---------------------------------------------------------------------------
# CPM normalisation
# ---------------------------------------------------------------------------

def cpm_normalize(adata):
    """
    Counts-per-million normalization on .X (in-place).

    Each micropool (row) is divided by its total counts and multiplied by 1e6.
    """
    X = adata.X if sparse.issparse(adata.X) else sparse.csr_matrix(adata.X)
    row_sums = np.array(X.sum(axis=1)).ravel()
    row_sums[row_sums == 0] = 1
    scaling = sparse.diags(1e6 / row_sums)
    adata.X = scaling.dot(X)
    return adata


# ---------------------------------------------------------------------------
# Statistical testing
# ---------------------------------------------------------------------------

def run_regulon_regression(auc_mtx, metadata, phenotype):
    """
    OLS regression of regulon AUCell scores against a phenotype.

    Model: AUCell ~ phenotype + age_death + pmi + niareagansc

    Parameters
    ----------
    auc_mtx : DataFrame
        Micropool x regulon AUCell matrix.
    metadata : DataFrame
        Per-micropool metadata with phenotype and covariates.
    phenotype : str
        Column name for the phenotype of interest.

    Returns
    -------
    DataFrame
        One row per regulon with coef, stderr, pvalue, n_targets.
    """
    results = []
    for regulon in auc_mtx.columns:
        data = pd.DataFrame({
            "AUCell": auc_mtx[regulon].values,
            "Condition": pd.to_numeric(metadata[phenotype].values, errors="coerce"),
            "Age": metadata["age_death"].values,
            "PMI": metadata["pmi"].values,
            "NIA": metadata["niareagansc"].values,
        })
        data = data.dropna()

        if len(data) < 5:
            continue

        try:
            model = smf.ols("AUCell ~ Condition + Age + PMI + NIA", data=data)
            fit = model.fit()
            term = next(
                (k for k in fit.params.keys() if "Condition" in k), None
            )
            if term is None:
                continue
            results.append({
                "regulon": regulon,
                "coef": fit.params[term],
                "stderr": fit.bse[term],
                "pvalue": fit.pvalues[term],
            })
        except Exception as exc:
            log.warning("OLS failed for regulon %s: %s", regulon, exc)
            continue

    if not results:
        return pd.DataFrame(
            columns=["regulon", "coef", "stderr", "pvalue", "padj"]
        )

    results_df = pd.DataFrame(results)
    _, padj, _, _ = multipletests(results_df["pvalue"], method="fdr_bh")
    results_df["padj"] = padj

    return results_df


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="ACE SCENIC analysis for one cell-type x sex combination"
    )
    p.add_argument("--cell-type", required=True, help="Cell type label (e.g. Mic)")
    p.add_argument("--sex", required=True, choices=["Male", "Female"],
                   help="Sex to analyse")
    p.add_argument("--phenotype", default="tot_adverse_exp",
                   help="Phenotype column (default: tot_adverse_exp)")
    p.add_argument("--integration", default="derived_batch",
                   help="Integration batch key (default: derived_batch)")
    p.add_argument("--input-h5ad", required=True,
                   help="Path to cell-type h5ad file")
    p.add_argument("--pheno-csv", required=True,
                   help="Path to ACE phenotype CSV")
    p.add_argument("--output-dir", required=True,
                   help="Output directory for this cell-type x sex")
    p.add_argument("--ranking-dir", required=True,
                   help="Directory containing SCENIC ranking databases")
    p.add_argument("--tf-list", required=True,
                   help="Path to TF names file (hg.txt)")
    p.add_argument("--pool-size", type=int, default=50,
                   help="Cells per micropool (default: 50)")
    p.add_argument("--num-workers", type=int, default=4,
                   help="Workers for motif enrichment (default: 4)")
    p.add_argument("--smoke", action="store_true",
                   help="Smoke-test mode: 500 genes, 2 workers, skip pruning")
    return p.parse_args()


def main():
    args = parse_args()

    cell_type = args.cell_type
    sex = args.sex
    phenotype = args.phenotype
    msex_val = SEX_MAP[sex]

    os.makedirs(args.output_dir, exist_ok=True)

    log.info("=== ACE SCENIC: %s / %s ===", cell_type, sex)
    log.info("Phenotype: %s", phenotype)
    log.info("Input h5ad: %s", args.input_h5ad)
    log.info("Output dir: %s", args.output_dir)
    log.info("Smoke mode: %s", args.smoke)

    # ------------------------------------------------------------------
    # 1. Load data and phenotype metadata
    # ------------------------------------------------------------------
    log.info("Loading h5ad: %s", args.input_h5ad)
    adata = ad.read_h5ad(args.input_h5ad)

    # Tsai h5ads carry `projid` (and `sample_id`); DeJager carries `patient_id`.
    # Normalize to `patient_id` so the rest of this script stays cohort-agnostic.
    if "patient_id" not in adata.obs.columns:
        for candidate in ("projid", "sample_id"):
            if candidate in adata.obs.columns:
                adata.obs["patient_id"] = adata.obs[candidate].astype(str)
                log.info("Mapped obs['%s'] -> obs['patient_id']", candidate)
                break
        else:
            raise KeyError(
                "h5ad obs has none of patient_id/projid/sample_id; "
                f"found columns: {list(adata.obs.columns)}"
            )
    adata.obs["patient_id"] = adata.obs["patient_id"].astype(str)

    log.info("Loading phenotype CSV: %s", args.pheno_csv)
    # Use the shared loader to drop duplicate projids from the pooled CSV
    # before merging — otherwise 1:many join blows up the obs table.
    cohort_projids = adata.obs["patient_id"].unique()
    pheno = load_for_projids(args.pheno_csv, cohort_projids).reset_index()

    pheno_cols = ["projid", phenotype, "msex", "age_death", "pmi", "niareagansc"]
    pheno_cols = [c for c in pheno_cols if c in pheno.columns]
    adata.obs = adata.obs.merge(
        pheno[pheno_cols],
        how="left",
        left_on="patient_id",
        right_on="projid",
    ).set_index(adata.obs.index)

    # ------------------------------------------------------------------
    # 2. Filter by sex (derive from metadata, not hardcoded lists)
    # ------------------------------------------------------------------
    sex_mask = adata.obs["msex"] == msex_val
    adata = adata[sex_mask].copy()
    n_patients = adata.obs["patient_id"].nunique()
    log.info("After sex filter (%s): %d cells from %d patients",
             sex, adata.n_obs, n_patients)

    if n_patients < 10:
        log.warning("Fewer than 10 patients (%d) -- skipping.", n_patients)
        sys.exit(0)

    # ------------------------------------------------------------------
    # 3. Micropool aggregation
    # ------------------------------------------------------------------
    log.info("Micropooling (pool_size=%d) ...", args.pool_size)
    adata = fixed_micropool(adata, patient_col="patient_id",
                            pool_size=args.pool_size)
    if adata is None or adata.n_obs == 0:
        log.warning("No micropools produced -- skipping.")
        sys.exit(0)

    log.info("Micropools: %d", adata.n_obs)

    # ------------------------------------------------------------------
    # 4. Micropool obs already carries phenotype cols (aggregated by mean
    #    during fixed_micropool). Skip re-merge; it would duplicate columns
    #    with _x/_y suffixes and break downstream column lookups.
    # ------------------------------------------------------------------
    meta = adata.obs.copy().reset_index()
    meta["patient_id"] = meta["patient_id"].astype(str)

    # Z-score continuous covariates for stable regression
    for col in ["age_death", "pmi"]:
        if col in meta.columns:
            vals = pd.to_numeric(meta[col], errors="coerce")
            mu, sd = vals.mean(), vals.std()
            if sd > 0:
                meta[col] = (vals - mu) / sd
            else:
                meta[col] = 0.0

    meta = meta.set_index("pool_id")

    # ------------------------------------------------------------------
    # 5. CPM normalize
    # ------------------------------------------------------------------
    log.info("CPM normalizing ...")
    adata = cpm_normalize(adata)

    # ------------------------------------------------------------------
    # 6. Prepare dense expression DataFrame
    # ------------------------------------------------------------------
    X = adata.X
    if sparse.issparse(X):
        X = X.toarray()
    df = pd.DataFrame(X, index=adata.obs_names, columns=adata.var_names)

    # Smoke mode: restrict to top 500 variable genes
    if args.smoke:
        gene_var = df.var(axis=0)
        top_genes = gene_var.nlargest(500).index
        df = df[top_genes]
        log.info("Smoke mode: restricted to %d genes", df.shape[1])

    log.info("Expression matrix: %d micropools x %d genes", *df.shape)

    # ------------------------------------------------------------------
    # 7. Load SCENIC ranking databases and TF list
    # ------------------------------------------------------------------
    log.info("Loading ranking databases from %s", args.ranking_dir)
    db_fnames = [
        os.path.join(
            args.ranking_dir,
            "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
        ),
        os.path.join(
            args.ranking_dir,
            "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather",
        ),
    ]
    dbs = [RankingDatabase(fname=f, name=os.path.basename(f)) for f in db_fnames]

    log.info("Loading TF names from %s", args.tf_list)
    tf_names = load_tf_names(args.tf_list)
    log.info("Loaded %d TF names", len(tf_names))

    # ------------------------------------------------------------------
    # 8. GRNBoost2 -- gene regulatory network inference
    # ------------------------------------------------------------------
    log.info("Running GRNBoost2 ...")
    adjacencies = grnboost2(expression_data=df, tf_names=tf_names, verbose=True)
    log.info("Adjacencies: %d edges", len(adjacencies))

    adj_path = os.path.join(args.output_dir, "adjacencies.csv.gz")
    adjacencies.to_csv(adj_path, index=False, compression="gzip")
    log.info("Saved adjacencies -> %s", adj_path)

    # ------------------------------------------------------------------
    # 9. Extract modules from adjacencies
    # ------------------------------------------------------------------
    log.info("Extracting modules ...")
    modules = list(modules_from_adjacencies(adjacencies, df))
    log.info("Modules: %d", len(modules))

    # ------------------------------------------------------------------
    # 10. Motif enrichment pruning (cisTarget)
    # ------------------------------------------------------------------
    num_workers = 2 if args.smoke else args.num_workers

    if args.smoke:
        log.info("Smoke mode: skipping motif pruning, using raw modules as regulons")
        # In smoke mode, convert modules directly to a minimal regulon-like
        # structure so we can still test AUCell and regression downstream.
        from pyscenic.transform import module2regulon

        regulons = []
        for m in modules[:20]:  # limit to 20 for speed
            try:
                regulons.append(m)
            except Exception:
                pass
        # AUCell can accept modules directly (they are Regulon-like)
    else:
        motif_annotation = os.path.join(
            args.ranking_dir,
            "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl",
        )
        log.info("Running motif enrichment pruning (num_workers=%d) ...", num_workers)
        enriched = prune2df(
            dbs,
            modules,
            motif_annotation,
            num_workers=num_workers,
        )
        log.info("Enriched motifs: %d rows", len(enriched))

        # Convert to regulons
        regulons = df2regulons(enriched)
        log.info("Regulons: %d", len(regulons))

    if len(regulons) == 0:
        log.warning("No regulons found -- exiting.")
        sys.exit(0)

    # Save regulons
    regulon_pkl = os.path.join(args.output_dir, "regulons.pkl")
    with open(regulon_pkl, "wb") as fh:
        pickle.dump(regulons, fh)
    log.info("Saved regulons -> %s", regulon_pkl)

    # Human-readable regulon list
    regulon_list = []
    for r in regulons:
        name = str(r)
        n_targets = len(r.genes) if hasattr(r, "genes") else len(r)
        regulon_list.append({"regulon": name, "n_targets": n_targets})
    regulon_list_df = pd.DataFrame(regulon_list)
    regulon_list_path = os.path.join(args.output_dir, "regulons_list.csv")
    regulon_list_df.to_csv(regulon_list_path, index=False)
    log.info("Saved regulon list -> %s", regulon_list_path)

    # ------------------------------------------------------------------
    # 11. AUCell scoring
    # ------------------------------------------------------------------
    log.info("Running AUCell (num_workers=%d) ...", num_workers)
    auc_mtx = aucell(df, regulons, num_workers=num_workers)
    log.info("AUCell matrix: %d micropools x %d regulons", *auc_mtx.shape)

    auc_path = os.path.join(args.output_dir, "auc_matrix.csv")
    auc_mtx.to_csv(auc_path)
    log.info("Saved AUCell matrix -> %s", auc_path)

    # ------------------------------------------------------------------
    # 12. OLS regression: regulon activity ~ phenotype + covariates
    # ------------------------------------------------------------------
    log.info("Running OLS regression (%s ~ %s + covariates) ...",
             "AUCell", phenotype)

    # Align metadata index with auc_mtx index
    meta_aligned = meta.loc[meta.index.isin(auc_mtx.index)].reindex(auc_mtx.index)

    results_df = run_regulon_regression(auc_mtx, meta_aligned, phenotype)

    # Add regulon target counts
    regulon_name_to_ntargets = {
        str(r): (len(r.genes) if hasattr(r, "genes") else len(r))
        for r in regulons
    }
    results_df["n_targets"] = results_df["regulon"].map(regulon_name_to_ntargets)

    results_path = os.path.join(args.output_dir, "regression_results.csv")
    results_df.to_csv(results_path, index=False)
    log.info("Saved regression results -> %s", results_path)

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    n_sig = (results_df["padj"] < 0.05).sum() if "padj" in results_df.columns else 0
    log.info("Significant regulons (padj < 0.05): %d / %d", n_sig, len(results_df))
    log.info("=== Done: %s / %s ===", cell_type, sex)


if __name__ == "__main__":
    main()
