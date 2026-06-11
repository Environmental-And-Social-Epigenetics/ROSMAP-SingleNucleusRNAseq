#!/usr/bin/env python3
"""
Cell-type "confusion" QC for ACE DEG hit lists.

Detects whether a canonical marker gene for cell type A turns up as a STRONG DEG
hit in a DIFFERENT cell type B -- a signal of cell-type annotation
contamination/ambiguity. Operates on the per-cell-type DEG result CSVs (which
genes ACE significantly moves), NOT on raw expression (that is marker_qc/).

Two marker sources:
  (A) Collaborator canonical panel  -- small, lineage-specific. PRIMARY evidence
      (independent of how cells were annotated).
  (B) Mohammadi 2020 marker list    -- the annotation basis itself, so cross-firing
      is partly CIRCULAR; corroborating only. Restricted to specific markers.

Outputs (under <output-dir>/tables and <output-dir>/figures):
  - confusion_events_<source>_<arm>.csv     one row per (cell_type, gene) event
  - confusion_summary_<source>_<arm>.csv    per cell type + hypergeometric enrichment
  - confusion_matrix_<source>_<arm>.csv     cell-type x marker-category counts
  - figures/confusion_heatmap_<source>_<arm>.{png,pdf}
  - run_metadata.json
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
try:
    sys.stdout.reconfigure(line_buffering=True)
except Exception:
    pass

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import hypergeom  # noqa: E402

# ---------------------------------------------------------------------------
# Reused constants. Source of truth:
#   marker_qc/make_marker_dotplots.py  (MARKERS, EXPECTED_CATEGORY, CELLTYPE_ORDER)
#   prep_celltype_splits.py            (clean_ct_name, broad_group)
# Inlined here to keep this QC script dependency-light (no scanpy import). Keep
# in sync if the marker panel changes.
# ---------------------------------------------------------------------------

# Collaborator canonical marker panel (category -> genes).
MARKERS = {
    "PV_inhibitory": ["PVALB", "GAD1", "GAD2", "SLC32A1", "SOX6", "ERBB4", "KCNS3"],
    "Microglia": ["CX3CR1", "P2RY12", "TMEM119", "AIF1", "CSF1R", "C3", "TYROBP", "TREM2", "APOE"],
    "Oligo_myelin": ["PLP1", "MBP", "MOG", "MAG", "MOBP", "CNP", "MYRF", "CLDN11"],
    "Astrocyte": ["AQP4", "GFAP", "SLC1A2", "SLC1A3", "ALDH1L1", "ATP1A2", "KCNJ10"],
    "Excitatory": ["SLC17A7", "CAMK2A", "SATB2", "CUX2", "RORB", "TBR1"],
}

# Which panel category each DEG cell-type label "should" match (cleaned form).
EXPECTED_CATEGORY = {
    "Ex-L2_3": "Excitatory",
    "Ex-L4": "Excitatory",
    "Ex-L4_5": "Excitatory",
    "Ex-L5": "Excitatory",
    "Ex-L5_6": "Excitatory",
    "Ex-L5_6-CC": "Excitatory",
    "Ex-NRGN": "Excitatory",
    "Inh": "PV_inhibitory",
    "In-PV_Basket": "PV_inhibitory",
    "In-PV_Chandelier": "PV_inhibitory",
    "In-Rosehip": "PV_inhibitory",
    "In-SST": "PV_inhibitory",
    "In-VIP": "PV_inhibitory",
    "Mic": "Microglia",
    "Ast": "Astrocyte",
    "Oli": "Oligo_myelin",
    "OPC": "Oligo_myelin",
}

# Lineage of each panel category (for cross-lineage scoring).
PANEL_CATEGORY_LINEAGE = {
    "Excitatory": "Excitatory",
    "PV_inhibitory": "Inhibitory",
    "Microglia": "Mic",
    "Astrocyte": "Ast",
    "Oligo_myelin": "Oli",  # Oli + OPC share the oligo lineage in this map
}

# Deterministic cell-type order for plots / matrix rows (cleaned labels only).
CELLTYPE_ORDER = [
    "Ex-L2_3", "Ex-L4", "Ex-L4_5", "Ex-L5", "Ex-L5_6", "Ex-L5_6-CC", "Ex-NRGN",
    "Inh", "In-PV_Basket", "In-PV_Chandelier", "In-Rosehip", "In-SST", "In-VIP",
    "Ast", "Oli", "OPC", "Mic", "Endo",
]


def clean_ct_name(ct: str) -> str:
    """Match prep_celltype_splits.clean_ct_name."""
    return ct.replace("/", "_").replace(" ", "_").replace("(", "").replace(")", "")


def lineage_of(ct_clean: str) -> str:
    """Broad lineage of a DEG cell-type label. Glial/vascular types are their own
    lineage; Oli and OPC are treated as one oligo lineage (oligodendrocyte axis)."""
    if ct_clean.startswith("Ex-"):
        return "Excitatory"
    if ct_clean.startswith("In-") or ct_clean == "Inh":
        return "Inhibitory"
    if ct_clean in ("Oli", "OPC"):
        return "Oli"
    return ct_clean  # Ast, Mic, Endo


# ---------------------------------------------------------------------------
# Arm -> result dir + per-contrast filename suffix.
# Mirrors report/summarize_ACE_AD_models.R conventions.
# ---------------------------------------------------------------------------
def arm_config(arm: str, contrast: str) -> tuple[str, str]:
    """Return (results_subdir, file_suffix) for a given arm + contrast."""
    subdir = f"results_derived_batch_{arm}_AllCellTypes"
    single_contrast_arms = {"MaleNoADadj", "MaleContAD", "MaleNiaReagan"}
    if arm in single_contrast_arms:
        suffix = arm
    else:
        # MaleBinaryAD / MaleAceByAD / MaleAncova_<v> -> "<arm>_<contrast>"
        suffix = f"{arm}_{contrast}"
    return subdir, suffix


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Cell-type confusion QC for ACE DEG hits.")
    here = Path(__file__).resolve().parent
    p.add_argument("--arm", default="MaleContAD")
    p.add_argument("--contrast", default="ACEmain")
    p.add_argument("--phenotype", default="tot_adverse_exp")
    p.add_argument("--padj", type=float, default=0.05, help="padj cutoff for a strong hit")
    p.add_argument("--lfc", type=float, default=1.0, help="|log2FC| cutoff for a strong hit")
    p.add_argument("--spec-min", type=float, default=0.5,
                   help="min Mohammadi specificity (1/#types) to count a marker; 0.5 = <=2 types")
    p.add_argument("--source", choices=["panel", "mohammadi", "both"], default="both")
    p.add_argument("--results-root", type=Path, default=None,
                   help="Defaults to ${ANALYSIS_OUTPUT_ROOT}/ACE/DEG/Tsai")
    p.add_argument("--mohammadi-csv", type=Path,
                   default=here.parent / "donor_qc" / "tables" / "mohammadi_markers.csv")
    p.add_argument("--output-dir", type=Path, default=here)
    p.add_argument("--timestamp", default="", help="ISO timestamp for run_metadata (optional)")
    return p.parse_args()


def resolve_results_root(args: argparse.Namespace) -> Path:
    if args.results_root is not None:
        return args.results_root.resolve()
    root = os.environ.get("ANALYSIS_OUTPUT_ROOT")
    if not root:
        raise SystemExit("ERROR: --results-root not given and ANALYSIS_OUTPUT_ROOT not set.")
    return (Path(root) / "ACE" / "DEG" / "Tsai").resolve()


# ---------------------------------------------------------------------------
# Marker ownership tables.
# Each returns: owner_map  (gene -> set of owner cell-type labels [cleaned]),
#               owner_lineage_map (gene -> set of owner lineages),
#               specificity (gene -> 1/#owner_types),
#               and a column-ordering for the confusion matrix.
# ---------------------------------------------------------------------------
def build_panel_tables() -> dict:
    gene_to_cat: dict[str, str] = {}
    for cat, genes in MARKERS.items():
        for g in genes:
            gene_to_cat[g.upper()] = cat
    # specificity is 1.0 for all panel genes EXCEPT those listed in >1 category.
    # (APOE is the only panel gene canonical to 2 lineages in Mohammadi, but in the
    # collaborator panel each gene is listed under exactly one category, so panel
    # specificity is 1.0 by construction.)
    cat_lineage = PANEL_CATEGORY_LINEAGE
    return {
        "gene_to_owner": {g: {gene_to_cat[g]} for g in gene_to_cat},  # owner = category
        "gene_lineage": {g: {cat_lineage[gene_to_cat[g]]} for g in gene_to_cat},
        "specificity": {g: 1.0 for g in gene_to_cat},
        "matrix_cols": list(MARKERS.keys()),
        "owner_kind": "category",
    }


def build_mohammadi_tables(csv_path: Path, spec_min: float) -> dict:
    df = pd.read_csv(csv_path)
    df["gene"] = df["gene"].astype(str).str.upper()
    df["ct"] = df["source"].astype(str).map(clean_ct_name)
    gene_to_owner: dict[str, set] = {}
    for g, ct in zip(df["gene"], df["ct"]):
        gene_to_owner.setdefault(g, set()).add(ct)
    specificity = {g: 1.0 / len(cts) for g, cts in gene_to_owner.items()}
    gene_lineage = {g: {lineage_of(ct) for ct in cts} for g, cts in gene_to_owner.items()}
    # Matrix columns = the 17 Mohammadi cell types, in CELLTYPE_ORDER where present.
    cols = [c for c in CELLTYPE_ORDER if c in set(df["ct"])]
    return {
        "gene_to_owner": gene_to_owner,
        "gene_lineage": gene_lineage,
        "specificity": specificity,
        "matrix_cols": cols,
        "owner_kind": "source_celltype",
        "spec_min": spec_min,
    }


# ---------------------------------------------------------------------------
# DEG reading.
# ---------------------------------------------------------------------------
def discover_deg_files(results_root: Path, subdir: str, suffix: str, phenotype: str):
    base = results_root / subdir / phenotype
    if not base.is_dir():
        raise SystemExit(f"ERROR: DEG dir not found: {base}")
    out = {}
    prefix = f"deseqAnalysisACE_{phenotype}_"
    tail = f"_{suffix}.csv"
    for f in sorted(base.glob(f"{prefix}*{tail}")):
        name = f.name
        # exclude filtered / summary / smoke variants
        if any(tok in name for tok in ("_padj", "smoke_summary", "filter_")):
            continue
        ct = name[len(prefix):-len(tail)]
        out[ct] = f
    return out


def read_deg(path: Path) -> pd.DataFrame | None:
    try:
        df = pd.read_csv(path)
    except Exception:
        return None
    if "gene" not in df.columns:
        return None
    df["gene"] = df["gene"].astype(str).str.upper()
    return df


# ---------------------------------------------------------------------------
# Core confusion analysis for one marker source.
# ---------------------------------------------------------------------------
def analyze_source(source_name: str, tables: dict, deg_files: dict, args: argparse.Namespace):
    gene_to_owner = tables["gene_to_owner"]
    gene_lineage = tables["gene_lineage"]
    specificity = tables["specificity"]
    spec_min = tables.get("spec_min", 0.0)
    owner_kind = tables["owner_kind"]
    matrix_cols = tables["matrix_cols"]

    # marker universe (uppercase) and the specific subset
    all_markers = set(gene_to_owner)
    specific_markers = {g for g in all_markers if specificity[g] >= spec_min}

    events = []
    summary_rows = []
    matrix = pd.DataFrame(0, index=[c for c in CELLTYPE_ORDER if c in deg_files],
                          columns=matrix_cols, dtype=int)

    for ct in [c for c in CELLTYPE_ORDER if c in deg_files] + \
              [c for c in deg_files if c not in CELLTYPE_ORDER]:
        df = read_deg(deg_files[ct])
        if df is None:
            continue
        df = df.dropna(subset=["padj"]).copy()
        df = df.sort_values(["padj", "pvalue"]).reset_index(drop=True)
        df["rank"] = np.arange(1, len(df) + 1)
        tested = set(df["gene"])
        ct_lineage = lineage_of(ct)

        # strong hits
        strong = df[(df["padj"] < args.padj) & (df["log2FoldChange"].abs() >= args.lfc)]
        n_strong = len(strong)

        # own ownership for this DEG cell type (for "cross" determination)
        if owner_kind == "category":
            own_owner = EXPECTED_CATEGORY.get(ct)  # a category string or None
            own_owners = {own_owner} if own_owner else set()
        else:
            # Mohammadi: own = this ct; for Inh (no Mohammadi key) own = all In-* types
            if ct == "Inh":
                own_owners = {c for c in matrix_cols if c.startswith("In-")}
            else:
                own_owners = {ct}

        # --- confusion matrix: count strong hits that are markers of each col ---
        for g in strong["gene"]:
            owners = gene_to_owner.get(g)
            if not owners:
                continue
            if owner_kind == "category":
                for cat in owners:
                    if cat in matrix.columns:
                        matrix.loc[ct, cat] += 1
            else:
                if specificity[g] < spec_min:
                    continue
                for oc in owners:
                    if oc in matrix.columns:
                        matrix.loc[ct, oc] += 1

        # --- events: strong hits that are cross-OWNER specific markers ---
        ct_events = []
        for _, row in strong.iterrows():
            g = row["gene"]
            owners = gene_to_owner.get(g)
            if not owners:
                continue
            if owner_kind != "category" and specificity[g] < spec_min:
                continue
            # cross-owner: this ct is not among the gene's owners
            if owner_kind == "category":
                is_cross_owner = (not own_owners) or (not owners & own_owners)
            else:
                is_cross_owner = not (owners & own_owners)
            if not is_cross_owner:
                continue
            owner_lins = gene_lineage.get(g, set())
            cross_lineage = ct_lineage not in owner_lins
            sev = specificity[g] * (2.0 if cross_lineage else 0.5)
            ev = {
                "deg_cell_type": ct,
                "deg_lineage": ct_lineage,
                "gene": g,
                "owner": "|".join(sorted(owners)),
                "owner_lineage": "|".join(sorted(owner_lins)),
                "marker_source": source_name,
                "rank": int(row["rank"]),
                "padj": float(row["padj"]),
                "log2FC": float(row["log2FoldChange"]),
                "baseMean": float(row["baseMean"]),
                "cross_lineage": bool(cross_lineage),
                "specificity": float(specificity[g]),
                "severity": float(sev),
                "hit_reason": f"padj<{args.padj}&|lfc|>={args.lfc}",
            }
            events.append(ev)
            ct_events.append(ev)

        # --- hypergeometric enrichment of cross-LINEAGE specific markers ---
        # population M = tested genes; successes K = cross-lineage specific markers
        # present in tested pool; draws n = n_strong; observed k = strong cross-lineage events.
        cross_specific_pool = {
            g for g in (specific_markers & tested)
            if ct_lineage not in gene_lineage.get(g, set())
        }
        M = len(tested)
        K = len(cross_specific_pool)
        n = n_strong
        cross_events = [e for e in ct_events if e["cross_lineage"]]
        k = len(cross_events)
        if n_strong < 5 or K == 0 or M == 0:
            expected = np.nan
            fold = np.nan
            pval = np.nan
            note = "too_few_hits" if n_strong < 5 else ("no_cross_markers_in_pool" if K == 0 else "")
        else:
            expected = n * K / M
            fold = (k / expected) if expected > 0 else np.nan
            pval = float(hypergeom.sf(k - 1, M, K, n))
            note = ""

        top_cross = sorted(cross_events, key=lambda e: e["padj"])[:10]
        summary_rows.append({
            "deg_cell_type": ct,
            "deg_lineage": ct_lineage,
            "n_top_hits": int(n_strong),
            "n_events": int(len(ct_events)),
            "n_cross_lineage": int(k),
            "frac_cross_lineage": float(k / n_strong) if n_strong else np.nan,
            "mean_severity": float(np.mean([e["severity"] for e in ct_events])) if ct_events else 0.0,
            "hyper_M": int(M),
            "hyper_K": int(K),
            "hyper_n": int(n),
            "observed_k": int(k),
            "expected_k": expected,
            "fold_enrichment": fold,
            "pvalue": pval,
            "top_cross_genes": ", ".join(f"{e['gene']}({e['owner']})" for e in top_cross),
            "note": note,
        })

    events_df = pd.DataFrame(events)
    summary_df = pd.DataFrame(summary_rows)
    # BH adjust the hypergeom p-values across cell types (ignoring NA)
    summary_df["padj"] = bh_adjust(summary_df["pvalue"].values)
    return events_df, summary_df, matrix


def bh_adjust(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    out = np.full_like(p, np.nan, dtype=float)
    mask = ~np.isnan(p)
    if mask.sum() == 0:
        return out
    pv = p[mask]
    order = np.argsort(pv)
    ranked = pv[order]
    m = len(pv)
    adj = ranked * m / (np.arange(1, m + 1))
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    adj = np.clip(adj, 0, 1)
    res = np.empty(m)
    res[order] = adj
    out[mask] = res
    return out


def draw_heatmap(matrix: pd.DataFrame, title: str, out_base: Path):
    if matrix.empty or matrix.values.sum() == 0:
        print(f"  (heatmap skipped: empty matrix for {title})")
        return
    fig_h = max(3.0, 0.42 * matrix.shape[0] + 1.5)
    fig_w = max(4.0, 0.6 * matrix.shape[1] + 2.0)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    data = matrix.values.astype(float)
    im = ax.imshow(data, aspect="auto", cmap="magma")
    ax.set_xticks(range(matrix.shape[1]))
    ax.set_xticklabels(matrix.columns, rotation=90, fontsize=7)
    ax.set_yticks(range(matrix.shape[0]))
    ax.set_yticklabels(matrix.index, fontsize=7)
    ax.set_xlabel("marker owner (category / source cell type)")
    ax.set_ylabel("DEG cell type")
    ax.set_title(title, fontsize=9)
    # annotate counts
    vmax = data.max() if data.max() > 0 else 1
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            v = int(data[i, j])
            if v > 0:
                ax.text(j, i, str(v), ha="center", va="center", fontsize=6,
                        color="white" if data[i, j] < 0.6 * vmax else "black")
    fig.colorbar(im, ax=ax, label="# strong DEG hits that are this owner's markers")
    fig.tight_layout()
    fig.savefig(out_base.with_suffix(".png"), dpi=150)
    fig.savefig(out_base.with_suffix(".pdf"))
    plt.close(fig)
    print(f"  Wrote: {out_base.name}.{{png,pdf}}")


def main() -> None:
    args = parse_args()
    results_root = resolve_results_root(args)
    subdir, suffix = arm_config(args.arm, args.contrast)
    deg_files = discover_deg_files(results_root, subdir, suffix, args.phenotype)
    if not deg_files:
        raise SystemExit(f"ERROR: no DEG files for arm={args.arm} suffix={suffix} under {results_root/subdir}")

    tables_dir = args.output_dir / "tables"
    figs_dir = args.output_dir / "figures"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figs_dir.mkdir(parents=True, exist_ok=True)

    print(f"=== Cell-type confusion QC ===")
    print(f"Arm:         {args.arm}  contrast={args.contrast}  suffix={suffix}")
    print(f"Results dir: {results_root/subdir/args.phenotype}")
    print(f"DEG files:   {len(deg_files)} cell types")
    print(f"Strong hit:  padj<{args.padj} & |log2FC|>={args.lfc}")
    print(f"Mohammadi spec-min: {args.spec_min} (1/#types; 0.5 = <=2 types)\n")

    sources = ["panel", "mohammadi"] if args.source == "both" else [args.source]
    universe_sizes = {}
    for src in sources:
        if src == "panel":
            tables = build_panel_tables()
        else:
            tables = build_mohammadi_tables(args.mohammadi_csv, args.spec_min)
        universe_sizes[src] = len(tables["gene_to_owner"])
        events_df, summary_df, matrix = analyze_source(src, tables, deg_files, args)

        tag = f"{src}_{args.arm}"
        ev_path = tables_dir / f"confusion_events_{tag}.csv"
        sm_path = tables_dir / f"confusion_summary_{tag}.csv"
        mx_path = tables_dir / f"confusion_matrix_{tag}.csv"
        events_df.to_csv(ev_path, index=False)
        summary_df.to_csv(sm_path, index=False)
        matrix.to_csv(mx_path)
        print(f"[{src}] wrote: {ev_path.name} ({len(events_df)} events), "
              f"{sm_path.name}, {mx_path.name}")

        draw_heatmap(matrix,
                     f"Cell-type confusion ({src}, {args.arm})\n"
                     f"strong hits = padj<{args.padj} & |log2FC|>={args.lfc}",
                     figs_dir / f"confusion_heatmap_{tag}")

        # console highlight: cross-lineage enriched cell types
        if not summary_df.empty:
            flagged = summary_df[(summary_df["fold_enrichment"] > 1) &
                                 (summary_df["padj"] < 0.05)]
            if len(flagged):
                print(f"  [{src}] cross-lineage-enriched cell types (fold>1, padj<0.05):")
                for _, r in flagged.sort_values("padj").iterrows():
                    print(f"    {r['deg_cell_type']:>16}  k={r['observed_k']} "
                          f"exp={r['expected_k']:.2f} fold={r['fold_enrichment']:.1f} "
                          f"padj={r['padj']:.2e}  | {r['top_cross_genes']}")
            else:
                print(f"  [{src}] no cell types with significant cross-lineage enrichment.")

    meta = {
        "arm": args.arm,
        "contrast": args.contrast,
        "suffix": suffix,
        "phenotype": args.phenotype,
        "padj_cut": args.padj,
        "lfc_cut": args.lfc,
        "spec_min": args.spec_min,
        "sources": sources,
        "results_root": str(results_root),
        "n_celltypes": len(deg_files),
        "marker_universe_sizes": universe_sizes,
        "timestamp": args.timestamp,
    }
    with open(args.output_dir / "run_metadata.json", "w") as fh:
        json.dump(meta, fh, indent=2)
    print(f"\nDone. Outputs under: {args.output_dir}")


if __name__ == "__main__":
    main()
