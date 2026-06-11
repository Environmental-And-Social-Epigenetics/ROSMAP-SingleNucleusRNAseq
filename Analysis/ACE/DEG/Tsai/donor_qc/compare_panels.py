#!/usr/bin/env python3
"""
Side-by-side comparison of derived_batch cell-type label quality against two
marker panels:

  - Mohammadi 2020 (the panel the annotation pipeline actually used for ORA)
  - the mentor's canonical lineage panel (an independent sanity-check set)

Reads the per-cell-type fraction-expressing matrices already produced by
donor_qc.py (Mohammadi) and marker_qc/make_marker_dotplots.py (canonical) — no
H5AD reload. For each cell type, scores how well its OWN assigned signature's
markers express under each panel, reports the gap, and writes a per-gene long
table so the driver genes are visible.

Outputs (into <output-dir>, default ./panel_compare):
  tables/panel_quality_by_celltype.csv   one row per cell type, both panels
  tables/panel_per_gene_long.csv          (cell_type, panel, gene, fraction)
  figures/panel_gap_diverging.pdf/png     canonical - Mohammadi gap, ranked
  figures/panel_scatter.pdf/png           canonical vs Mohammadi own-panel mean
  panel_compare_report.md                 short narrative
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Mentor's canonical panel, grouped by lineage category.
CANONICAL_MARKERS = {
    "PV_inhibitory": ["PVALB", "GAD1", "GAD2", "SLC32A1", "SOX6", "ERBB4", "KCNS3"],
    "Microglia": ["CX3CR1", "P2RY12", "TMEM119", "AIF1", "CSF1R", "C3", "TYROBP", "TREM2", "APOE"],
    "Oligo_myelin": ["PLP1", "MBP", "MOG", "MAG", "MOBP", "CNP", "MYRF", "CLDN11"],
    "Astrocyte": ["AQP4", "GFAP", "SLC1A2", "SLC1A3", "ALDH1L1", "ATP1A2", "KCNJ10"],
    "Excitatory": ["SLC17A7", "CAMK2A", "SATB2", "CUX2", "RORB", "TBR1"],
}

# Each cell_type's expected canonical category (for the canonical "own-panel" score).
EXPECTED_CATEGORY = {
    "Ex-L2/3": "Excitatory", "Ex-L4": "Excitatory", "Ex-L4/5": "Excitatory",
    "Ex-L5": "Excitatory", "Ex-L5/6": "Excitatory", "Ex-L5/6-CC": "Excitatory",
    "Ex-NRGN": "Excitatory",
    "In-PV (Basket)": "PV_inhibitory", "In-PV (Chandelier)": "PV_inhibitory",
    "In-Rosehip": "PV_inhibitory", "In-SST": "PV_inhibitory", "In-VIP": "PV_inhibitory",
    "Mic": "Microglia",
    "Ast": "Astrocyte",
    "Oli": "Oligo_myelin", "OPC": "Oligo_myelin",
}


def parse_args() -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    p = argparse.ArgumentParser(description="Compare label quality vs Mohammadi & canonical panels.")
    p.add_argument("--mohammadi-frac",
                   type=Path,
                   default=here / "tables" / "marker_fraction_expressing_mohammadi_all.csv",
                   help="Cell_type x gene fraction-expressing matrix (Mohammadi genes).")
    p.add_argument("--mohammadi-markers",
                   type=Path,
                   default=here / "tables" / "mohammadi_markers.csv",
                   help="source,gene table for the Mohammadi panel.")
    p.add_argument("--canonical-frac",
                   type=Path,
                   default=here.parent / "marker_qc" / "tables" / "marker_fraction_expressing.csv",
                   help="Cell_type x gene fraction-expressing matrix (canonical genes).")
    p.add_argument("--output-dir", type=Path, default=here / "panel_compare")
    p.add_argument("--threshold", type=float, default=0.30,
                   help="Fraction cutoff for 'expressed' calls in the per-gene summary.")
    return p.parse_args()


def load_frac_matrix(path: Path) -> pd.DataFrame:
    """Row-indexed (cell_type) fraction matrix."""
    df = pd.read_csv(path, index_col=0)
    df.index = df.index.astype(str)
    return df


def main() -> None:
    args = parse_args()
    out = args.output_dir
    (out / "tables").mkdir(parents=True, exist_ok=True)
    (out / "figures").mkdir(parents=True, exist_ok=True)

    for pth in (args.mohammadi_frac, args.mohammadi_markers, args.canonical_frac):
        if not pth.exists():
            raise SystemExit(f"ERROR: required input not found: {pth}")

    moh_frac = load_frac_matrix(args.mohammadi_frac)        # cell_type x Mohammadi gene
    canon_frac = load_frac_matrix(args.canonical_frac)      # cell_type x canonical gene
    moh_markers = pd.read_csv(args.mohammadi_markers)       # source, gene
    # source labels in the markers CSV match cell_type labels (slash form).
    moh_by_source = {src: grp["gene"].tolist() for src, grp in moh_markers.groupby("source", sort=False)}

    canon_gene_to_cat = {g: cat for cat, genes in CANONICAL_MARKERS.items() for g in genes}

    summary_rows = []
    per_gene_rows = []

    cell_types = sorted(set(moh_frac.index) | set(canon_frac.index))
    for ct in cell_types:
        # ----- Mohammadi own-panel score -----
        moh_src_genes = moh_by_source.get(ct, [])
        moh_present = [g for g in moh_src_genes if g in moh_frac.columns]
        if ct in moh_frac.index and moh_present:
            moh_vals = moh_frac.loc[ct, moh_present].astype(float)
            moh_mean = float(moh_vals.mean())
            moh_median = float(moh_vals.median())
            moh_above = float((moh_vals > args.threshold).mean())
            for g, v in moh_vals.items():
                per_gene_rows.append({"cell_type": ct, "panel": "Mohammadi",
                                      "gene": g, "fraction_expressing": float(v)})
        else:
            moh_mean = moh_median = moh_above = np.nan

        # ----- canonical own-panel score -----
        cat = EXPECTED_CATEGORY.get(ct)
        canon_src_genes = CANONICAL_MARKERS.get(cat, []) if cat else []
        canon_present = [g for g in canon_src_genes if g in canon_frac.columns]
        if ct in canon_frac.index and canon_present:
            canon_vals = canon_frac.loc[ct, canon_present].astype(float)
            canon_mean = float(canon_vals.mean())
            canon_median = float(canon_vals.median())
            canon_above = float((canon_vals > args.threshold).mean())
            for g, v in canon_vals.items():
                per_gene_rows.append({"cell_type": ct, "panel": "Canonical",
                                      "gene": g, "fraction_expressing": float(v)})
        else:
            canon_mean = canon_median = canon_above = np.nan

        gap = (canon_mean - moh_mean) if (not np.isnan(canon_mean) and not np.isnan(moh_mean)) else np.nan

        # Verdict logic.
        thr = args.threshold
        if np.isnan(canon_mean) or np.isnan(moh_mean):
            verdict = "incomplete"
        elif moh_mean >= thr and canon_mean >= thr:
            verdict = "both_pass"
        elif moh_mean >= thr and canon_mean < thr:
            verdict = "mohammadi_only (canonical flag = panel mismatch)"
        elif moh_mean < thr and canon_mean >= thr:
            verdict = "canonical_only (unusual)"
        else:
            verdict = "both_fail (likely real annotation problem)"

        summary_rows.append({
            "cell_type": ct,
            "canonical_category": cat,
            "n_mohammadi_markers": len(moh_present),
            "n_canonical_markers": len(canon_present),
            "mohammadi_mean_frac": moh_mean,
            "mohammadi_median_frac": moh_median,
            "mohammadi_frac_markers_above_thr": moh_above,
            "canonical_mean_frac": canon_mean,
            "canonical_median_frac": canon_median,
            "canonical_frac_markers_above_thr": canon_above,
            "canonical_minus_mohammadi": gap,
            "verdict": verdict,
        })

    summary = pd.DataFrame(summary_rows).sort_values("canonical_minus_mohammadi")
    summary.to_csv(out / "tables" / "panel_quality_by_celltype.csv", index=False)
    per_gene = pd.DataFrame(per_gene_rows)
    per_gene.to_csv(out / "tables" / "panel_per_gene_long.csv", index=False)

    # ----- Diverging bar: canonical - Mohammadi gap -----
    plot_df = summary.dropna(subset=["canonical_minus_mohammadi"]).copy()
    fig, ax = plt.subplots(figsize=(9, max(4, 0.4 * len(plot_df))))
    colors = ["#d62728" if g < 0 else "#2ca02c" for g in plot_df["canonical_minus_mohammadi"]]
    ax.barh(plot_df["cell_type"], plot_df["canonical_minus_mohammadi"], color=colors)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("canonical mean frac − Mohammadi mean frac\n(< 0: canonical panel scores lower → likely panel mismatch)")
    ax.set_title("Own-panel marker expression: canonical vs Mohammadi, per cell type")
    fig.tight_layout()
    fig.savefig(out / "figures" / "panel_gap_diverging.pdf", bbox_inches="tight")
    fig.savefig(out / "figures" / "panel_gap_diverging.png", bbox_inches="tight", dpi=180)
    plt.close(fig)

    # ----- Scatter: canonical vs Mohammadi own-panel mean -----
    sc_df = summary.dropna(subset=["canonical_mean_frac", "mohammadi_mean_frac"])
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.scatter(sc_df["mohammadi_mean_frac"], sc_df["canonical_mean_frac"],
               s=40, color="steelblue")
    for _, r in sc_df.iterrows():
        ax.annotate(r["cell_type"], (r["mohammadi_mean_frac"], r["canonical_mean_frac"]),
                    fontsize=7, xytext=(3, 3), textcoords="offset points")
    lim = [0, 1]
    ax.plot(lim, lim, color="grey", linestyle="--", linewidth=0.8)
    ax.axhline(args.threshold, color="red", linestyle=":", linewidth=0.7)
    ax.axvline(args.threshold, color="red", linestyle=":", linewidth=0.7)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_xlabel("Mohammadi own-panel mean fraction expressing")
    ax.set_ylabel("Canonical own-panel mean fraction expressing")
    ax.set_title("Cell-type label quality: pipeline panel vs canonical panel\n"
                 "(red lines = threshold; below both = both panels fail)")
    fig.tight_layout()
    fig.savefig(out / "figures" / "panel_scatter.pdf", bbox_inches="tight")
    fig.savefig(out / "figures" / "panel_scatter.png", bbox_inches="tight", dpi=180)
    plt.close(fig)

    # ----- Narrative -----
    md = ["# Panel comparison: Mohammadi (pipeline) vs canonical (mentor)\n\n"]
    md.append("For each derived_batch cell-type label, we score how well that "
              "type's OWN marker set expresses, under each panel. 'Own' = the "
              "Mohammadi signature whose name matches the label, and the "
              "canonical lineage category the label belongs to.\n\n")
    md.append(f"_Threshold for 'expressed': fraction > {args.threshold}._\n\n")

    both_fail = summary[summary["verdict"].str.startswith("both_fail")]
    moh_only = summary[summary["verdict"].str.startswith("mohammadi_only")]
    both_pass = summary[summary["verdict"] == "both_pass"]

    md.append("## Verdict counts\n\n")
    md.append(f"- **both_pass**: {len(both_pass)} — "
              f"{', '.join(both_pass['cell_type']) or '—'}\n")
    md.append(f"- **mohammadi_only** (canonical flag is a panel-mismatch artifact, "
              f"not mis-annotation): {len(moh_only)} — "
              f"{', '.join(moh_only['cell_type']) or '—'}\n")
    md.append(f"- **both_fail** (likely a real annotation problem): {len(both_fail)} — "
              f"{', '.join(both_fail['cell_type']) or '—'}\n\n")

    md.append("## Per-cell-type table\n\n")
    show = summary[[
        "cell_type", "mohammadi_mean_frac", "canonical_mean_frac",
        "canonical_minus_mohammadi", "n_mohammadi_markers", "n_canonical_markers",
        "verdict",
    ]].copy()
    for c in ["mohammadi_mean_frac", "canonical_mean_frac", "canonical_minus_mohammadi"]:
        show[c] = show[c].round(3)
    md.append(show.to_markdown(index=False))
    md.append("\n\n")

    md.append("## How to read this\n\n")
    md.append("- **mohammadi_only**: the pipeline panel passes but the canonical "
              "panel fails. The canonical 'missing markers' flag is a "
              "panel-mismatch artifact — the two reference sets simply expect "
              "different genes. The label is internally consistent with the "
              "pipeline.\n")
    md.append("- **both_fail**: both panels show low own-marker expression. This "
              "is where the label is genuinely suspect.\n")
    md.append("- **canonical_minus_mohammadi** (the diverging bar): how much "
              "lower (negative) the canonical panel scores than the Mohammadi "
              "panel. Large negative = the canonical panel is the harsher / "
              "more mismatched reference for that cell type.\n\n")
    md.append("See `figures/panel_gap_diverging.png`, `figures/panel_scatter.png`, "
              "and `tables/panel_per_gene_long.csv` (per-gene drivers).\n")

    (out / "panel_compare_report.md").write_text("".join(md))

    print(f"Wrote {out / 'panel_compare_report.md'}")
    print(f"      {out / 'tables' / 'panel_quality_by_celltype.csv'}")
    print(f"      {out / 'tables' / 'panel_per_gene_long.csv'}")
    print(f"      {out / 'figures'}/panel_gap_diverging.{{pdf,png}}, panel_scatter.{{pdf,png}}")
    print()
    print(summary[["cell_type", "mohammadi_mean_frac", "canonical_mean_frac",
                   "canonical_minus_mohammadi", "verdict"]].to_string(index=False))


if __name__ == "__main__":
    main()
