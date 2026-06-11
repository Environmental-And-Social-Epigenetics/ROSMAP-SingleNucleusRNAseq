#!/usr/bin/env python3
"""
Assemble donor_qc_report.md from the artifacts produced by donor_qc.py.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    return p.parse_args()


def md_table(df: pd.DataFrame, max_rows: int = 25) -> str:
    if len(df) > max_rows:
        df = df.head(max_rows)
    return df.to_markdown(index=True if df.index.name else False)


def safe_read(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        return None
    try:
        return pd.read_csv(path)
    except Exception:
        return None


def safe_read_indexed(path: Path, index_col: int | str = 0) -> pd.DataFrame | None:
    if not path.exists():
        return None
    try:
        return pd.read_csv(path, index_col=index_col)
    except Exception:
        return None


def section(title: str) -> str:
    return f"\n## {title}\n\n"


def fig_md(rel: str, alt: str) -> str:
    return f"![{alt}]({rel})\n\n"


def main() -> None:
    args = parse_args()
    root = args.root
    fig = root / "figures"
    tab = root / "tables"

    md: list[str] = []
    md.append("# Donor-level + annotation QC report\n\n")
    meta_path = tab / "run_metadata.json"
    if meta_path.exists():
        meta = json.loads(meta_path.read_text())
        md.append(f"_Input:_ `{meta.get('input_h5ad', '?')}`  \n")
        md.append(f"_Cells:_ {meta.get('n_cells', '?'):,}  \n")
        md.append(f"_Run:_ {meta.get('timestamp', '?')}  \n\n")

    # --- Auto verdict ---
    md.append(section("1. Summary verdict"))
    bullets: list[str] = []
    out_df = safe_read(tab / "donor_outliers.csv")
    if out_df is not None and len(out_df) > 0:
        n_multi = int((out_df["n_flags"] >= 2).sum())
        bullets.append(f"**{len(out_df)} donors flagged on ≥1 QC metric** "
                       f"({n_multi} on ≥2). See Section 2.")
    skew = safe_read(tab / "cluster_donor_skew.csv")
    if skew is not None and len(skew) > 0:
        med_inv = float(skew["inverse_simpson"].median())
        n_lowdiv = int((skew["inverse_simpson"] < med_inv / 2).sum())
        bullets.append(f"**{n_lowdiv} Leiden clusters have inverse-Simpson < median/2** "
                       f"(donor-dominated). See Section 3.")
    panel = safe_read(tab / "panel_comparison.csv")
    if panel is not None and len(panel) > 0:
        both_fail = panel[(panel["canonical_mean_frac"] < 0.30)
                          & (panel["mohammadi_mean_frac"] < 0.30)]
        only_canon = panel[(panel["canonical_mean_frac"] < 0.30)
                           & (panel["mohammadi_mean_frac"] >= 0.30)]
        bullets.append(f"**{len(both_fail)} cell types fail BOTH panels** "
                       f"(canonical and Mohammadi mean < 30%): "
                       f"{', '.join(both_fail['cell_type'].tolist()) or '—'}.")
        bullets.append(f"**{len(only_canon)} cell types pass Mohammadi but fail canonical** "
                       f"(panel-mismatch only, not necessarily mis-annotation): "
                       f"{', '.join(only_canon['cell_type'].tolist()) or '—'}.")
    sens = safe_read(tab / "drop_one_donor_sensitivity.csv")
    if sens is not None and len(sens) > 0:
        for ct in sens["cell_type"].unique():
            sub = sens[sens["cell_type"] == ct].sort_values("donors_dropped")
            if len(sub) < 2:
                continue
            baseline = float(sub.iloc[0]["cohort_mean_frac"])
            after5 = float(sub[sub["donors_dropped"] <= 5]["cohort_mean_frac"].iloc[-1])
            after10 = float(sub[sub["donors_dropped"] <= 10]["cohort_mean_frac"].iloc[-1])
            bullets.append(
                f"**{ct}** drop-one-donor: baseline {baseline:.2f}, "
                f"after dropping 5 worst donors {after5:.2f}, after 10 {after10:.2f}."
            )
    for b in bullets:
        md.append(f"- {b}\n")

    # --- Section 2: cohort QC ---
    md.append(section("2. Cohort per-donor QC"))
    md.append(fig_md("figures/cohort_qc_violins.png",
                     "Per-donor QC distributions"))
    md.append("Each point above is one donor. Outliers (MAD-based, see "
              "`tables/donor_outliers.csv`) are donors whose median sits >3 robust "
              "SDs from the cohort median on at least one metric. "
              "High MT% → stressed/dying nuclei. Low n_genes → poor library or "
              "ambient-RNA-heavy nuclei. High doublet rate → under-dissociated "
              "samples.\n\n")
    if out_df is not None and len(out_df) > 0:
        md.append("**Top flagged donors:**\n\n")
        md.append(md_table(out_df.head(20)) + "\n\n")

    # --- Section 3: cluster x donor composition ---
    md.append(section("3. Cluster × donor composition"))
    md.append(fig_md("figures/cluster_donor_heatmap.png",
                     "Cluster × donor heatmap"))
    md.append("Rows = Leiden clusters (annotated with their cell_type), cols = "
              "donors ordered worst-QC first. A healthy cluster draws cells "
              "broadly from the cohort. Single-donor-dominated clusters (bright "
              "vertical stripe) are likely capturing one donor's technical "
              "state, not a real shared cell type. Inverse-Simpson < median/2 "
              "is the cohort's bottom-half donor diversity.\n\n")
    if skew is not None and len(skew) > 0:
        worst = skew.sort_values("inverse_simpson").head(15)
        md.append("**Most donor-skewed clusters:**\n\n")
        md.append(md_table(worst) + "\n\n")

    # --- Section 4: annotation confidence ---
    md.append(section("4. Annotation confidence"))
    md.append(fig_md("figures/per_cell_ora_gap_by_celltype.png",
                     "Per-cell ORA gap by cell_type"))
    md.append("Per cell, gap = ORA score of its assigned signature minus the "
              "next-best signature. Small gap → cell sits between two cell-type "
              "identities, which is where mis-annotation piles up at cluster "
              "boundaries.\n\n")
    cluster_gap = safe_read(tab / "cluster_annotation_gap.csv")
    if cluster_gap is not None and len(cluster_gap) > 0:
        worst_clusters = cluster_gap.head(15)
        md.append("**Clusters with smallest top1−top2 gap:**\n\n")
        md.append(md_table(worst_clusters) + "\n\n")

    # --- Section 5: panel comparison (D3) ---
    md.append(section("5. Panel comparison (Mohammadi vs canonical)"))
    if panel is not None and len(panel) > 0:
        md.append("Reads as: \"of the marker genes in each panel that are "
                  "present in the assay, what fraction of cells of this "
                  "cell_type express them (>0 counts)?\" Cell types where the "
                  "Mohammadi panel passes but the canonical panel fails are "
                  "internally consistent with the pipeline; the canonical-panel "
                  "flag is a panel-mismatch artifact. Cell types where BOTH "
                  "fail are likely genuinely mis-annotated or contaminated.\n\n")
        md.append(md_table(panel.sort_values("mohammadi_mean_frac")) + "\n\n")

    # --- Section 6: per-cell-type donor drill-down (E) ---
    md.append(section("6. Per-cell-type donor drill-down"))
    md.append("**Headline question: is the problem concentrated in a few "
              "donors?** The drop-one-donor curve answers it: a sharp jump in "
              "the first few drops means yes; a slow climb means the problem "
              "is cohort-wide and won't be fixed by donor exclusion alone.\n\n")
    if sens is not None and len(sens) > 0:
        for ct in sorted(sens["cell_type"].unique()):
            ct_safe = ct.replace("/", "_")
            md.append(f"### {ct}\n\n")
            md.append(fig_md(f"figures/drop_one_donor_sensitivity_{ct_safe}.png",
                             f"{ct}: drop-one-donor sensitivity"))
            corr_pdf = fig / f"donor_qc_vs_problem_marker_{ct_safe}.pdf"
            if corr_pdf.exists():
                md.append(f"_See `figures/donor_qc_vs_problem_marker_{ct_safe}.pdf` "
                          f"for per-donor QC correlations._\n\n")
            heatmap_pdf = fig / f"donor_problem_markers_{ct_safe}_heatmap.pdf"
            if heatmap_pdf.exists():
                md.append(f"_See `figures/donor_problem_markers_{ct_safe}_heatmap.pdf` "
                          f"for the donor × failing-gene heatmap._\n\n")

    # --- Section 7: doublet drill-down ---
    md.append(section("7. Doublet drill-down"))
    md.append(fig_md("figures/doublet_score_by_celltype.png",
                     "Doublet score by cell_type"))
    md.append("scDblFinder.score is a per-sample-relative call. Cell types "
              "with elevated median score are likely doublet-contaminated, "
              "which directly explains multi-lineage marker expression (e.g. "
              "OPC expressing both PV and astrocyte markers).\n\n")
    dbl = safe_read_indexed(tab / "doublet_score_by_celltype.csv")
    if dbl is not None and len(dbl) > 0:
        md.append(md_table(dbl.sort_values("median", ascending=False)) + "\n\n")

    # --- Section 8: next steps ---
    md.append(section("8. Suggested next steps"))
    md.append("- If §6 drop-one-donor curves climb sharply with the first "
              "few drops: rerun the affected DEG analyses excluding those donors.\n")
    md.append("- If §5 shows Mohammadi-passes/canonical-fails for a cell type: "
              "the canonical-panel flag is a labelling-convention artifact; the "
              "Mohammadi-pipeline label is internally consistent.\n")
    md.append("- If §7 shows OPC/Mic with elevated doublet medians: consider "
              "tighter scDblFinder thresholds in stage 2 or a per-cell-type "
              "filter at integration.\n")
    md.append("- For clusters with small top1−top2 ORA gap (§4): consider "
              "re-clustering at finer resolution or manual review.\n")

    out_path = root / "donor_qc_report.md"
    out_path.write_text("".join(md))
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
