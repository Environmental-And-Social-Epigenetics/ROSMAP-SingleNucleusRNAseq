#!/usr/bin/env python3
"""Summarize existing inhibitory cell-proportion context for the overlap report."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prop-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def collect_results(prop_root: Path) -> pd.DataFrame:
    rows = []
    for path in sorted(prop_root.glob("results_*/*/*_results.csv")):
        rel = path.relative_to(prop_root)
        integration = rel.parts[0].replace("results_", "")
        model_group = rel.parts[1]
        df = pd.read_csv(path)
        if "cell_group" not in df.columns:
            continue
        sub = df[df["cell_group"].astype(str).str.contains(r"^(?:Inh|In-|In-PV|Exc)$", regex=True)].copy()
        if sub.empty:
            continue
        sub["integration"] = integration
        sub["model_group"] = model_group
        sub["source_file"] = str(path)
        rows.append(sub)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / "inhibitory_proportion_context.csv"
    md_path = args.output_dir / "pathology_context.md"

    if not args.prop_root.exists():
        pd.DataFrame().to_csv(csv_path, index=False)
        md_path.write_text(
            "# Pathology Context\n\n"
            "Tsai cell-type proportion outputs were not found, so no sccomp context "
            "table was generated. Human DEG overlap models already adjust for "
            "`niareagansc`. A non-AD-only inhibitory proportion analysis remains "
            "follow-up work.\n"
        )
        print(f"WARNING: proportion root not found: {args.prop_root}")
        return

    context = collect_results(args.prop_root)
    context.to_csv(csv_path, index=False)

    if context.empty:
        summary_text = "No matching inhibitory or excitatory sccomp result rows were found."
    else:
        sig = context[pd.to_numeric(context.get("c_FDR"), errors="coerce") < 0.1] if "c_FDR" in context else pd.DataFrame()
        summary_text = (
            f"Collected {len(context)} sccomp context rows from `{args.prop_root}`. "
            f"{len(sig)} rows had composition FDR < 0.1."
        )

    md_path.write_text(
        "# Pathology Context\n\n"
        "The Tsai human ACE DEG models used for overlap include `niareagansc` as a "
        "covariate, so the gene-expression overlap is pathology-adjusted at the "
        "DEG-model level.\n\n"
        f"{summary_text}\n\n"
        "Because AD can reduce inhibitory neuron density, the inhibitory proportion "
        "question should ultimately be revisited with pathology control or a "
        "non-AD subset. That subset analysis is outside this v1 overlap module.\n"
    )

    print(f"Wrote: {csv_path}")
    print(f"Wrote: {md_path}")


if __name__ == "__main__":
    main()
