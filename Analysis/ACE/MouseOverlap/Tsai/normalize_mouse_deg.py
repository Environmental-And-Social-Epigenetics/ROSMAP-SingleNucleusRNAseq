#!/usr/bin/env python3
"""Normalize mouse LNB DEG tables and make the LNB direction explicit."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


MOUSE_FILES = {
    "PairedEnd_GeneLevel_PFC_NeuN_NGHvLNB_Dark_DESEQ2output.csv": {
        "mouse_dataset": "NeuN_Dark",
        "mouse_population": "NeuN",
        "condition": "Dark",
        "primary_human_target": "Exc",
    },
    "Paired_End_GeneLevel_PFC_NeuN_NGHvLNB_light_DESEQ2output.csv": {
        "mouse_dataset": "NeuN_Light",
        "mouse_population": "NeuN",
        "condition": "Light",
        "primary_human_target": "Exc",
    },
    "PairedEnd_PFC_pvalb_Dark_Differential.csv": {
        "mouse_dataset": "PV_Dark",
        "mouse_population": "PV",
        "condition": "Dark",
        "primary_human_target": "Inh",
    },
    "PairedEnd_PFC_pvalb_Light_Differential.csv": {
        "mouse_dataset": "PV_Light",
        "mouse_population": "PV",
        "condition": "Light",
        "primary_human_target": "Inh",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mouse-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--alpha", type=float, default=0.05)
    return parser.parse_args()


def read_mouse_file(path: Path, meta: dict[str, str], alpha: float) -> pd.DataFrame:
    df = pd.read_csv(path)
    first_col = df.columns[0]
    if first_col.startswith("Unnamed") or first_col == "":
        df = df.rename(columns={first_col: "mouse_ensembl_id"})
    elif first_col != "mouse_ensembl_id":
        df = df.rename(columns={first_col: "mouse_ensembl_id"})

    required = {"mouse_ensembl_id", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "gene.name"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")

    df = df.rename(
        columns={
            "gene.name": "mouse_symbol",
            "log2FoldChange": "mouse_original_log2FoldChange",
            "baseMean": "mouse_baseMean",
            "lfcSE": "mouse_lfcSE",
            "stat": "mouse_stat",
            "pvalue": "mouse_pvalue",
            "padj": "mouse_padj",
        }
    )
    df["mouse_lnb_log2fc"] = -pd.to_numeric(df["mouse_original_log2FoldChange"], errors="coerce")
    df["mouse_sig"] = pd.to_numeric(df["mouse_padj"], errors="coerce") < alpha
    df["mouse_direction"] = "zero"
    df.loc[df["mouse_lnb_log2fc"] > 0, "mouse_direction"] = "LNB_up"
    df.loc[df["mouse_lnb_log2fc"] < 0, "mouse_direction"] = "LNB_down"
    for key, value in meta.items():
        df[key] = value
    df["source_file"] = path.name
    return df[
        [
            "mouse_dataset",
            "mouse_population",
            "condition",
            "primary_human_target",
            "source_file",
            "mouse_ensembl_id",
            "mouse_symbol",
            "mouse_baseMean",
            "mouse_original_log2FoldChange",
            "mouse_lnb_log2fc",
            "mouse_lfcSE",
            "mouse_stat",
            "mouse_pvalue",
            "mouse_padj",
            "mouse_sig",
            "mouse_direction",
        ]
    ]


def main() -> None:
    args = parse_args()
    frames: list[pd.DataFrame] = []
    missing_files: list[str] = []

    for filename, meta in MOUSE_FILES.items():
        path = args.mouse_dir / filename
        if not path.exists():
            missing_files.append(filename)
            continue
        frames.append(read_mouse_file(path, meta, args.alpha))

    if not frames:
        missing = "\n  ".join(missing_files)
        raise SystemExit(f"No known mouse DEG files found in {args.mouse_dir}. Missing:\n  {missing}")

    out = pd.concat(frames, ignore_index=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)

    print(f"Wrote: {args.output}")
    for dataset, sub in out.groupby("mouse_dataset", sort=True):
        sig = sub[sub["mouse_sig"]]
        up = int((sig["mouse_direction"] == "LNB_up").sum())
        down = int((sig["mouse_direction"] == "LNB_down").sum())
        print(f"  {dataset}: tested={len(sub)} sig={len(sig)} LNB_up={up} LNB_down={down}")
    if missing_files:
        print("WARNING: skipped missing known files: " + ", ".join(missing_files))


if __name__ == "__main__":
    main()

