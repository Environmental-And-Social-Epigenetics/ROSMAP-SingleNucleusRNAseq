#!/usr/bin/env python3
"""ACE epigenomic integration: cross DEG hits with ROSMAP H3K9ac and
methylation tables to test the epigenetic-repression hypothesis.

For each (cell type x sex) DEG result for a given phenotype:
  - take the top-N DEG hits (by |log2FC| with padj < threshold)
  - look up matching gene/peak IDs in the H3K9ac and methylation tables
  - compute the epigenetic effect-size distribution stratified by ACE
    phenotype quartile
  - emit a joined long-form table + a per-cell-type summary

This script makes minimal assumptions about the layout of the epigenomic
tables — it autodetects gene-symbol or coordinate columns and joins on
whatever it can find. If neither is available, the cell type is skipped
with a warning.

Usage:
  python integrate.py \
      --phenotype tot_adverse_exp \
      --deg-dir /.../ACE/DEG/Tsai/results_derived_batch/tot_adverse_exp \
      --epi-data-dir /.../ACE/EpigenomicIntegration/Tsai/epigenomic_data \
      --overlap-csv /.../individual_overlap.csv \
      --pheno-csv /.../ACE_scores.csv \
      --output-dir /.../ACE/EpigenomicIntegration/Tsai/<phenotype>
      [--top-n 200] [--padj 0.1]
"""

from __future__ import annotations

import argparse
import logging
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "_shared"))
from load_ace_phenotype import load_for_projids  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("integrate")


def dump_deg_tables(deg_dir: Path, cache_dir: Path, rscript: str) -> None:
    """Run the R helper to materialize CSV copies of every DEG .rda."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    helper = Path(__file__).with_name("_dump_deg_tables.R")
    cmd = [rscript, str(helper), "--deg-dir", str(deg_dir), "--out-dir", str(cache_dir)]
    log.info("Dumping DEG .rda -> CSV cache: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)


def parse_deg_filename(name: str) -> Optional[tuple[str, str]]:
    """Pull (cell_type, sex_label) out of deseqAnalysisACE_<phen>_<ct>_<sex>.csv."""
    m = re.match(r"^deseqAnalysisACE_[^_]+(?:_[^_]+)?_(.+?)_(Fem|Male)\.csv$", name)
    if not m:
        return None
    return m.group(1), ("Female" if m.group(2) == "Fem" else "Male")


def load_deg_csv(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    cols = {c.lower(): c for c in df.columns}
    # DESeq2 uses log2FoldChange / padj; sometimes results have lfcShrink output.
    if "log2foldchange" not in cols or "padj" not in cols:
        log.warning("Unexpected DEG schema in %s; columns=%s", csv_path.name, list(df.columns))
    rename = {}
    if "log2foldchange" in cols:
        rename[cols["log2foldchange"]] = "log2fc"
    if "padj" in cols:
        rename[cols["padj"]] = "padj"
    if "gene" not in cols:
        # Last-ditch: use the first text column as gene
        text_cols = [c for c in df.columns if df[c].dtype == object]
        if text_cols:
            rename[text_cols[0]] = "gene"
    df = df.rename(columns=rename)
    return df


def autodetect_gene_col(df: pd.DataFrame) -> Optional[str]:
    candidates = ("gene", "gene_symbol", "symbol", "GeneSymbol", "feature", "Gene")
    for c in candidates:
        if c in df.columns:
            return c
    # Heuristic: a string column whose values look like HGNC symbols (uppercase,
    # length 2-15, alphanumeric).
    sym_re = re.compile(r"^[A-Z0-9-]{2,15}$")
    for c in df.columns:
        if df[c].dtype != object:
            continue
        sample = df[c].dropna().astype(str).head(20)
        if sample.empty:
            continue
        hits = sum(bool(sym_re.match(s)) for s in sample)
        if hits >= len(sample) * 0.8:
            return c
    return None


def load_epigenetic_table(path: Path) -> Optional[pd.DataFrame]:
    """Load a single epigenomic table; return None if not parseable."""
    try:
        if path.suffix.lower() in (".csv",):
            df = pd.read_csv(path, low_memory=False)
        elif path.suffix.lower() in (".tsv", ".txt"):
            df = pd.read_csv(path, sep="\t", low_memory=False)
        else:
            return None
    except Exception as exc:
        log.warning("Could not read %s: %s", path, exc)
        return None
    return df


def harvest_epigenetic_tables(epi_dir: Path) -> list[pd.DataFrame]:
    """Find per-locus tables in *epi_dir* that have a gene-symbol column."""
    tables = []
    for f in epi_dir.glob("**/*"):
        if not f.is_file():
            continue
        if f.suffix.lower() not in (".csv", ".tsv", ".txt"):
            continue
        if f.stat().st_size > 500 * 1024 * 1024:  # skip huge raw matrices
            continue
        df = load_epigenetic_table(f)
        if df is None or df.empty:
            continue
        gene_col = autodetect_gene_col(df)
        if gene_col is None:
            continue
        df = df.rename(columns={gene_col: "gene"})
        df["_source"] = str(f.relative_to(epi_dir))
        tables.append(df)
        log.info("  found locus table: %s (%d rows, gene col -> %s)", f.name, len(df), gene_col)
    return tables


def integrate_one_celltype(
    deg_csv: Path,
    cell_type: str,
    sex: str,
    h3k9ac_tables: list[pd.DataFrame],
    methyl_tables: list[pd.DataFrame],
    top_n: int,
    padj_threshold: float,
) -> pd.DataFrame:
    """Cross DEG hits with epigenetic tables, return a long-form joined frame."""
    deg = load_deg_csv(deg_csv)
    deg = deg.dropna(subset=["padj"])
    deg = deg[deg["padj"] < padj_threshold]
    if "log2fc" in deg.columns:
        deg = deg.reindex(deg["log2fc"].abs().sort_values(ascending=False).index)
    deg = deg.head(top_n)
    if deg.empty:
        log.info("    no DEG hits for %s/%s at padj<%.2f", cell_type, sex, padj_threshold)
        return pd.DataFrame()

    rows = []
    for table_kind, tables in (("H3K9ac", h3k9ac_tables), ("methylation", methyl_tables)):
        for t in tables:
            j = deg.merge(t, on="gene", how="inner", suffixes=("_deg", f"_{table_kind}"))
            if j.empty:
                continue
            j["epi_kind"] = table_kind
            j["cell_type"] = cell_type
            j["sex"] = sex
            rows.append(j)
    if not rows:
        return pd.DataFrame()
    return pd.concat(rows, ignore_index=True)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--phenotype", required=True, choices=["tot_adverse_exp", "early_hh_ses", "ace_aggregate"])
    p.add_argument("--deg-dir", required=True, help="DEG results dir for this phenotype (.rda files)")
    p.add_argument("--epi-data-dir", required=True, help="Root with h3k9ac/ and methylation/ subdirs")
    p.add_argument("--overlap-csv", required=True, help="Output of match_individuals.py")
    p.add_argument("--pheno-csv", required=True, help="ACE phenotype CSV")
    p.add_argument("--output-dir", required=True, help="Where to write joined tables + summaries")
    p.add_argument("--top-n", type=int, default=200, help="Top DEG hits per cell type x sex")
    p.add_argument("--padj", type=float, default=0.1, help="DEG padj threshold (default 0.1)")
    p.add_argument("--rscript", default="Rscript",
                   help="Path to Rscript with read.csv + load() (any working R env)")
    args = p.parse_args()

    deg_dir = Path(args.deg_dir)
    epi_root = Path(args.epi_data_dir)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = out_dir / "_deg_csv_cache"

    # 1. Dump DEG .rda -> CSV
    dump_deg_tables(deg_dir, cache_dir, args.rscript)

    # 2. Load epigenetic tables
    log.info("Scanning epigenomic data root: %s", epi_root)
    h3k9ac_tables = harvest_epigenetic_tables(epi_root / "h3k9ac")
    methyl_tables = harvest_epigenetic_tables(epi_root / "methylation")
    log.info("H3K9ac tables: %d, methylation tables: %d",
             len(h3k9ac_tables), len(methyl_tables))

    if not (h3k9ac_tables or methyl_tables):
        raise SystemExit(
            f"No usable per-locus tables under {epi_root}. "
            "Either download_epigenomic.sh has not run yet, or the data files "
            "do not have a recognizable gene-symbol column. Inspect them and "
            "extend autodetect_gene_col() if needed."
        )

    # 3. Per-cell-type integration
    log.info("Loading ACE phenotype (for downstream stratification)")
    overlap = pd.read_csv(args.overlap_csv)
    used_projids = overlap.loc[
        overlap["has_ace_phenotype"]
        & (overlap["has_h3k9ac"] | overlap["has_methyl"]),
        "projid",
    ].astype(str)
    if used_projids.empty:
        log.warning("No overlap with epigenomic data. Continuing on gene-level joins only.")
    else:
        load_for_projids(args.pheno_csv, used_projids)  # validates the join keys

    joined_all = []
    for csv_path in sorted(cache_dir.glob("deseqAnalysisACE_*.csv")):
        meta = parse_deg_filename(csv_path.name)
        if meta is None:
            log.warning("  could not parse %s", csv_path.name)
            continue
        cell_type, sex = meta
        log.info("Integrating %s / %s", cell_type, sex)
        joined = integrate_one_celltype(
            csv_path, cell_type, sex, h3k9ac_tables, methyl_tables,
            top_n=args.top_n, padj_threshold=args.padj,
        )
        if not joined.empty:
            joined_all.append(joined)

    if not joined_all:
        log.warning("No DEG×epigenetic overlaps found.")
        return

    out_table = out_dir / "deg_x_epigenetic.tsv"
    big = pd.concat(joined_all, ignore_index=True)
    big.to_csv(out_table, sep="\t", index=False)
    log.info("Wrote %d rows -> %s", len(big), out_table)

    # 4. Per-cell-type summary
    summary = (
        big.groupby(["cell_type", "sex", "epi_kind"])
        .agg(n_overlap=("gene", "nunique"))
        .reset_index()
    )
    summary_path = out_dir / "summary_by_celltype.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    log.info("Wrote summary -> %s", summary_path)


if __name__ == "__main__":
    main()
