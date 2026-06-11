#!/usr/bin/env python3
"""Identify projid overlap between Tsai snRNA-seq, ACE phenotypes, and
ROSMAP epigenomic datasets (H3K9ac, methylation).

Walks the downloaded epigenomic data trees, harvests every projid it can
find from sample-sheet-like files, intersects with the snRNA-seq cohort,
and writes a summary table.

Usage:
  python match_individuals.py \
      --tsai-h5ad /path/to/tsai_annotated.h5ad \
      --pheno-csv /path/to/ACE_scores.csv \
      --epi-data-dir /path/to/epigenomic_data \
      --output-csv /path/to/individual_overlap.csv
"""

from __future__ import annotations

import argparse
import logging
import os
import re
import sys
from pathlib import Path

import anndata as ad
import pandas as pd

# Shared ACE phenotype loader
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "_shared"))
from load_ace_phenotype import load_for_projids  # noqa: E402

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("match_individuals")

PROJID_RE = re.compile(r"\b\d{7,9}\b")  # ROSMAP projids are 7-9 digit ints


def harvest_projids_from_tree(root: Path, file_glob: str = "**/*") -> set[str]:
    """Walk *root* and pull projid tokens out of any text-like file <50MB.

    ROSMAP sample sheets vary in format (TSV, CSV, manifest TXTs); this
    grabs anything that looks like a projid rather than committing to a
    specific schema.
    """
    if not root.exists():
        return set()

    projids: set[str] = set()
    for f in root.glob(file_glob):
        if not f.is_file():
            continue
        if f.suffix.lower() not in {".csv", ".tsv", ".txt", ".manifest", ""}:
            continue
        try:
            if f.stat().st_size > 50 * 1024 * 1024:
                continue
            with f.open("r", errors="ignore") as fh:
                text = fh.read()
            for m in PROJID_RE.findall(text):
                projids.add(m)
        except (OSError, UnicodeDecodeError):
            continue
    return projids


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--tsai-h5ad", required=True, help="Annotated Tsai snRNA-seq h5ad")
    p.add_argument("--pheno-csv", required=True, help="ACE phenotype CSV")
    p.add_argument("--epi-data-dir", required=True,
                   help="Root of downloaded epigenomic data (contains h3k9ac/ and methylation/)")
    p.add_argument("--output-csv", required=True, help="Output overlap table")
    args = p.parse_args()

    # ----- snRNA-seq projids -----
    log.info("Loading Tsai h5ad: %s", args.tsai_h5ad)
    a = ad.read_h5ad(args.tsai_h5ad, backed="r")
    pid_col = next((c for c in ("projid", "patient_id", "sample_id") if c in a.obs.columns), None)
    if pid_col is None:
        raise SystemExit(f"No projid/patient_id/sample_id column in obs: {list(a.obs.columns)}")
    snrna_projids = set(a.obs[pid_col].astype(str).unique())
    log.info("snRNA-seq cohort: %d projids", len(snrna_projids))

    # ----- ACE phenotypes (dedup'd) -----
    log.info("Loading ACE phenotype CSV: %s", args.pheno_csv)
    pheno = load_for_projids(args.pheno_csv, snrna_projids)
    pheno_projids = set(pheno.dropna(how="all").index)
    log.info("ACE phenotype overlap: %d projids", len(pheno_projids))

    # ----- Epigenomic projids -----
    epi_root = Path(args.epi_data_dir)
    h3k9ac_projids = harvest_projids_from_tree(epi_root / "h3k9ac")
    methyl_projids = harvest_projids_from_tree(epi_root / "methylation")
    log.info("H3K9ac projids harvested:       %d", len(h3k9ac_projids))
    log.info("Methylation projids harvested:  %d", len(methyl_projids))

    if not h3k9ac_projids and not methyl_projids:
        log.warning("No projids harvested from %s. "
                    "Either download_epigenomic.sh has not run, or the sample "
                    "sheets use a non-numeric ID and need a custom parser.", epi_root)

    # ----- Build the overlap table -----
    universe = sorted(snrna_projids | pheno_projids | h3k9ac_projids | methyl_projids)
    rows = []
    for pid in universe:
        rows.append({
            "projid": pid,
            "has_snrna": pid in snrna_projids,
            "has_ace_phenotype": pid in pheno_projids,
            "has_h3k9ac": pid in h3k9ac_projids,
            "has_methyl": pid in methyl_projids,
        })
    df = pd.DataFrame(rows)
    df["n_modalities"] = df[["has_snrna", "has_ace_phenotype", "has_h3k9ac", "has_methyl"]].sum(axis=1)

    out = Path(args.output_csv)
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)

    log.info("Wrote %d rows -> %s", len(df), out)
    log.info("Triple overlap (snRNA + ACE + H3K9ac):       %d",
             ((df.has_snrna) & (df.has_ace_phenotype) & (df.has_h3k9ac)).sum())
    log.info("Triple overlap (snRNA + ACE + methyl):       %d",
             ((df.has_snrna) & (df.has_ace_phenotype) & (df.has_methyl)).sum())
    log.info("Quad overlap (snRNA + ACE + H3K9ac + methyl): %d",
             ((df.has_snrna) & (df.has_ace_phenotype) & (df.has_h3k9ac) & (df.has_methyl)).sum())


if __name__ == "__main__":
    main()
