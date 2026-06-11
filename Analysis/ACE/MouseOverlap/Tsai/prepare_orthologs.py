#!/usr/bin/env python3
"""Prepare mouse-human ortholog pairs from the MGI homology report."""

from __future__ import annotations

import argparse
import hashlib
import json
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, help="Local path or URL to HOM_MouseHumanSequence.rpt.")
    parser.add_argument("--raw-output", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--metadata-output", type=Path, required=True)
    return parser.parse_args()


def is_url(source: str) -> bool:
    return source.startswith("http://") or source.startswith("https://")


def fetch_source(source: str, raw_output: Path) -> bytes:
    raw_output.parent.mkdir(parents=True, exist_ok=True)
    if is_url(source):
        with urllib.request.urlopen(source, timeout=120) as response:
            content = response.read()
        raw_output.write_bytes(content)
        return content

    path = Path(source).expanduser()
    if not path.exists():
        raise FileNotFoundError(f"Ortholog source does not exist: {path}")
    content = path.read_bytes()
    if path.resolve() != raw_output.resolve():
        raw_output.write_bytes(content)
    return content


def prepare_pairs(raw_path: Path) -> pd.DataFrame:
    hom = pd.read_csv(raw_path, sep="\t", dtype=str).fillna("")
    required = {
        "DB Class Key",
        "Common Organism Name",
        "NCBI Taxon ID",
        "Symbol",
        "EntrezGene ID",
        "Mouse MGI ID",
        "HGNC ID",
    }
    missing = required - set(hom.columns)
    if missing:
        raise ValueError("MGI report missing columns: " + ", ".join(sorted(missing)))

    mouse = hom[hom["Common Organism Name"] == "mouse, laboratory"].copy()
    human = hom[hom["Common Organism Name"] == "human"].copy()

    mouse_counts = mouse.groupby("DB Class Key")["Symbol"].nunique().rename("n_mouse_in_class")
    human_counts = human.groupby("DB Class Key")["Symbol"].nunique().rename("n_human_in_class")

    merged = mouse.merge(
        human,
        on="DB Class Key",
        suffixes=("_mouse", "_human"),
        how="inner",
    )
    merged = merged.merge(mouse_counts, on="DB Class Key", how="left")
    merged = merged.merge(human_counts, on="DB Class Key", how="left")

    out = pd.DataFrame(
        {
            "ortholog_source": "MGI_HOM_MouseHumanSequence",
            "db_class_key": merged["DB Class Key"],
            "mouse_symbol": merged["Symbol_mouse"],
            "human_symbol": merged["Symbol_human"],
            "mouse_entrez": merged["EntrezGene ID_mouse"],
            "human_entrez": merged["EntrezGene ID_human"],
            "mouse_mgi_id": merged["Mouse MGI ID_mouse"],
            "hgnc_id": merged["HGNC ID_human"],
            "n_mouse_in_class": merged["n_mouse_in_class"].astype(int),
            "n_human_in_class": merged["n_human_in_class"].astype(int),
        }
    )
    out["is_one_to_one"] = (out["n_mouse_in_class"] == 1) & (out["n_human_in_class"] == 1)
    out = out.drop_duplicates().sort_values(["db_class_key", "mouse_symbol", "human_symbol"])
    return out


def main() -> None:
    args = parse_args()
    content = fetch_source(args.source, args.raw_output)
    sha256 = hashlib.sha256(content).hexdigest()
    pairs = prepare_pairs(args.raw_output)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    pairs.to_csv(args.output, index=False)

    metadata = {
        "source": args.source,
        "raw_output": str(args.raw_output),
        "prepared_output": str(args.output),
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        "sha256": sha256,
        "n_pairs": int(len(pairs)),
        "n_one_to_one_pairs": int(pairs["is_one_to_one"].sum()),
        "n_homology_classes": int(pairs["db_class_key"].nunique()),
    }
    args.metadata_output.write_text(json.dumps(metadata, indent=2) + "\n")

    print(f"Wrote: {args.output}")
    print(f"Wrote: {args.metadata_output}")
    print(f"Pairs: {metadata['n_pairs']} one_to_one={metadata['n_one_to_one_pairs']}")


if __name__ == "__main__":
    main()

