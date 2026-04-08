#!/usr/bin/env python3
"""Generate deterministic ACE smoke-test fixtures for Tsai and DeJager."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


ACE_COMPONENTS = [
    "emotional_neglect",
    "family_pro_sep",
    "financial_need",
    "parental_intimidation",
    "parental_violence",
]
GENES = ["GeneA", "GeneB", "GeneC"]


def build_adata(obs_names: list[str], obs_df: pd.DataFrame) -> ad.AnnData:
    matrix = np.arange(1, len(obs_names) * len(GENES) + 1, dtype=np.int32).reshape(len(obs_names), len(GENES))
    return ad.AnnData(
        X=sp.csr_matrix(matrix),
        obs=obs_df.copy(),
        var=pd.DataFrame(index=GENES),
    )


def categoricalize(obs_df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    result = obs_df.copy()
    for column in columns:
        result[column] = pd.Categorical(result[column].astype(str))
    return result


def phenotype_rows() -> list[dict[str, float | int | str]]:
    rows = []
    all_ids = ["1001", "1002", "1003", "1004", "2001", "2002", "2003", "2004"]
    for idx, projid in enumerate(all_ids, start=1):
        rows.append(
            {
                "projid": projid,
                "tot_adverse_exp": float(idx),
                "early_hh_ses": float(10 - idx),
                "msex": idx % 2,
                "age_death": float(70 + idx),
                "pmi": float(4 + idx / 10),
                "niareagansc": float(idx % 4),
                "emotional_neglect": float(idx % 3),
                "family_pro_sep": float((idx + 1) % 3),
                "financial_need": float((idx + 2) % 3),
                "parental_intimidation": float(idx % 2),
                "parental_violence": float((idx + 1) % 2),
            }
        )
    return rows


def write_tsai(root: Path) -> dict[str, str]:
    annotated_dir = root / "tsai" / "annotated"
    doublet_dir = root / "tsai" / "doublet_removed"
    annotated_dir.mkdir(parents=True, exist_ok=True)
    doublet_dir.mkdir(parents=True, exist_ok=True)

    samples = {
        "1001": ["Ex-L2_IT", "Ex-L2_IT", "Ast"],
        "1002": ["Ex-L2_IT", "In-L2", "Ast"],
        "1003": ["Ex-L3_IT", "Ex-L3_IT", "Ast"],
        "1004": ["Ex-L3_IT", "In-L2", "Ast"],
    }

    annotated_obs = []
    annotated_names = []
    for sample_id, cell_types in samples.items():
        singlet_names = []
        singlet_obs = []
        for idx, cell_type in enumerate(cell_types, start=1):
            base = f"tsai_{sample_id}_cell{idx}"
            annotated_names.append(f"{base}-1")
            singlet_names.append(base)
            singlet_obs.append({"sample_id": sample_id, "projid": sample_id})
        build_adata(singlet_names, pd.DataFrame(singlet_obs, index=singlet_names)).write_h5ad(
            doublet_dir / f"{sample_id}_singlets.h5ad"
        )

    annotated_obs = []
    for sample_id, cell_types in samples.items():
        for idx, cell_type in enumerate(cell_types, start=1):
            annotated_obs.append({"cell_type": cell_type, "sample_id": sample_id, "projid": sample_id})
    annotated_df = categoricalize(
        pd.DataFrame(annotated_obs, index=annotated_names),
        ["cell_type", "sample_id", "projid"],
    )
    annotated = build_adata(annotated_names, annotated_df)
    annotated_path = annotated_dir / "tsai_annotated.h5ad"
    annotated.write_h5ad(annotated_path)
    return {"annotated_h5ad": str(annotated_path), "doublet_dir": str(doublet_dir)}


def write_dejager(root: Path) -> dict[str, str]:
    annotated_dir = root / "dejager" / "annotated"
    doublet_dir = root / "dejager" / "doublet_removed"
    annotated_dir.mkdir(parents=True, exist_ok=True)
    doublet_dir.mkdir(parents=True, exist_ok=True)

    libraries = {
        "LIB001": ("2001", ["Ex-L2_IT", "Ex-L2_IT", "Ast"]),
        "LIB002": ("2002", ["Ex-L2_IT", "In-L2", "Ast"]),
        "LIB003": ("2003", ["Ex-L3_IT", "Ex-L3_IT", "Ast"]),
        "LIB004": ("2004", ["Ex-L3_IT", "In-L2", "Ast"]),
    }

    annotated_obs = []
    annotated_names = []
    for library_id, (patient_id, cell_types) in libraries.items():
        singlet_names = []
        singlet_obs = []
        for idx, cell_type in enumerate(cell_types, start=1):
            base = f"dj_{library_id}_cell{idx}"
            annotated_names.append(f"{base}-1")
            singlet_names.append(base)
            annotated_obs.append(
                {
                    "cell_type": cell_type,
                    "library_id": library_id,
                    "patient_id": patient_id,
                    "projid": patient_id,
                }
            )
            singlet_obs.append({"library_id": library_id, "patient_id": patient_id, "projid": patient_id})
        build_adata(singlet_names, pd.DataFrame(singlet_obs, index=singlet_names)).write_h5ad(
            doublet_dir / f"{library_id}_singlets.h5ad"
        )

    annotated_df = categoricalize(
        pd.DataFrame(annotated_obs, index=annotated_names),
        ["cell_type", "library_id", "patient_id", "projid"],
    )
    annotated = build_adata(annotated_names, annotated_df)
    annotated_path = annotated_dir / "dejager_annotated.h5ad"
    annotated.write_h5ad(annotated_path)
    return {"annotated_h5ad": str(annotated_path), "doublet_dir": str(doublet_dir)}


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate ACE smoke-test fixtures.")
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--cohort", choices=["tsai", "dejager", "both"], default="both")
    args = parser.parse_args()

    root = args.output_root.resolve()
    root.mkdir(parents=True, exist_ok=True)

    phenotype_path = root / "ace_scores.csv"
    pd.DataFrame(phenotype_rows()).to_csv(phenotype_path, index=False)

    manifest: dict[str, object] = {"ace_scores_csv": str(phenotype_path)}
    if args.cohort in {"tsai", "both"}:
        manifest["tsai"] = write_tsai(root)
    if args.cohort in {"dejager", "both"}:
        manifest["dejager"] = write_dejager(root)

    manifest_path = root / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Wrote ACE smoke fixtures to {root}")
    print(f"Manifest: {manifest_path}")


if __name__ == "__main__":
    main()
