#!/usr/bin/env python3
"""
Prepare sample-by-cell-type count tables for ACE sccomp analysis on DeJager data.

The script reads only the AnnData obs table from the integrated annotated h5ad,
validates the ACE phenotype inputs, and writes the generated count tables under
ANALYSIS_OUTPUT_ROOT rather than back into the repo tree.
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.stats import zscore

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
sys.stdout.reconfigure(line_buffering=True)


ACE_COMPONENTS = [
    "emotional_neglect",
    "family_pro_sep",
    "financial_need",
    "parental_intimidation",
    "parental_violence",
]
REQUIRED_PHENOTYPE_COLUMNS = [
    "projid",
    "tot_adverse_exp",
    "early_hh_ses",
    "msex",
    "age_death",
    "pmi",
    "niareagansc",
    *ACE_COMPONENTS,
]


def repo_root() -> Path:
    return Path(__file__).resolve().parents[4]


def default_analysis_root(cohort: str) -> Path:
    configured_root = os.environ.get("ANALYSIS_OUTPUT_ROOT")
    if configured_root:
        return Path(configured_root) / "ACE" / "CellTypeProportion" / cohort
    return repo_root().parent / "Analysis_Outputs" / "ACE" / "CellTypeProportion" / cohort


def default_h5ad_paths() -> dict[str, Path]:
    repo = repo_root()
    dejager_processing = repo / "Data" / "Transcriptomics" / "DeJager" / "Processing_Outputs"
    return {
        "library_id": Path(
            os.environ.get("DEJAGER_INTEGRATED", str(dejager_processing / "03_Integrated"))
        )
        / "dejager_annotated.h5ad",
        "derived_batch": Path(
            os.environ.get(
                "DEJAGER_INTEGRATED_DERIVED_BATCH",
                str(dejager_processing / "03_Integrated_derived_batch"),
            )
        )
        / "dejager_annotated.h5ad",
        "patient_id": Path(
            os.environ.get(
                "DEJAGER_INTEGRATED_PATIENT_ID",
                str(dejager_processing / "03_Integrated_patient_id"),
            )
        )
        / "dejager_annotated.h5ad",
        "pool_batch": Path(
            os.environ.get(
                "DEJAGER_INTEGRATED_POOL_BATCH",
                str(dejager_processing / "03_Integrated_pool_batch"),
            )
        )
        / "dejager_annotated.h5ad",
    }


def default_pheno_csv() -> Path:
    return Path(
        os.environ.get(
            "ACE_SCORES_CSV",
            str(repo_root() / "Data" / "Phenotypes" / "TSAI_DEJAGER_all_patients_wACEscores.csv"),
        )
    )


def normalize_integration(name: str) -> str:
    aliases = {
        "batch": "library_id",
        "projid": "patient_id",
        "library_id": "library_id",
        "patient_id": "patient_id",
        "pool_batch": "pool_batch",
        "derived_batch": "derived_batch",
    }
    if name not in aliases:
        choices = ", ".join(sorted(aliases))
        raise SystemExit(f"ERROR: unsupported integration '{name}'. Expected one of: {choices}")
    return aliases[name]


def ensure_output_root(output_root: Path) -> None:
    repo = repo_root().resolve()
    resolved = output_root.resolve()
    if resolved.is_relative_to(repo):
        raise SystemExit(
            f"ERROR: output root must live outside the repo tree.\n"
            f"  repo:   {repo}\n"
            f"  output: {resolved}"
        )


def broad_group(cell_type: str) -> str:
    if cell_type.startswith(("Ex-", "Ex_")):
        return "Exc"
    if cell_type.startswith(("In-", "In_")):
        return "Inh"
    return cell_type


def decode_values(values) -> list[str]:
    return [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values]


def read_obs_column(obs_group: h5py.Group, column: str) -> pd.Series:
    if column not in obs_group:
        raise KeyError(column)
    node = obs_group[column]
    if isinstance(node, h5py.Group) and "codes" in node and "categories" in node:
        codes = node["codes"][:]
        categories = decode_values(node["categories"][:])
        return pd.Series(pd.Categorical.from_codes(codes, categories=categories))
    return pd.Series(decode_values(node[:]), dtype="string")


def read_obs_h5py(h5ad_path: Path) -> pd.DataFrame:
    with h5py.File(h5ad_path, "r") as handle:
        obs = handle["obs"]
        missing = [field for field in ["cell_type"] if field not in obs]
        if missing:
            raise SystemExit(f"ERROR: required obs fields missing from {h5ad_path}: {', '.join(missing)}")

        sample_field = None
        for candidate in ["patient_id", "sample_id", "projid"]:
            if candidate in obs:
                sample_field = candidate
                break
        if sample_field is None:
            raise SystemExit("ERROR: h5ad obs must contain one of: patient_id, sample_id, projid")

        return pd.DataFrame(
            {
                "cell_type": read_obs_column(obs, "cell_type"),
                "sample_id": read_obs_column(obs, sample_field).astype("string"),
            }
        )


def load_phenotypes(pheno_csv: Path) -> pd.DataFrame:
    if not pheno_csv.exists():
        raise SystemExit(f"ERROR: phenotype CSV not found: {pheno_csv}")
    phenotypes = pd.read_csv(pheno_csv)
    missing = [column for column in REQUIRED_PHENOTYPE_COLUMNS if column not in phenotypes.columns]
    if missing:
        raise SystemExit(
            f"ERROR: phenotype CSV is missing required columns: {', '.join(missing)}"
        )
    phenotypes["projid"] = phenotypes["projid"].astype(str)
    phenotypes["ace_aggregate"] = np.nan
    ace_mask = phenotypes["tot_adverse_exp"].notna() & phenotypes["early_hh_ses"].notna()
    if ace_mask.any():
        phenotypes.loc[ace_mask, "ace_aggregate"] = (
            zscore(phenotypes.loc[ace_mask, "tot_adverse_exp"])
            + zscore(phenotypes.loc[ace_mask, "early_hh_ses"])
        )
    return phenotypes


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare ACE sccomp count tables for DeJager")
    parser.add_argument(
        "--integration",
        default="library_id",
        choices=["batch", "projid", "library_id", "patient_id", "pool_batch", "derived_batch"],
        help="Integrated annotated object to summarize.",
    )
    parser.add_argument(
        "--h5ad",
        type=Path,
        default=None,
        help="Override the annotated h5ad path for the chosen integration.",
    )
    parser.add_argument(
        "--pheno-csv",
        type=Path,
        default=default_pheno_csv(),
        help="ACE phenotype CSV. Defaults to ACE_SCORES_CSV or the tracked phenotype table.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=default_analysis_root("DeJager"),
        help="Output root under ANALYSIS_OUTPUT_ROOT/ACE/CellTypeProportion/DeJager by default.",
    )
    args = parser.parse_args()

    integration = normalize_integration(args.integration)
    h5ad_path = args.h5ad or default_h5ad_paths()[integration]
    output_root = args.output_root
    ensure_output_root(output_root)

    if not h5ad_path.exists():
        raise SystemExit(f"ERROR: annotated h5ad not found: {h5ad_path}")

    data_dir = output_root / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reading obs metadata from {h5ad_path}")
    obs_df = read_obs_h5py(h5ad_path)
    print(f"  Total cells: {len(obs_df):,}")
    print(f"  Unique samples: {obs_df['sample_id'].nunique()}")

    fine_wide = pd.crosstab(obs_df["sample_id"], obs_df["cell_type"])
    fine_wide.index.name = "sample"
    fine_wide["total_cells"] = fine_wide.sum(axis=1)

    obs_df["broad_type"] = obs_df["cell_type"].astype(str).map(broad_group)
    broad_wide = pd.crosstab(obs_df["sample_id"], obs_df["broad_type"])
    broad_wide.index.name = "sample"
    broad_wide["total_cells"] = broad_wide.sum(axis=1)

    fine_long = fine_wide.drop(columns=["total_cells"]).reset_index().melt(
        id_vars="sample",
        var_name="cell_group",
        value_name="count",
    )
    broad_long = broad_wide.drop(columns=["total_cells"]).reset_index().melt(
        id_vars="sample",
        var_name="cell_group",
        value_name="count",
    )

    phenotypes = load_phenotypes(args.pheno_csv)
    sample_totals = fine_wide[["total_cells"]].reset_index().rename(columns={"sample": "projid"})
    sample_totals["projid"] = sample_totals["projid"].astype(str)
    metadata = sample_totals.merge(phenotypes, on="projid", how="left")

    print(f"Loading phenotypes from {args.pheno_csv}")
    print(f"  Samples in h5ad: {len(metadata)}")
    print(f"  Samples with ACE data: {metadata['tot_adverse_exp'].notna().sum()}")
    print(f"  Samples with ace_aggregate: {metadata['ace_aggregate'].notna().sum()}")

    fine_path = data_dir / f"cell_counts_fine_{integration}.csv"
    broad_path = data_dir / f"cell_counts_broad_{integration}.csv"
    meta_path = data_dir / f"metadata_{integration}.csv"

    fine_long.to_csv(fine_path, index=False)
    broad_long.to_csv(broad_path, index=False)
    metadata.to_csv(meta_path, index=False)

    print(f"Outputs written under {data_dir}")
    print(f"  {fine_path.name}")
    print(f"  {broad_path.name}")
    print(f"  {meta_path.name}")


if __name__ == "__main__":
    main()
