"""Shared ACE phenotype loader for Tsai / DeJager workflows.

The canonical CSV (TSAI_DEJAGER_all_patients_wACEscores.csv) pools both
cohorts and contains duplicate projids. Any workflow that does
pheno.set_index("projid").reindex(cohort_ids) must filter to a single cohort
first or pandas raises ValueError("cannot reindex on an axis with duplicate
labels"). Use load_for_projids() as the single entry point.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Union

import numpy as np
import pandas as pd

PathLike = Union[str, Path]

REQUIRED_COLUMNS = [
    "projid",
    "tot_adverse_exp",
    "early_hh_ses",
    "msex",
    "age_death",
    "pmi",
    "niareagansc",
]

ACE_COMPONENTS = [
    "emotional_neglect",
    "family_pro_sep",
    "financial_need",
    "parental_intimidation",
    "parental_violence",
]


def _add_ace_aggregate(pheno: pd.DataFrame) -> pd.DataFrame:
    if "ace_aggregate" in pheno.columns:
        return pheno
    if not {"tot_adverse_exp", "early_hh_ses"}.issubset(pheno.columns):
        return pheno
    z_adv = (pheno["tot_adverse_exp"] - pheno["tot_adverse_exp"].mean()) / pheno["tot_adverse_exp"].std()
    z_ses = (pheno["early_hh_ses"] - pheno["early_hh_ses"].mean()) / pheno["early_hh_ses"].std()
    pheno = pheno.copy()
    pheno["ace_aggregate"] = z_adv - z_ses
    return pheno


def load_raw(pheno_csv: PathLike) -> pd.DataFrame:
    """Load the ACE phenotype CSV with schema validation and ace_aggregate derived."""
    pheno = pd.read_csv(pheno_csv)
    # Drop the R-written rowname column if present ("" header → "Unnamed: 0").
    # Carrying it makes drop_duplicates() ineffective because every row has a
    # unique rownumber, leaving spurious projid duplicates in the pooled CSV.
    if "Unnamed: 0" in pheno.columns:
        pheno = pheno.drop(columns=["Unnamed: 0"])
    missing = set(REQUIRED_COLUMNS) - set(pheno.columns)
    if missing:
        raise ValueError(f"Phenotype CSV missing columns: {sorted(missing)}")
    pheno["projid"] = pheno["projid"].astype(str)
    return _add_ace_aggregate(pheno)


def load_for_projids(
    pheno_csv: PathLike,
    projids: Iterable[str],
    require_all: bool = False,
) -> pd.DataFrame:
    """Load phenotypes restricted to *projids* (e.g. pseudobulk obs_names).

    Filters to rows whose projid is in *projids*, drops exact duplicate rows,
    and raises if any projid still appears more than once (non-exact
    duplicate — indicates real data conflict, not pool-artifact duplication).

    Parameters
    ----------
    pheno_csv
        Path to the ACE phenotype CSV.
    projids
        Iterable of projid strings to keep.
    require_all
        If True, raise when any projid in *projids* is missing from the CSV.

    Returns
    -------
    DataFrame indexed by projid (str), ordered to match *projids*.
    """
    wanted = pd.Index([str(p) for p in projids], name="projid")
    pheno = load_raw(pheno_csv)

    sub = pheno[pheno["projid"].isin(wanted)].copy()
    sub = sub.drop_duplicates()

    dup_mask = sub["projid"].duplicated(keep=False)
    if dup_mask.any():
        offenders = sorted(sub.loc[dup_mask, "projid"].unique().tolist())
        raise ValueError(
            "Phenotype CSV has conflicting rows for projid(s): "
            f"{offenders[:10]}{' ...' if len(offenders) > 10 else ''} "
            "(exact duplicates were dropped; these remain after dedup)."
        )

    sub = sub.set_index("projid")
    sub = sub.reindex(wanted)

    if require_all and sub.isna().all(axis=1).any():
        missing = wanted[sub.isna().all(axis=1)].tolist()
        raise ValueError(
            f"{len(missing)} projid(s) from cohort not found in phenotype CSV: "
            f"{missing[:10]}{' ...' if len(missing) > 10 else ''}"
        )
    return sub
