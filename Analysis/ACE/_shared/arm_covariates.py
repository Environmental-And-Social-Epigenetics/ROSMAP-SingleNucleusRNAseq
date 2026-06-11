"""Single source of truth for the non-ANCOVA male ACE AD-model arms.

Defines, per arm, the extra AD covariates (added to the baseline `age_death + pmi`),
the interaction term (if any), the DEG `.rda` object name and filename suffix, and the
niareagansc cohort filter. Imported by the downstream TF / SCENIC analyses so the per-arm
association model matches each DEG arm exactly. Keep in sync with arm_covariates.sh / .R.

All arms are males-only (msex == 1), phenotype = tot_adverse_exp, integration = derived_batch.
`AD_binary` = 1 if niareagansc in {1,2} (AD pathology present) else 0; defined on the fly from
niareagansc when the arm needs it.
"""

from __future__ import annotations

# arm -> spec
ARM_SPECS = {
    "MaleNoADadj": {
        "ad_covars": [],                       # ~ age_death + pmi + phenotype
        "interaction": False,
        "needs_ad_binary": False,
        "deg_obj": "res",
        "deg_suffix": "MaleNoADadj",
        "nia_filter": [1, 2, 3, 4],
    },
    "MaleNiaReagan": {
        "ad_covars": ["niareagansc"],
        "interaction": False,
        "needs_ad_binary": False,
        "deg_obj": "res",
        "deg_suffix": "MaleNiaReagan",
        "nia_filter": [1, 2, 3, 4],
    },
    "MaleBinaryAD": {
        "ad_covars": ["AD_binary"],
        "interaction": False,
        "needs_ad_binary": True,
        "deg_obj": "res_ace",
        "deg_suffix": "MaleBinaryAD",
        "nia_filter": [1, 2, 3, 4],
    },
    "MaleContAD": {
        "ad_covars": ["amylsqrt", "tangsqrt"],
        "interaction": False,
        "needs_ad_binary": False,
        "deg_obj": "res",
        "deg_suffix": "MaleContAD",
        "nia_filter": [1, 2, 3, 4],
    },
    "MaleAceByAD": {
        "ad_covars": ["AD_binary"],
        "interaction": True,                    # AD_binary:phenotype
        "needs_ad_binary": True,
        "deg_obj": "res_ace",
        "deg_suffix": "MaleAceByAD",
        "nia_filter": [1, 2, 3, 4],
    },
}

NON_ANCOVA_ARMS = list(ARM_SPECS.keys())

# Cell types (signal-bearing) used across all downstream analyses. Note: the pooled
# inhibitory label is "Inh" in DEG outputs but the h5ad / SCENIC dir is "broad_Inh".
SIGNAL_CELLTYPES = ["In-PV_Basket", "Ast", "Mic", "Oli", "OPC", "Inh"]

# pheno-CSV columns needed to build any arm's covariates
PHENO_COLS = ["projid", "msex", "age_death", "pmi", "niareagansc",
              "amylsqrt", "tangsqrt", "tot_adverse_exp"]


def get_spec(arm: str) -> dict:
    if arm not in ARM_SPECS:
        raise KeyError(f"Unknown arm '{arm}'. Known: {NON_ANCOVA_ARMS}")
    return ARM_SPECS[arm]


def add_ad_binary(df, nia_col: str = "niareagansc", out_col: str = "AD_binary"):
    """Add the AD_binary column (1 if niareagansc in {1,2}) to a DataFrame in place."""
    df[out_col] = df[nia_col].isin([1, 2]).astype(int)
    return df


def ols_formula(arm: str, phenotype: str, response: str = "y") -> str:
    """Build a patsy/statsmodels OLS formula for an arm:
        y ~ phenotype + age_death + pmi [+ ad_covars...] [+ AD_binary:phenotype]
    """
    spec = get_spec(arm)
    terms = [phenotype, "age_death", "pmi", *spec["ad_covars"]]
    if spec["interaction"]:
        terms.append(f"AD_binary:{phenotype}")
    return f"{response} ~ " + " + ".join(terms)
