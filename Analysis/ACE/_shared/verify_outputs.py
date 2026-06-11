#!/usr/bin/env python3
"""Verify ACE Tsai analysis outputs are complete.

Walks the expected output directory tree under ${ANALYSIS_OUTPUT_ROOT}/ACE
and prints a one-line PASS/PARTIAL/MISSING verdict per workflow.

Usage:
  python verify_outputs.py [--workflow DEG] [--phenotype tot_adverse_exp]
                           [--integration derived_batch] [--strict]

Without filters, runs every workflow.
--strict: exit 1 if any workflow is not PASS.
"""

from __future__ import annotations

import argparse
import os
import sys
from dataclasses import dataclass
from pathlib import Path

DEFAULT_ROOT = os.environ.get(
    "ACE_OUTPUT_ROOT",
    str(Path(__file__).resolve().parents[3].parent / "Analysis_Outputs" / "ACE"),
)

PHENOTYPES = ("tot_adverse_exp", "early_hh_ses", "ace_aggregate")
SEXES = ("Female", "Male")

# Cell-type lists per workflow are best-effort. Workflows that operate on
# all 19 cell types use CELLTYPES_FULL; SCENIC uses a curated subset.
CELLTYPES_FULL = (
    "Ast", "Endo", "Ex-L2_3", "Ex-L4", "Ex-L4_5", "Ex-L5", "Ex-L5_6",
    "Ex-L5_6-CC", "Ex-NRGN", "In-LAMP5", "In-PV_Basket", "In-PV_Chandelier",
    "In-SST", "In-VIP", "Inh", "Mic", "Oli", "OPC", "Per",
)


@dataclass
class WorkflowResult:
    name: str
    expected: int
    found: int
    examples_missing: list[str]

    @property
    def status(self) -> str:
        if self.expected == 0:
            return "N/A"
        if self.found == self.expected:
            return "PASS"
        if self.found == 0:
            return "MISSING"
        return "PARTIAL"

    def __str__(self) -> str:
        ratio = f"{self.found}/{self.expected}"
        if self.status == "PASS":
            return f"  [PASS]    {self.name:<24} {ratio}"
        marker = "[MISS]   " if self.status == "MISSING" else "[PARTIAL]"
        ex = ""
        if self.examples_missing:
            ex = " e.g. " + ", ".join(self.examples_missing[:3])
        return f"  {marker} {self.name:<24} {ratio}{ex}"


def check(label: str, expected_paths: list[Path]) -> WorkflowResult:
    missing = [str(p.relative_to(p.anchor)) for p in expected_paths if not p.exists() or p.stat().st_size == 0]
    found = len(expected_paths) - len(missing)
    return WorkflowResult(label, len(expected_paths), found, missing)


def deg_paths(root: Path, integration: str, phenotype: str) -> list[Path]:
    base = root / "DEG" / "Tsai" / f"results_{integration}" / phenotype
    paths = []
    for sex in ("Fem", "Male"):
        for ct in CELLTYPES_FULL:
            paths.append(base / f"deseqAnalysisACE_{phenotype}_{ct}_{sex}.rda")
    return paths


def gsea_paths(root: Path, integration: str, phenotype: str) -> list[Path]:
    # GSEA writes per-celltype × per-database RDS. We only check that at least
    # one RDS exists per (celltype, sex) — DBs vary by what survives filtering.
    base = root / "GSEA" / "Tsai" / f"results_{integration}" / phenotype
    paths = []
    for sex in SEXES:
        for ct in CELLTYPES_FULL:
            # require ranked_genes.csv as the most reliable proxy
            paths.append(base / sex / f"{ct}_ranked_genes.csv")
    return paths


def celltype_proportion_paths(root: Path, integration: str, phenotype: str) -> list[Path]:
    base = root / "CellTypeProportion" / "Tsai" / f"results_{integration}"
    if phenotype == "tot_adverse_exp":
        sub = "primary"
    elif phenotype == "early_hh_ses":
        sub = "ses"
    elif phenotype == "ace_aggregate":
        sub = "aggregate"
    else:
        return []
    return [
        base / sub / f"sccomp_{phenotype}_all_fine.rds",
        base / sub / f"sccomp_{phenotype}_female_fine.rds",
        base / sub / f"sccomp_{phenotype}_male_fine.rds",
    ]


def tfactivity_paths(root: Path, integration: str, phenotype: str) -> list[Path]:
    base = root / "TFActivity" / "Tsai" / f"results_{integration}" / phenotype
    return [base / "tf_summary.csv"]


def scenic_paths(root: Path, integration: str, phenotype: str) -> list[Path]:
    base = root / "SCENIC" / "Tsai" / f"results_{integration}" / phenotype
    # SCENIC submits per-celltype × sex; cell-type list lives in the launcher.
    male_cts = ("broad_Inh", "Mic", "Ast", "In-PV_Basket", "Oli", "broad_Exc")
    female_cts = ("Oli", "Ast", "Ex-L2_3", "broad_Exc", "OPC", "broad_Inh", "Mic")
    paths = []
    for ct in male_cts:
        paths.append(base / f"Male_{ct}" / "regression_results.csv")
    for ct in female_cts:
        paths.append(base / f"Female_{ct}" / "regression_results.csv")
    return paths


def cellchat_paths(root: Path, integration: str, phenotype: str) -> list[Path]:
    base = root / "CellChat" / "Tsai" / f"results_{integration}" / phenotype
    # CellChat writes a mergedObj per sex; require a marker file.
    return [base / sex / "differential_summary.csv" for sex in SEXES]


def hdwgcna_paths(root: Path, integration: str, phenotype: str) -> list[Path]:
    base = root / "hdWGCNA" / "Tsai" / f"results_{integration}" / phenotype
    male_cts = ("broad_Inh", "Mic", "Ast", "In-PV_Basket", "Oli", "broad_Exc")
    female_cts = ("Oli", "Ex-L2_3", "Ast", "broad_Exc", "OPC")
    paths = []
    for ct in male_cts:
        paths.append(base / f"Male_{ct}" / "module_assignments.csv")
    for ct in female_cts:
        paths.append(base / f"Female_{ct}" / "module_assignments.csv")
    return paths


def micstate_paths(root: Path, integration: str, phenotype: str) -> list[Path]:
    base = root / "MicState" / "Tsai" / f"results_{integration}" / phenotype
    return [base / sex / "state_regression_results.csv" for sex in SEXES]


def epi_paths(root: Path, integration: str, phenotype: str) -> list[Path]:
    base = root / "EpigenomicIntegration" / "Tsai" / phenotype
    return [base / "deg_x_epigenetic.tsv", base / "summary_by_celltype.tsv"]


WORKFLOWS = {
    "DEG": deg_paths,
    "CellTypeProportion": celltype_proportion_paths,
    "GSEA": gsea_paths,
    "TFActivity": tfactivity_paths,
    "SCENIC": scenic_paths,
    "CellChat": cellchat_paths,
    "hdWGCNA": hdwgcna_paths,
    "MicState": micstate_paths,
    "EpigenomicIntegration": epi_paths,
}


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--root", default=DEFAULT_ROOT, help="Analysis output root containing ACE/")
    p.add_argument("--workflow", choices=list(WORKFLOWS), help="Single workflow to check")
    p.add_argument("--phenotype", choices=PHENOTYPES, help="Single phenotype to check")
    p.add_argument("--integration", default="derived_batch")
    p.add_argument("--strict", action="store_true",
                   help="Exit non-zero if any workflow is not PASS")
    args = p.parse_args()

    root = Path(args.root)
    workflows = [args.workflow] if args.workflow else list(WORKFLOWS)
    phenotypes = [args.phenotype] if args.phenotype else list(PHENOTYPES)

    print(f"=== ACE Tsai output verification ({args.integration}) ===")
    print(f"Root: {root}")
    print()

    any_not_pass = False
    for wf in workflows:
        for phen in phenotypes:
            paths = WORKFLOWS[wf](root, args.integration, phen)
            if not paths:
                continue
            r = check(f"{wf}/{phen}", paths)
            print(r)
            if r.status not in ("PASS", "N/A"):
                any_not_pass = True

    if args.strict and any_not_pass:
        sys.exit(1)


if __name__ == "__main__":
    main()
