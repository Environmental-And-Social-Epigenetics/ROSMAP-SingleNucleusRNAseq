# {Phenotype Name}

> Copy this template directory to create a new phenotype analysis.
> Rename the directory and fill in the sections below.
>
> **New here?** Read [`HOWTO_START.md`](HOWTO_START.md) first — it walks through
> copying the template, the `<AnalysisType>/<Cohort>` layout, the launcher
> pattern, where outputs go, and how to add a smoke test. A ready-to-edit
> launcher skeleton lives at `DEG/Tsai/run_analysis.sh.template`.

## Phenotype Definition

{Describe what this phenotype measures and its relevance to AD research.}

## Patient Selection

{How are patients classified? Binary (exposed/unexposed)? Continuous score? Multiple groups?}

| Group | Criteria | N (Tsai) | N (DeJager) |
|-------|----------|----------|-------------|
| {Group 1} | {criteria} | - | - |
| {Group 2} | {criteria} | - | - |

## Key Comparisons

- {Group 1} vs {Group 2}, stratified by cell type
- {Additional comparisons}

## Phenotype Data

Source: `Data/Phenotypes/{relevant_file}.csv`
Key columns: `projid`, {phenotype columns}
