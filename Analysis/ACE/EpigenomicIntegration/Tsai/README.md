# ACE Epigenomic Integration - Tsai

Canonical entry point:

```bash
bash aceEpiIntT.sh
```

## Status: Scaffold

This directory contains the scaffolded structure for epigenomic integration.
Implementation requires confirming data availability on Synapse and individual
overlap with the ACE cohort.

## Method (Planned)

1. Download ROSMAP H3K9ac ChIP-seq and/or DNA methylation data from Synapse
2. Match individuals by projid between transcriptomic and epigenomic datasets
3. For each ACE-DEG locus:
   - Extract H3K9ac signal at TSS ± 2kb
   - Extract DNA methylation at promoter CpGs
   - Test: regression of epigenetic mark vs ACE with covariates
4. Enrichment: are ACE-DEG loci enriched for differential epigenetic marks
   compared to non-DEG genes? (permutation test)

## Prerequisites

- Synapse credentials configured
- Sufficient individual overlap (target: n > 30 with both ACE + epigenomic data)
- DEG results: `${ACE_OUTPUT_ROOT}/DEG/Tsai/results_${INTEGRATION}/`
- Python environment with pyBigWig, pybedtools, pandas, scipy

## Scripts (To Be Implemented)

- `aceEpiIntT.sh` — Main launcher
- `download_epigenomic.sh` — Synapse download
- `match_individuals.py` — Individual matching
- `epi_integration.py` — Core analysis
- `epi_visualize.py` — Visualization
