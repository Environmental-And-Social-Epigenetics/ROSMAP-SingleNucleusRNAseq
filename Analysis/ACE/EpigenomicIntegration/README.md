# ACE Epigenomic Integration

Multi-modal integration of ACE transcriptomic results with published ROSMAP
epigenomic data (H3K9ac ChIP-seq, DNA methylation) to directly test the
epigenetic repression hypothesis.

## Biological Motivation

The DEG analysis shows 1,613 downregulated genes in males with high ACE
exposure, suggesting widespread transcriptional repression. The epigenomic
integration tests whether these same loci show reduced H3K9ac (an activating
histone mark) or increased DNA methylation — providing direct evidence that
the transcriptional repression operates through epigenetic mechanisms.

## Data Sources

| Data Type | Synapse ID | Description |
|-----------|------------|-------------|
| H3K9ac ChIP-seq | syn7428143 (approx.) | Histone acetylation marks, bulk DLPFC |
| DNA methylation | syn3157275 (approx.) | 450K array or RRBS, DLPFC |

**Note:** Exact Synapse IDs should be verified against the current ROSMAP data
dictionary before downloading. Individual overlap with the ACE phenotype
cohort must be confirmed.

## Status

This analysis is scaffolded with directory structure and documentation.
Full implementation requires:

1. Confirming Synapse data availability and individual overlap (n > 30)
2. Downloading epigenomic data
3. Developing locus-level integration scripts
4. Running enrichment analyses

## Entry Point

```bash
bash Tsai/aceEpiIntT.sh
```

## Expected Workflow

1. `download_epigenomic.sh` — Download from Synapse (requires synapse credentials)
2. `match_individuals.py` — Identify overlapping projids
3. `epi_integration.py` — Test DEG loci for differential epigenetic marks
4. `epi_visualize.py` — Generate integration figures
