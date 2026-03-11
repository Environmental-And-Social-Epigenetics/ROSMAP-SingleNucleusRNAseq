# Demuxlet Parameter Tuning Results

This document records the parameter sweep performed during demuxlet optimization.
The experiments used library **191121-B6** as a test case, comparing demuxlet
assignments against known cell annotations to compute accuracy (Adjusted Rand
Index).

## Test Setup

- **Library**: 191121-B6 (8 patients)
- **Pileup**: plpDemux10 with `filteredFixedVCF.vcf.gz`
- **Barcodes**: CellBender filtered barcodes
- **Metric**: Adjusted Rand Index (ARI) against known cell-annotation.csv

## Parameter Sweep Results

| Run | geno-error-offset | doublet-prior | cap-BQ | alpha | min-mac | min-total | min-callrate | ARI |
|-----|-------------------|---------------|--------|-------|---------|-----------|--------------|-----|
| Baseline | 0.2 | 0.35 | - | default | default | - | - | 0.841 |
| 2 | 0.2 | 0.25 | - | default | default | - | - | 0.841 |
| 3 | 0.3 | 0.35 | - | default | default | - | - | 0.857 |
| 4 | 0.45 | 0.35 | - | default | default | - | - | 0.863 |
| 5 | 0.4 | 0.35 | - | default | default | - | - | 0.862 |
| 6 | 0.5 | 0.35 | - | default | default | - | - | 0.863 |
| 7a | 0.45 | 0.35 | - | default | default | 1000 | - | 0.886 |
| 7b | 0.45 | 0.35 | - | default | default | 1200 | - | 0.916 |
| 8a | 0.5 | 0.35 | 15 | default | default | - | - | 0.864 |
| 8b | 0.5 | 0.35 | 15 | default | default | 1400 | - | 0.945 |
| 8c | 0.5 | 0.35 | 15 | 0.01 | default | - | - | 0.864 |
| 8d | 0.5 | 0.35 | 15 | 0.8 | default | - | - | 0.823 |
| 8e | 0.5 | 0.35 | 15 | 0.1 | default | - | - | 0.864 |
| 8f | 0.5 | 0.35 | 15 | 0.1 | 5 | - | - | 0.610 |
| 8g | 0.5 | 0.35 | 15 | 0.1 | default | - | 0.4 | 0.864 |
| 8h | 0.5 | 0.35 | 15 | 0.1 | default | 900 | - | 0.874 |
| 9a | 0.25 | 0.35 | 15 | 0.1 | 1 | 900 | - | 0.874 |
| 9b | default | 0.1 | - | 0.05 | 1 | - | - | 0.857 |
| 9c | 0.1 | 0.15 | - | 0.05 | 1 | - | - | - |

## Key Findings

1. **`--min-total` has the largest positive effect**: Increasing the minimum total
   read depth from default to 1000-1400 dramatically improved accuracy (0.841 ->
   0.945). However, this drops cells with low read counts, so it trades coverage
   for accuracy.

2. **`--min-mac 5` is catastrophic**: Setting minimum minor allele count to 5
   dropped accuracy from 0.864 to 0.610. Keep `--min-mac 1`.

3. **`--alpha` (ambient RNA prior) has moderate effect**: Very high values (0.8)
   reduce accuracy. Values between 0.01-0.1 perform similarly. Default behavior
   and `--alpha 0.05` are both good.

4. **`--geno-error-offset` and `--cap-BQ` have marginal effects**: Increasing
   geno-error-offset from 0.2 to 0.5 improved accuracy slightly (0.841 -> 0.863).
   Adding `--cap-BQ 15` had negligible additional effect.

5. **`--min-callrate`** had no measurable effect at 0.4.

6. **VCF quality matters more than most parameters**: Switching from the initial
   VCF to the SNP-only filtered version (`filteredSNPFreq2.vcf.gz`) with proper
   indexing had a larger effect than any single parameter change.

## Production Parameters

After tuning, the following parameters were chosen for production processing of
all 60+ libraries:

```bash
popscle demuxlet \
    --alpha 0.05 \
    --min-mac 1 \
    --doublet-prior 0.1 \
    --field "PL"
```

**Rationale**: These parameters prioritize retaining all cells (no `--min-total`
filter) while using conservative ambient RNA and doublet priors. The `--min-mac 1`
keeps all informative SNPs. The `--doublet-prior 0.1` (10%) matches the expected
doublet rate for 10x Chromium experiments. The `--field "PL"` uses Phred-scaled
genotype likelihoods from the WGS VCF.

The `--min-total` filter was not used in production because dropping low-coverage
cells would bias the dataset. Instead, low-confidence assignments are handled
downstream during QC filtering in the Processing pipeline (Stage 1).

## Original Experiments

The raw parameter experiments are preserved in `archive/demuxTest.sh` (lines 86-125).
