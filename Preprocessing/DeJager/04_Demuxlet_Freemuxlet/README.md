# 04_Demuxlet_Freemuxlet

Assign cells to patients using genotype-based demultiplexing.

## Overview

The DeJager libraries contain cells from multiple patients (multiplexed). This step uses WGS genotype data to assign each cell to its patient of origin using Demuxlet or Freemuxlet.

## Methods

| Method | Description | Requirement |
|--------|-------------|-------------|
| **Demuxlet** | Reference-based demultiplexing | Requires WGS VCF |
| **Freemuxlet** | Reference-free demultiplexing | No genotype data needed |

## Scripts

| Script | Description |
|--------|-------------|
| `demuxTest.sh` | Main demuxlet workflow with parameter testing |
| `demux.sh` | Simplified demuxlet script |
| `freemux.sh` | Freemuxlet workflow |
| `pileup.sh` | Generate pileup files from BAM |
| `demuxValCode.py` | Validate demuxlet results |
| `freemuxletValidateCode.py` | Validate freemuxlet results |

## Prerequisites

### Singularity Container

```bash
module load openmind/singularity/3.10.4
# Demuxafy container at: /path/to/Demuxafy.sif
```

### Input Files

1. **BAM file**: Possorted BAM from Cell Ranger
2. **VCF file**: Filtered WGS variants for patients in the library
3. **Barcodes**: Cell barcodes from CellBender output
4. **Sample list**: Patient IDs expected in the library

## Workflow

### 1. Generate Pileup

```bash
singularity exec Demuxafy.sif popscle_pileup.py \
    --sam <BAM> \
    --vcf <VCF> \
    --group-list <BARCODES> \
    --out <OUTPUT_PREFIX> \
    --sm-list <PATIENT_IDS>
```

### 2. Run Demuxlet

```bash
singularity exec Demuxafy.sif popscle demuxlet \
    --plp <PILEUP_PREFIX> \
    --vcf <VCF> \
    --field "PL" \
    --group-list <BARCODES> \
    --sm-list <PATIENT_IDS> \
    --out <OUTPUT_PREFIX> \
    --alpha 0.05 \
    --min-mac 1 \
    --doublet-prior 0.1
```

### 3. Validate Results

```bash
python demuxValCode.py --input <DEMUXLET_OUTPUT>
```

## Key Parameters

| Parameter | Description | Recommended |
|-----------|-------------|-------------|
| `--alpha` | Prior for ambient RNA | 0.05 |
| `--min-mac` | Minimum minor allele count | 1 |
| `--doublet-prior` | Prior probability of doublets | 0.1 |
| `--field` | Genotype field in VCF | "PL" |

## Output

```
<OUTPUT_PREFIX>.best    # Cell-patient assignments
<OUTPUT_PREFIX>.single  # Singlet calls
<OUTPUT_PREFIX>.sing2   # Alternative singlet format
```

### Output Format (.best)

| Column | Description |
|--------|-------------|
| BARCODE | Cell barcode |
| BEST | Best assignment (SNG/DBL/AMB) |
| SNG.1ST | Most likely patient |
| SNG.LLK1 | Log-likelihood for best assignment |

## Resource Requirements

| Parameter | Value |
|-----------|-------|
| Cores | 80 |
| Memory | 400GB |
| Time | 48 hours |

## Known Issues

> **⚠️ demuxTest.sh Cleanup**: This script contains extensive commented-out parameter tuning experiments (lines 1-52, 70-151). The active commands are at lines 53-54 and 69. Consider cleaning up for production use.

> **⚠️ VCF Preparation**: The VCF must be filtered for common variants (MAF > 0.01) and indexed with tabix.

> **⚠️ Memory Usage**: Pileup generation is memory-intensive. Use 400GB+ for large libraries.

> **⚠️ Singularity Binds**: Ensure proper bind mounts for data directories.

## VCF Preparation

The WGS VCF should be prepared as follows:

1. Subset to patients in the library
2. Filter for common variants (MAF > 0.01)
3. Remove multiallelic sites
4. Index with tabix

```bash
bcftools view -s <PATIENT_IDS> input.vcf.gz | \
    bcftools filter -i 'MAF>0.01' | \
    bcftools norm -m -any | \
    bgzip > filtered.vcf.gz
tabix filtered.vcf.gz
```

## Troubleshooting

### No cells assigned

- Check VCF format (must have PL field)
- Verify patient IDs match between VCF and sample list
- Ensure barcodes are in correct format

### Too many doublets

- Lower `--doublet-prior`
- Check library quality

### Low assignment rate

- Use more SNPs (lower MAF filter)
- Check BAM/VCF chromosome naming consistency

