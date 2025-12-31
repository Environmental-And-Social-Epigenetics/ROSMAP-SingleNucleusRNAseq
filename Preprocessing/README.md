# Preprocessing

This directory contains scripts for the initial processing of raw snRNA-seq data, from FASTQ files through ambient RNA correction.

## Overview

The preprocessing phase converts raw sequencing data into cleaned, ambient RNA-corrected count matrices suitable for downstream analysis.

```
FASTQs → Cell Ranger (alignment/counting) → CellBender (ambient RNA removal) → Clean matrices
```

## Directory Structure

```
Preprocessing/
├── DeJager/
│   ├── 01_FASTQ_Download/    # Download FASTQs from Synapse
│   ├── 02_Cellranger_Counts/ # Run Cell Ranger count
│   ├── 03_Cellbender/        # Remove ambient RNA
│   └── 04_Demuxlet_Freemuxlet/ # Assign cells to patients
└── Tsai/
    ├── 01_FASTQ_Location/    # Locate FASTQs on Engaging
    ├── 02_Cellranger_Counts/ # Run Cell Ranger count
    └── 03_Cellbender/        # Remove ambient RNA
```

## Dataset Differences

### DeJager Dataset

- **Source**: Downloaded from Synapse (syn21438684)
- **Sample Assignment**: Requires Demuxlet/Freemuxlet using WGS data to assign cells to patients (multiplexed libraries)
- **Extra Step**: Step 04 (Demuxlet/Freemuxlet) is unique to DeJager

### Tsai Dataset

- **Source**: Located on MIT Engaging cluster (`/nfs/picower*`)
- **Sample Assignment**: Patient assignments are known from sequencing metadata
- **Cohorts**: ACE, Resilient, SocIsl (Social Isolation)

## Workflow

### 1. Data Acquisition

| Dataset | Method |
|---------|--------|
| DeJager | `synapse get` via Python SDK |
| Tsai | Enumerate FASTQs from known paths |

### 2. Cell Ranger Count

Alignment and counting with Cell Ranger v8.0.0:

```bash
cellranger count \
    --include-introns true \
    --nosecondary \
    --r1-length 26 \
    --transcriptome=/orcd/data/lhtsai/001/om2/mabdel03/yard/references/human/refdata-gex-GRCh38-2020-A \
    --sample <SAMPLE_ID> \
    --fastqs <FASTQ_DIR> \
    --output-dir=<OUTPUT_DIR>
```

Key flags:
- `--include-introns true`: Include intronic reads (important for nuclear RNA)
- `--nosecondary`: Skip secondary analysis to save time
- `--r1-length 26`: Trim R1 to 26bp (library-specific)

### 3. CellBender

Ambient RNA removal using CellBender:

```bash
cellbender remove-background \
    --cuda \
    --input <RAW_MATRIX.h5> \
    --fpr 0 \
    --output <OUTPUT.h5>
```

Key flags:
- `--cuda`: Use GPU acceleration
- `--fpr 0`: False positive rate (stringent setting)

### 4. Demuxlet/Freemuxlet (DeJager only)

Assign cells to patients using genotype data:

1. Generate pileup from BAM and VCF
2. Run demuxlet (with WGS) or freemuxlet (without WGS)
3. Validate assignments

## Outputs

Each preprocessing step produces:

| Step | Output |
|------|--------|
| Cell Ranger | `outs/raw_feature_bc_matrix.h5`, `outs/filtered_feature_bc_matrix.h5`, BAM |
| CellBender | `processed_feature_bc_matrix.h5` |
| Demuxlet | `.best` file with cell-patient assignments |

## Next Steps

After preprocessing, proceed to the `Processing/` directory for QC and batch correction.

