# 02_Cellranger_Counts

Run Cell Ranger `count` on Tsai FASTQ files.

## Overview

This step aligns reads to the human reference genome and generates gene expression count matrices using Cell Ranger v8.0.0.

## Scripts

| Script | Description |
|--------|-------------|
| `Count_TsaiACE.ipynb` | Generate batch scripts for ACE cohort |
| `Count_TsaiResilient.ipynb` | Generate batch scripts for Resilient cohort |
| `example_count_Resilient.sh` | Example batch script (Resilient) |
| `example_count_SocIsl.sh` | Example batch script (SocIsl) |

## Prerequisites

### Cell Ranger

```bash
export PATH=/om2/user/$USER/apps/yard/cellranger-8.0.0:$PATH
```

### Reference Genome

```
/orcd/data/lhtsai/001/om2/mabdel03/yard/references/human/refdata-gex-GRCh38-2020-A
```

## Usage

### 1. Generate batch scripts

Open and run the appropriate Jupyter notebook for your cohort:
- `Count_TsaiACE.ipynb` for ACE
- `Count_TsaiResilient.ipynb` for Resilient

### 2. Submit jobs

```bash
sbatch example_count_Resilient.sh
```

## Key Parameters

```bash
cellranger count \
    --create-bam false \          # No BAM needed (no demuxlet)
    --include-introns true \      # Include intronic reads (nuclear RNA)
    --nosecondary \               # Skip secondary analysis
    --r1-length 26 \              # Trim R1 to 26bp
    --id <SAMPLE_ID> \
    --transcriptome <REF_PATH> \
    --sample <SAMPLE_ID> \
    --fastqs <FASTQ_DIRS> \       # Comma-separated if multiple
    --output-dir=<OUTPUT_DIR>
```

### Key Difference from DeJager

- `--create-bam false`: No BAM generation needed since patient assignment is known

### Multiple FASTQ Directories

Some samples have FASTQs split across directories:
```bash
--fastqs /path/to/sample_1,/path/to/sample_2
```

## Input

FASTQs located via Step 01:
```
# Locations vary by cohort, indexed in All_ROSMAP_FASTQs.csv
```

## Output

Cell Ranger outputs:
```
/om/scratch/Mon/mabdel03/Tsai/{cohort}/Counts/{projid}/
├── outs/
│   ├── raw_feature_bc_matrix.h5      # Raw counts
│   ├── filtered_feature_bc_matrix.h5 # Filtered counts
│   └── ...
└── ...
```

## Cohort Paths

| Cohort | Output Path |
|--------|-------------|
| ACE | `/om/scratch/Mon/mabdel03/Tsai/ACE/Counts/` |
| Resilient | `/om/scratch/Mon/mabdel03/Tsai/Resilient/Counts/` |
| SocIsl | `/om/scratch/Mon/mabdel03/Tsai/SocIsl/Counts/` |

## Resource Requirements

| Parameter | Value |
|-----------|-------|
| Cores | 32 |
| Memory | 128GB |
| Time | 47 hours |

## Known Issues

> **⚠️ Multiple FASTQ Directories**: Ensure comma-separated paths have no spaces.

> **⚠️ Hardcoded Paths**: Update paths in notebooks for your environment.

> **⚠️ Sample ID Matching**: The `--sample` parameter must match the FASTQ filename prefix.

