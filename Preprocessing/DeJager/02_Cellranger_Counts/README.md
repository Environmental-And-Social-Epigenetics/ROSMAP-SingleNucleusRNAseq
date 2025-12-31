# 02_Cellranger_Counts

Run Cell Ranger `count` on DeJager FASTQ files.

## Overview

This step aligns reads to the human reference genome and generates gene expression count matrices using Cell Ranger v8.0.0.

## Scripts

| Script | Description |
|--------|-------------|
| `Count_DeJager.py` | Generates batch scripts for all libraries and submits them |
| `Count_DeJager.sh` | SLURM wrapper for the Python script |
| `example_count.sh` | Example batch script for a single library |

## Prerequisites

### Cell Ranger

```bash
export PATH=/orcd/data/lhtsai/001/om2/mabdel03/apps/yard/cellranger-8.0.0:$PATH
```

### Reference Genome

```
/orcd/data/lhtsai/001/om2/mabdel03/yard/references/human/refdata-gex-GRCh38-2020-A
```

## Usage

### Option 1: Generate and submit all jobs

```bash
python Count_DeJager.py
```

This script:
1. Reads library IDs from the FASTQ CSV
2. Creates output directories
3. Generates batch scripts for each library
4. Submits jobs to SLURM

### Option 2: Submit single library

```bash
sbatch example_count.sh
```

## Key Parameters

```bash
cellranger count \
    --create-bam true \           # Generate BAM (needed for demuxlet)
    --include-introns true \      # Include intronic reads (nuclear RNA)
    --nosecondary \               # Skip secondary analysis
    --r1-length 26 \              # Trim R1 to 26bp
    --id <LIBRARY_ID> \
    --transcriptome <REF_PATH> \
    --sample <SAMPLE_ID> \
    --fastqs <FASTQ_DIR> \
    --output-dir=<OUTPUT_DIR>
```

## Input

FASTQs from Step 01:
```
/om/scratch/Mon/mabdel03/FASTQs/{LibraryID}/
```

## Output

Cell Ranger outputs:
```
/om/scratch/Mon/mabdel03/Counts/{LibraryID}/
├── outs/
│   ├── raw_feature_bc_matrix.h5      # Raw counts
│   ├── filtered_feature_bc_matrix.h5 # Filtered counts
│   ├── possorted_genome_bam.bam      # Aligned reads
│   └── ...
└── ...
```

## Resource Requirements

| Parameter | Value |
|-----------|-------|
| Cores | 32 |
| Memory | 128GB |
| Time | 47 hours |

## Known Issues

> **⚠️ Hardcoded Paths**: The `Count_DeJager.py` script has hardcoded input/output paths. Update lines 10-19 for your environment.

> **⚠️ Sample Naming**: Some libraries have two samples (A and B lanes). The script handles this automatically.

> **⚠️ Scratch Space**: Cell Ranger generates large intermediate files. Ensure sufficient scratch space.

