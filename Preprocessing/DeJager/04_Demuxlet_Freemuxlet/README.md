# 04_Demuxlet_Freemuxlet

Assign cells to patients using WGS genotype-based demultiplexing.

## Overview

The DeJager snRNA-seq libraries are **multiplexed**: each library contains cells
from multiple patients pooled before sequencing. This step uses Whole Genome
Sequencing (WGS) genotype data to computationally assign each cell to its patient
of origin using **Demuxlet** (reference-based demultiplexing).

### Why Demuxlet Is Needed

Unlike the Tsai dataset (where patient identity is known from metadata), the
DeJager dataset was sequenced with multiple patients per 10x channel. The only way
to determine which cell belongs to which patient is to compare the SNPs observed
in each cell's RNA reads against the patients' known WGS genotypes.

### Methods

| Method | Type | Requirement | Used In Production |
|--------|------|-------------|-------------------|
| **Demuxlet** | Reference-based | WGS VCF required | Yes |
| **Freemuxlet** | Reference-free | No genotypes needed | Testing only |

Demuxlet was chosen for production because WGS data is available for all ROSMAP
patients, providing much higher accuracy than reference-free approaches.

## Workflow

The demuxlet pipeline has three steps, each submitted as separate SLURM jobs
(or combined into one job for single libraries):

```
Cell Ranger BAM --> [A] Filter BAM --> [B] Pileup + Demuxlet --> [C] Aggregate
                        (3 hrs)              (36 hrs)              (minutes)
```

### Step A: BAM Filtering

Filters the Cell Ranger possorted BAM to keep only reads that:
- Overlap with SNP positions in the VCF file
- Have a valid cell barcode from CellBender output

This dramatically reduces BAM file size (often 10-50x), making pileup generation
tractable.

**Tool**: `popscle_helper_tools/filter_bam_file_for_popscle_dsc_pileup.sh`
**Dependencies**: samtools >= 1.10, bedtools

### Step B: Pileup + Demuxlet

Generates allele count matrices at SNP positions for each cell, then runs demuxlet
to probabilistically assign cells to patients using genotype likelihoods.

**Tool**: Demuxafy Singularity container (`Demuxafy.sif`) with popscle toolkit
**Dependencies**: Singularity

### Step C: Post-processing

Aggregates all per-library `demux1.best` files into a single CSV mapping each cell
barcode to its assigned patient.

**Tool**: `postprocess_assignments.py`

## Scripts

| Script | Description |
|--------|-------------|
| `Demuxlet_DeJager.py` | Generates and submits batch scripts for all libraries |
| `Demuxlet_DeJager.sh` | SLURM wrapper for the Python script |
| `example_demuxlet.sh` | Example script for a single library (all 3 steps) |
| `postprocess_assignments.py` | Aggregate demux results into patient assignment CSV |
| `validate_demuxlet.py` | Validate results against ground truth annotations |

### Helper Tools

The `popscle_helper_tools/` directory contains utilities for BAM and VCF
manipulation (see its own `README.md`):

| Tool | Description |
|------|-------------|
| `filter_bam_file_for_popscle_dsc_pileup.sh` | Filter BAM to SNP-overlapping reads with valid barcodes |
| `filter_vcf_file_for_popscle.sh` | VCF filtering functions (SNP-only, MAF, subsetting) |
| `sort_vcf_same_as_bam.sh` | Sort VCF contigs to match BAM header ordering |
| `popscle_dsc_pileup_merge_splitted.py` | Merge split pileup outputs |

## Prerequisites

### 1. Singularity Container

The Demuxafy container provides the popscle toolkit (pileup and demuxlet):

```bash
# Default location (see config/paths.sh):
ls "${DEMUXAFY_SIF}"

# To build from scratch:
module load ${SINGULARITY_MODULE}
singularity build Demuxafy.sif docker://drneavin/demuxafy:latest
```

### 2. WGS VCF

The SNP-only, GRCh38-lifted VCF must be available:

```bash
ls "${DEJAGER_DEMUX_VCF}"
# Expected: snp_fixedconcatenated_liftedROSMAP.vcf.gz (~84GB)
```

See `docs/VCF_PREPARATION.md` for how this VCF was prepared from per-chromosome
ROSMAP WGS files.

### 3. Patient ID Lists

Each library needs a text file listing the WGS sample IDs of patients in that
library (one SM-* ID per line):

```bash
ls "${DEJAGER_PATIENT_IDS_DIR}/individPat${LIBRARY_ID}.txt"
```

### 4. Conda Environments

```bash
source config/paths.sh
init_conda
# For BAM filtering:
conda activate "${BCFTOOLS_ENV}"
# For post-processing and batch generation:
conda activate "${PYTHON_ENV}"
```

### 5. Previous Steps Complete

- Cell Ranger BAMs: `${DEJAGER_COUNTS}/{LibraryID}/outs/possorted_genome_bam.bam`
- CellBender barcodes: `${DEJAGER_PREPROCESSED}/{LibraryID}/processed_feature_bc_matrix_cell_barcodes.csv`

## Usage

### Option 1: Run all libraries (batch)

```bash
source config/paths.sh
sbatch Demuxlet_DeJager.sh
```

This generates per-library scripts and submits them to SLURM.

### Option 2: Run a single library

Edit `LIBRARY_ID` in `example_demuxlet.sh`, then:

```bash
source config/paths.sh
sbatch example_demuxlet.sh
```

### Option 3: Step-by-step

```bash
source config/paths.sh

# Generate BAM filter scripts only:
python Demuxlet_DeJager.py --bam-only --submit

# After BAM jobs complete, generate demuxlet scripts:
python Demuxlet_DeJager.py --demux-only --submit

# After demuxlet jobs complete, aggregate results:
python postprocess_assignments.py \
    --wgs-dir "${DEJAGER_WGS_DIR}" \
    --output cell_to_patient_assignments.csv
```

### Option 4: Generate without submitting (dry run)

```bash
python Demuxlet_DeJager.py --all
# Scripts are written to ${DEJAGER_WGS_DIR}/generated_scripts/ but not submitted
```

## Key Parameters

| Parameter | Value | Rationale |
|-----------|-------|-----------|
| `--alpha` | 0.05 | Prior for ambient RNA fraction; balances sensitivity vs specificity |
| `--min-mac` | 1 | Minimum minor allele count; using 5 drops accuracy to 0.61 |
| `--doublet-prior` | 0.1 | Prior probability of doublets; 10% matches 10x expected rate |
| `--field` | "PL" | Use Phred-scaled genotype likelihoods from WGS VCF |

See `docs/PARAMETER_TUNING.md` for the full parameter sweep results and rationale.

## Input

| Input | Source | Path |
|-------|--------|------|
| Possorted BAM | Step 02 (Cell Ranger) | `${DEJAGER_COUNTS}/{LibID}/outs/possorted_genome_bam.bam` |
| Cell barcodes | Step 03 (CellBender) | `${DEJAGER_PREPROCESSED}/{LibID}/processed_feature_bc_matrix_cell_barcodes.csv` |
| SNP VCF | WGS data (external) | `${DEJAGER_DEMUX_VCF}` |
| Patient IDs | WGS metadata (external) | `${DEJAGER_PATIENT_IDS_DIR}/individPat{LibID}.txt` |

## Output

### Per-library (in `${DEJAGER_WGS_DIR}/{LibraryID}/`)

```
BAMOutput1.bam      # Filtered BAM (Step A output)
BAMOutput1.bam.csi  # BAM index
plpDemux1.cel.gz    # Pileup cell file
plpDemux1.plp.gz    # Pileup data (allele counts per SNP per cell)
plpDemux1.var.gz    # Pileup variant information
demux1.best         # Cell-patient assignments (Step B output)
```

### Aggregated (Step C output)

```
cell_to_patient_assignments.csv
    Columns: Cell Barcode, Assigned Patient, Library
```

### Output Format (demux1.best)

| Column | Description |
|--------|-------------|
| BARCODE | Cell barcode |
| NUM.SNPS | Number of SNPs covering this cell |
| NUM.READS | Number of reads at SNP positions |
| DROPLET.TYPE | SNG (singlet), DBL (doublet), AMB (ambiguous) |
| BEST.GUESS | Best patient assignment (comma-separated pair if doublet) |
| BEST.LLK | Log-likelihood for best assignment |
| DIFF.LLK.BEST.NEXT | Log-likelihood difference between best and next-best |

## Resource Requirements

| Step | Cores | Memory | Time | Scope |
|------|-------|--------|------|-------|
| A. BAM Filter | 45 | 400GB | 3h | Per library |
| B. Pileup + Demuxlet | 10 | 500GB | 36h | Per library |
| C. Post-processing | 1 | 8GB | 5min | All libraries |
| Batch generator | 1 | 4GB | <1h | Once |

## Validation

To validate demuxlet results against known cell annotations (if available):

```bash
python validate_demuxlet.py \
    --demux-file ${DEJAGER_WGS_DIR}/191121-B6/demux1.best \
    --annotation-file cell-annotation.csv \
    --library-id 191121-B6 \
    --output-dir results/
```

This computes ARI, NMI, and generates a confusion matrix heatmap.

## Known Issues

> **Singularity Bind Mounts**: The `--bind` flag maps the entire WGS directory
> to `/mnt` inside the container. Cell barcode files must be accessible at
> `/mnt/processed_feature_bc_matrix_cell_barcodes_{LibID}.csv`.

> **Memory**: Pileup generation on large libraries can approach 500GB. Monitor
> jobs and increase `--mem` if needed.

> **"alone" Libraries**: Libraries containing "alone" in their name (e.g.,
> `191122-B6-R1969233-alone`) are single-patient libraries created by splitting.
> They are automatically excluded from batch processing. Their patient IDs are
> handled via `Processing/DeJager/Pipeline/Resources/patient_id_overrides.json`.

> **Path Configuration**: All paths are configured in `config/paths.sh`. Run
> `source config/paths.sh && check_paths` to verify your setup.

## Troubleshooting

### No cells assigned

- Check VCF format (must have PL field): `bcftools query -l your.vcf.gz`
- Verify patient IDs match between VCF sample names and `individPat*.txt`
- Ensure barcodes are in correct format (no header, one per line)

### Too many doublets

- Check library quality metrics from Cell Ranger
- Consider lowering `--doublet-prior`

### Low assignment confidence

- Check BAM/VCF chromosome naming consistency
- Use the SNP-only VCF (not the full variant VCF)
- Verify VCF is properly indexed: `tabix -p vcf your.vcf.gz`

### Pileup runs out of memory

- Ensure BAM filtering (Step A) completed successfully
- Increase `--mem` in the SLURM script
- Check filtered BAM size; if still very large, the VCF may have too many variants

## Documentation

- `docs/PARAMETER_TUNING.md` - Full parameter sweep results with accuracy metrics
- `docs/VCF_PREPARATION.md` - How the WGS VCF was prepared (liftover, filtering, indexing)
- `popscle_helper_tools/README.md` - Helper tool documentation
- `archive/` - Original experimental scripts (preserved for reference)
