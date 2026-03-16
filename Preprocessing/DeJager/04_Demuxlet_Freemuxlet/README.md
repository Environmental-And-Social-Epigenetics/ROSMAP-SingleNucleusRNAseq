# 04_Demuxlet_Freemuxlet

Assign cells to patients using WGS genotype-based demultiplexing.

> **New to this project?** See `BACKGROUND.md` at the repository root for context
> on the ROSMAP study, multiplexed snRNA-seq, and where this step fits in the
> overall pipeline.

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

### References

- **Demuxlet paper**: Kang et al. (2018) "Multiplexed droplet single-cell
  RNA-sequencing using natural genetic variation." *Nature Biotechnology* 36:89-94.
  https://doi.org/10.1038/nbt.4042
- **popscle toolkit** (implements demuxlet/freemuxlet): https://github.com/statgen/popscle
- **Demuxafy** (container with popscle + helper tools): https://demuxafy.readthedocs.io/
  Docker image: `drneavin/demuxafy:latest`

## Data Access and Provenance

### ROSMAP Data Access

ROSMAP data is controlled-access and requires a Data Use Agreement (DUA) through
the AD Knowledge Portal on Synapse (https://adknowledgeportal.synapse.org/).
To access any ROSMAP data:

1. Create a Synapse account at synapse.org
2. Complete a Data Use Agreement for the relevant ROSMAP study
3. Request access to the specific Synapse datasets

If you have not previously set up Synapse credentials on this cluster, see
`Preprocessing/DeJager/01_FASTQ_Download/README.md` for authentication setup.

### Where the WGS Genotype Data Comes From

The Whole Genome Sequencing (WGS) VCF used for demultiplexing was generated as
part of the ROSMAP genomics collection on Synapse (study: `syn3219045`).
The WGS variant calls are available under the ROSMAP WGS data on the
AD Knowledge Portal — navigate to the study page and look for the WGS
per-chromosome VCFs. The original files are per-chromosome VCFs aligned to
GRCh37 (hg19):

```
DEJ_11898_B01_GRM_WGS_2017-05-15_{chr}.recalibrated_variants_chr.vcf.gz
```

These were processed through a 5-step pipeline (liftover to GRCh38,
chromosome name fixing, concatenation, SNP filtering, indexing) to produce the
production VCF at `${DEJAGER_DEMUX_VCF}` (~84GB). See `docs/VCF_PREPARATION.md`
for the full preparation pipeline.

### Where Patient ID Files Come From

The per-library patient ID files (`individPat{LibID}.txt`) were created by
cross-referencing:
- The **library pooling manifest** from the DeJager lab, which records which
  patients were loaded into each 10x Chromium channel
- The **WGS sample manifest**, which maps WGS sample identifiers (SM-* format)
  to ROSMAP participant IDs

These manifests are lab-internal records provided with the dataset metadata.
The resulting individPat files are stored at `${DEJAGER_PATIENT_IDS_DIR}/`.

### Patient Identifier Systems

This pipeline uses several identifier systems. Understanding the mapping is
important for connecting demuxlet output to clinical metadata:

| ID Type | Example | Where Used |
|---------|---------|------------|
| SM-* (WGS sample ID) | SM-CJGLI | VCF header, individPat files, demuxlet output |
| projid (ROSMAP participant ID) | 482428 | Clinical phenotype data, downstream analysis |
| individualID | R3978789 | Some metadata files |

**Mapping between IDs**:
- SM-* to projid: handled by the aggregated assignment CSV (`${DEJAGER_PATIENT_MAP}`)
- projid to individualID: `Data/Phenotypes/DeJager_ID_Map.csv` (`${DEJAGER_ID_MAP}`)

The SM-* to projid conversion happens during downstream Processing (Stage 1 QC
filtering), not during this preprocessing step. Demuxlet outputs use SM-* IDs.

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

All path variables below (e.g., `${DEMUXAFY_SIF}`, `${DEJAGER_DEMUX_VCF}`) are
defined in `config/paths.sh` at the repository root. Source it before any
operation:

```bash
source config/paths.sh
```

If you need to customize paths for your environment (e.g., different cluster),
create `config/paths.local.sh` (gitignored) with your overrides.

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

The BAM filtering step requires `samtools` and `bedtools`; the batch generation
and post-processing steps require Python 3 with standard libraries.

```bash
source config/paths.sh
init_conda

# For BAM filtering (needs samtools >= 1.10, bedtools):
conda activate "${BCFTOOLS_ENV}"
# If creating from scratch:
#   conda create -n bcftools_env -c bioconda -c conda-forge samtools>=1.10 bedtools bcftools

# For post-processing and batch generation (needs Python 3):
conda activate "${PYTHON_ENV}"
# If creating from scratch:
#   conda create -n python_env python=3.9
```

### 5. Previous Steps Complete

- Cell Ranger BAMs: `${DEJAGER_COUNTS}/{LibraryID}/outs/possorted_genome_bam.bam`
- CellBender barcodes: `${DEJAGER_PREPROCESSED}/{LibraryID}/processed_feature_bc_matrix_cell_barcodes.csv`

### 6. Barcode Files in WGS Directory

The Singularity container bind-mounts `${DEJAGER_WGS_DIR}` as `/mnt`, so the
pileup and demuxlet commands expect barcode files to be present **in the WGS
directory root** with a specific naming convention:

```
${DEJAGER_WGS_DIR}/processed_feature_bc_matrix_cell_barcodes_{LibraryID}.csv
```

You must copy (or symlink) each library's CellBender barcode file there before
running Step B:

```bash
for LIB in $(ls "${DEJAGER_PATIENT_IDS_DIR}" | sed 's/individPat//;s/.txt//'); do
    cp "${DEJAGER_PREPROCESSED}/${LIB}/processed_feature_bc_matrix_cell_barcodes.csv" \
       "${DEJAGER_WGS_DIR}/processed_feature_bc_matrix_cell_barcodes_${LIB}.csv"
done
```

The batch script (`Demuxlet_DeJager.py`) does **not** do this automatically —
it assumes the barcode files are already in place.

### 7. Preflight Check (Recommended)

Before running any demuxlet scripts, verify that all prerequisites are in place:

```bash
source config/paths.sh
bash config/preflight.sh dejager-04-demuxlet
```

This checks that config variables are set, conda environments exist, the
Singularity module is loadable, the Demuxafy container exists, the WGS VCF is
accessible, patient ID files are present, and helper tool scripts are in place.
Fix any `FAIL` items before proceeding.

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

### Example Output

A typical `demux1.best` file (tab-separated):

```
BARCODE               NUM.SNPS  NUM.READS  DROPLET.TYPE  BEST.GUESS          BEST.LLK    DIFF.LLK.BEST.NEXT
AAACCCAAGTCGTACT-1    847       1203       SNG           SM-CJGLI,SM-CJGLI   -312.45     87.23
AAACCCATCAGTTGAC-1    523       689        DBL           SM-CTDRI,SM-CJGGZ   -287.11     12.56
AAACGAACAGCTGTAT-1    34        41         AMB           SM-CJJ1R,SM-CJJ1R   -89.02      1.34
```

- **SNG** (singlet): cell assigned to a single patient (BEST.GUESS shows same ID twice)
- **DBL** (doublet): droplet contains cells from two patients (BEST.GUESS shows the pair)
- **AMB** (ambiguous): assignment is uncertain

### Interpreting DIFF.LLK.BEST.NEXT

`DIFF.LLK.BEST.NEXT` is the log-likelihood ratio between the best and
second-best patient assignment. Higher values indicate more confident
assignments:

- **> 20**: High confidence -- strong genotype match, unambiguous assignment
- **5-20**: Moderate confidence -- usually correct but worth noting
- **< 5**: Low confidence -- the cell may be mislabeled; often flagged as AMB

This pipeline does NOT filter on DIFF.LLK.BEST.NEXT. All cells (including AMB)
are passed through to downstream processing, where low-quality cells are
removed during QC filtering (Processing Stage 1).

### Expected Distribution

For a typical well-prepared multiplexed library with 6-8 patients:
- **Singlets (SNG)**: 80-90% of droplets
- **Doublets (DBL)**: 5-15% (depends on loading concentration)
- **Ambiguous (AMB)**: 1-5%

If you see a very high doublet rate (>20%) or high AMB rate (>10%), check:
- Library quality metrics from Cell Ranger (Step 02)
- VCF/BAM chromosome naming consistency
- Whether the correct patient ID file was used for the library

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

> **"alone" Libraries**: Some DeJager libraries have names containing "alone"
> (e.g., `191122-B6-R1969233-alone`). These are **single-patient libraries**
> where a patient's cells were isolated into a separate 10x channel rather than
> being pooled with other patients. Because there is only one patient, demuxlet
> cannot and need not be run -- there is nothing to demultiplex.
>
> These libraries are **automatically excluded** from batch processing by
> `Demuxlet_DeJager.py` (the `discover_libraries()` function skips any library
> ID containing "alone") and by `postprocess_assignments.py`.
>
> Their patient identity is known from the library name itself: the R-number
> in the name (e.g., R1969233) is mapped to a projid via the overrides file at
> `Processing/DeJager/Pipeline/Resources/patient_id_overrides.json`. This
> mapping is applied during Processing Stage 1 (QC filtering), not during this
> preprocessing step.
>
> **You do not need to do anything special for "alone" libraries** -- they are
> handled automatically.

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
