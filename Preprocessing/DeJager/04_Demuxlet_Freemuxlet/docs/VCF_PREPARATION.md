# VCF Preparation for Demuxlet

This document describes how the WGS VCF file used for demuxlet was prepared.
The final production VCF is `snp_fixedconcatenated_liftedROSMAP.vcf.gz` (~84GB).

## Source Data

The ROSMAP Whole Genome Sequencing (WGS) data was obtained as per-chromosome VCF
files from the ROSMAP study. These files contained variants for all sequenced
ROSMAP subjects.

Original files (per chromosome):
```
DEJ_11898_B01_GRM_WGS_2017-05-15_{chr}.recalibrated_variants_chr.vcf.gz
```

## Preparation Pipeline

### 1. Liftover from GRCh37 to GRCh38

The original WGS VCFs were aligned to GRCh37 (hg19). Since Cell Ranger aligns
reads to GRCh38 (hg38), the VCFs needed to be lifted over to match.

Tools: Picard LiftoverVcf or CrossMap

### 2. Chromosome Name Fixing

After liftover, chromosome names were standardized to match the BAM file contig
naming convention (e.g., `chr1`, `chr2`, ..., `chrX`, `chrY`).

The `popscle_helper_tools/sort_vcf_same_as_bam.sh` utility can be used to ensure
VCF contig ordering matches the BAM header.

### 3. Concatenation

Per-chromosome VCFs were concatenated into a single file:

```bash
bcftools concat \
    chr1_lifted.vcf.gz chr2_lifted.vcf.gz ... chrY_lifted.vcf.gz \
    -Oz -o concatenated_liftedROSMAP.vcf.gz
tabix concatenated_liftedROSMAP.vcf.gz
```

Output: `concatenated_liftedROSMAP.vcf.gz` (~276GB)

### 4. SNP-Only Filtering

For demuxlet, only SNPs are informative (not indels or structural variants).
Filtering to SNPs dramatically reduces file size:

```bash
bcftools view -v snps \
    fixedconcatenated_liftedROSMAP.vcf.gz \
    -Oz -o snp_fixedconcatenated_liftedROSMAP.vcf.gz
tabix snp_fixedconcatenated_liftedROSMAP.vcf.gz
```

Output: `snp_fixedconcatenated_liftedROSMAP.vcf.gz` (~84GB)

### 5. Indexing

The final VCF must be indexed with tabix for random access:

```bash
tabix -p vcf snp_fixedconcatenated_liftedROSMAP.vcf.gz
```

## Per-Library Patient ID Files

For each library, a text file lists the WGS sample IDs (SM-* format) of the
patients pooled in that library. These are located at:

```
${DEJAGER_PATIENT_IDS_DIR}/individPat{LibraryID}.txt
```

Example (`individPat190403-B4-A.txt`):
```
SM-CJGLI
SM-CTDRI
SM-CJGGZ
SM-CJJ1R
SM-CJIYJ
SM-CTEEG
```

These files were created by cross-referencing the library pooling records with
the WGS sample manifest.

## Additional VCF Filtering (Optional)

The `popscle_helper_tools/filter_vcf_file_for_popscle.sh` script provides
functions for additional VCF filtering:

- **Subset to specific patients**: `subset_samples_from_vcf`
- **Filter by MAF**: useful for keeping only common variants
- **Remove missing genotypes**: `filter_out_mutations_missing_genotype_for_one_or_more_samples`
- **Recalculate AF/AC/AN**: `calculate_AF_AC_AN_values_based_on_genotype_info`

These filters were applied during early testing but the production workflow uses
the full SNP VCF with per-library patient subsetting handled by the `--sm-list`
argument to popscle.

## File Sizes

| File | Size | Description |
|------|------|-------------|
| `concatenated_liftedROSMAP.vcf.gz` | ~276GB | All variants, all chromosomes |
| `fixedconcatenated_liftedROSMAP.vcf.gz` | ~275GB | Chromosome names fixed |
| `snp_fixedconcatenated_liftedROSMAP.vcf.gz` | ~84GB | **Production VCF** (SNPs only) |
