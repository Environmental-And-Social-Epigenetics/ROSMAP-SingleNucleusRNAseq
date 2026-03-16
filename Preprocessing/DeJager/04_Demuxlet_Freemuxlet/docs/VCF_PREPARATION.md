# VCF Preparation for Demuxlet

This document describes how the WGS VCF file used for demuxlet was prepared.
The final production VCF is `snp_fixedconcatenated_liftedROSMAP.vcf.gz` (~84GB).

## Source Data

The ROSMAP Whole Genome Sequencing (WGS) data was obtained from the
AD Knowledge Portal on Synapse (https://adknowledgeportal.synapse.org/) as part
of the ROSMAP genomics collection. These are per-chromosome VCF files containing
variants for all sequenced ROSMAP subjects.

**Access requirements**: ROSMAP WGS data is controlled-access. You must have a
Synapse account and an approved Data Use Agreement (DUA) for the ROSMAP study.
Navigate to the ROSMAP study page on the AD Knowledge Portal to locate the WGS
VCF files.

Original files (per chromosome):
```
DEJ_11898_B01_GRM_WGS_2017-05-15_{chr}.recalibrated_variants_chr.vcf.gz
```

## Preparation Pipeline

### 1. Liftover from GRCh37 to GRCh38

The original WGS VCFs were aligned to GRCh37 (hg19). Since Cell Ranger aligns
reads to GRCh38 (hg38), the VCFs needed to be lifted over to match.

**Tools**: Either Picard LiftoverVcf or CrossMap can be used. Example with Picard:

```bash
java -jar picard.jar LiftoverVcf \
    I=input_hg19.vcf.gz \
    O=output_hg38.vcf.gz \
    CHAIN=hg19ToHg38.over.chain.gz \
    REJECT=rejected_variants.vcf.gz \
    R=hg38.fa
```

**Chain file**: `hg19ToHg38.over.chain.gz` is available from the UCSC Genome
Browser downloads (https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/).

**Reference FASTA**: The GRCh38 reference must match the Cell Ranger reference.
If using the 10x-provided reference (`refdata-gex-GRCh38-2020-A`), the FASTA
is at `fasta/genome.fa` within the reference directory.

Some variants will fail liftover (typically <1%). These are written to the
REJECT file and excluded from the final VCF.

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

These files were created by cross-referencing:
- The **library pooling records** from the DeJager lab, which document which
  patients (by WGS sample ID) were loaded into each 10x Chromium channel.
  These records are lab-internal and were provided with the dataset metadata.
- The **WGS sample manifest**, which maps SM-* identifiers to ROSMAP
  participant IDs.

The SM-* IDs in these files must exactly match the sample names in the VCF
header. You can verify this with:

```bash
bcftools query -l ${DEJAGER_DEMUX_VCF} | grep "SM-CJGLI"
```

If you are applying this pipeline to a different multiplexed dataset, you will
need to create equivalent patient ID files for your libraries by determining
which patients were pooled in each library.

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
