#!/bin/bash
#
# Example demuxlet workflow for a single DeJager library.
#
# Performs the full 3-step workflow:
#   A) Filter BAM to reads overlapping SNPs with valid cell barcodes
#   B) Generate pileup at SNP positions
#   C) Run demuxlet to assign cells to patients
#
# Paths are resolved from config/paths.sh environment variables.
# Edit LIBRARY_ID below to process a different library.
#
#SBATCH -n 45
#SBATCH -t 40:00:00
#SBATCH --mem=500G
#SBATCH --mail-type=FAIL

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

# ---- Configuration ----
LIBRARY_ID="190403-B4-A"

OUT_DIR="${DEJAGER_WGS_DIR}/${LIBRARY_ID}"
BARCODES="${DEJAGER_PREPROCESSED}/${LIBRARY_ID}/processed_feature_bc_matrix_cell_barcodes.csv"
PATIENT_IDS="${DEJAGER_PATIENT_IDS_DIR}/individPat${LIBRARY_ID}.txt"
VCF_BASENAME="$(basename "${DEJAGER_DEMUX_VCF}")"

mkdir -p "${OUT_DIR}"

# ===========================================================================
# Step A: Filter BAM
# ===========================================================================
# Keeps only reads that overlap SNP positions in the VCF and have a valid
# cell barcode from the CellBender output.  This dramatically reduces BAM
# file size, making pileup generation tractable.

echo "=== Step A: Filtering BAM for ${LIBRARY_ID} ==="

source "${CONDA_INIT_SCRIPT}"
conda activate "${BCFTOOLS_ENV}"

bash "${SCRIPT_DIR}/popscle_helper_tools/filter_bam_file_for_popscle_dsc_pileup.sh" \
    "${DEJAGER_COUNTS}/${LIBRARY_ID}/outs/possorted_genome_bam.bam" \
    "${BARCODES}" \
    "${DEJAGER_DEMUX_VCF}" \
    "${OUT_DIR}/BAMOutput1.bam"

conda deactivate

echo "=== Step A complete ==="

# ===========================================================================
# Step B: Pileup generation
# ===========================================================================
# Generates allele count matrices at SNP positions for each cell barcode.

echo "=== Step B: Generating pileup for ${LIBRARY_ID} ==="

source /etc/profile.d/modules.sh
module load "${SINGULARITY_MODULE}"
unset SINGULARITY_VERIFY_CHECKS

singularity exec \
    --pwd "${DEJAGER_WGS_DIR}" \
    --bind "${DEJAGER_WGS_DIR}:/mnt" \
    "${DEMUXAFY_SIF}" \
    popscle_pileup.py \
        --sam "/mnt/${LIBRARY_ID}/BAMOutput1.bam" \
        --vcf "/mnt/${VCF_BASENAME}" \
        --group-list "/mnt/processed_feature_bc_matrix_cell_barcodes_${LIBRARY_ID}.csv" \
        --out "/mnt/${LIBRARY_ID}/plpDemux1" \
        --sm-list "/mnt/individ/individPat${LIBRARY_ID}.txt"

echo "=== Step B complete ==="

# ===========================================================================
# Step C: Demuxlet
# ===========================================================================
# Assigns each cell to its most likely patient of origin using genotype
# likelihoods from the WGS VCF.

echo "=== Step C: Running demuxlet for ${LIBRARY_ID} ==="

singularity exec \
    --pwd "${DEJAGER_WGS_DIR}" \
    --bind "${DEJAGER_WGS_DIR}:/mnt" \
    "${DEMUXAFY_SIF}" \
    popscle demuxlet \
        --plp "/mnt/${LIBRARY_ID}/plpDemux1" \
        --vcf "/mnt/${VCF_BASENAME}" \
        --field "PL" \
        --group-list "/mnt/processed_feature_bc_matrix_cell_barcodes_${LIBRARY_ID}.csv" \
        --sm-list "/mnt/individ/individPat${LIBRARY_ID}.txt" \
        --out "/mnt/${LIBRARY_ID}/demux1" \
        --alpha 0.05 \
        --min-mac 1 \
        --doublet-prior 0.1

echo "=== Step C complete: ${OUT_DIR}/demux1.best ==="
