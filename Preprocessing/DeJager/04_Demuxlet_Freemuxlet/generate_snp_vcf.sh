#!/bin/bash
#SBATCH -J dej_snp_vcf
#SBATCH -t 6:00:00
#SBATCH -n 4
#SBATCH --mem=64G
#SBATCH --partition=mit_normal
#SBATCH --mail-type=FAIL

# Generate SNP-only VCF from the full WGS VCF for demuxlet.
# The full VCF contains all variant types; demuxlet needs only SNPs.
#
# Input:  fixedconcatenated_liftedROSMAP.vcf.gz  (~295 GB)
# Output: snp_fixedconcatenated_liftedROSMAP.vcf.gz (SNPs only)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
source "${REPO_ROOT}/config/paths.sh"

set +u
source "${CONDA_INIT_SCRIPT}"
conda activate "${BCFTOOLS_ENV}"
set -u

INPUT_VCF="${DEJAGER_FULL_WGS_VCF}"
OUTPUT_VCF="${DEJAGER_DEMUX_VCF}"

if [[ "${INPUT_VCF}" == *"__UNCONFIGURED__"* || "${OUTPUT_VCF}" == *"__UNCONFIGURED__"* ]]; then
    echo "ERROR: DeJager WGS paths are not configured. Set DEJAGER_WGS_DIR or DEJAGER_FULL_WGS_VCF/DEJAGER_DEMUX_VCF in config/paths.local.sh."
    exit 1
fi

if [[ ! -f "${INPUT_VCF}" ]]; then
    echo "ERROR: Input VCF not found: ${INPUT_VCF}"
    exit 1
fi

echo "Filtering SNPs from ${INPUT_VCF}..."
bcftools view -v snps "${INPUT_VCF}" -Oz -o "${OUTPUT_VCF}" --threads 4

echo "Indexing output VCF..."
tabix -p vcf "${OUTPUT_VCF}"

echo "Done. Output: ${OUTPUT_VCF}"
ls -lh "${OUTPUT_VCF}" "${OUTPUT_VCF}.tbi"
