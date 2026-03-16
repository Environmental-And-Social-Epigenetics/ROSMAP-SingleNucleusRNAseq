#!/bin/bash
#SBATCH --job-name=fix_20151388
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=fix_20151388_%j.out
#SBATCH --error=fix_20151388_%j.err
#SBATCH --mail-type=FAIL

# Fix chemistry conflict for patient 20151388 by running libraries separately
# - D19-2458: SC3Pv2 chemistry
# - D19-4144: SC3Pv3 chemistry

set -euo pipefail

# Source central path config
source "$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../.." && pwd)/config/paths.sh"

export PATH="${CELLRANGER_PATH}:$PATH"
export HDF5_USE_FILE_LOCKING=FALSE

OUTPUT_BASE="${TSAI_CELLRANGER_OUTPUT}"
FASTQ_BASE="${TSAI_FASTQS_DIR}"
TRANSCRIPTOME="${CELLRANGER_REF}"
TRACKING_DIR="${REPO_ROOT}/Preprocessing/Tsai/02_Cellranger_Counts/Tracking"

PATIENT="20151388"

echo "=========================================="
echo "Processing patient ${PATIENT}"
echo "=========================================="

# Clean previous failed attempt
rm -rf ${OUTPUT_BASE}/${PATIENT}
mkdir -p ${OUTPUT_BASE}/${PATIENT}

# Run D19-2458 with SC3Pv2 chemistry
echo "Running D19-2458 with SC3Pv2..."
cellranger count \
    --create-bam=false \
    --include-introns=true \
    --nosecondary \
    --r1-length=26 \
    --chemistry=SC3Pv2 \
    --id=${PATIENT}_D19-2458 \
    --transcriptome=${TRANSCRIPTOME} \
    --sample=D19-2458 \
    --fastqs=${FASTQ_BASE}/${PATIENT}/D19-2458 \
    --output-dir=${OUTPUT_BASE}/${PATIENT}/${PATIENT}_D19-2458

# Run D19-4144 with SC3Pv3 chemistry
echo "Running D19-4144 with SC3Pv3..."
cellranger count \
    --create-bam=false \
    --include-introns=true \
    --nosecondary \
    --r1-length=26 \
    --chemistry=SC3Pv3 \
    --id=${PATIENT}_D19-4144 \
    --transcriptome=${TRANSCRIPTOME} \
    --sample=D19-4144 \
    --fastqs=${FASTQ_BASE}/${PATIENT}/D19-4144/10x-3862H,${FASTQ_BASE}/${PATIENT}/D19-4144/HL32WDMXX \
    --output-dir=${OUTPUT_BASE}/${PATIENT}/${PATIENT}_D19-4144

# Create aggregation CSV
echo "Creating aggregation CSV..."
cat > ${OUTPUT_BASE}/${PATIENT}/aggr.csv << EOF
sample_id,molecule_h5
${PATIENT}_D19-2458,${OUTPUT_BASE}/${PATIENT}/${PATIENT}_D19-2458/outs/molecule_info.h5
${PATIENT}_D19-4144,${OUTPUT_BASE}/${PATIENT}/${PATIENT}_D19-4144/outs/molecule_info.h5
EOF

# Run aggregation
echo "Running cellranger aggr..."
cellranger aggr \
    --id=${PATIENT} \
    --csv=${OUTPUT_BASE}/${PATIENT}/aggr.csv \
    --normalize=none \
    --nosecondary \
    --output-dir=${OUTPUT_BASE}/${PATIENT}/aggregated

# Copy per-library raw matrix to expected CellBender input path
# Use D19-4144 (SC3Pv3 = newer chemistry, more cells)
mkdir -p ${OUTPUT_BASE}/${PATIENT}/outs
cp ${OUTPUT_BASE}/${PATIENT}/${PATIENT}_D19-4144/outs/raw_feature_bc_matrix.h5 \
   ${OUTPUT_BASE}/${PATIENT}/outs/raw_feature_bc_matrix.h5

# Also copy aggregated filtered matrix
cp ${OUTPUT_BASE}/${PATIENT}/aggregated/outs/count/filtered_feature_bc_matrix.h5 \
   ${OUTPUT_BASE}/${PATIENT}/filtered_feature_bc_matrix.h5 || true

echo "${PATIENT}" >> ${TRACKING_DIR}/cellranger_completed.txt
echo "Cell Ranger completed for ${PATIENT}"
