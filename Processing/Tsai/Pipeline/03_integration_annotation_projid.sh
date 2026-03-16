#!/bin/bash
#SBATCH -J tsai_integrate_projid
#SBATCH -t 48:00:00
#SBATCH -n 32
#SBATCH --mem=500G
#SBATCH --mail-type=BEGIN,END,FAIL

set -euo pipefail

# Hardcoded for standalone sbatch submission (BASH_SOURCE resolves to
# SLURM's temp copy, not the original path).
REPO_ROOT="/om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq"
SCRIPT_DIR="${REPO_ROOT}/Processing/Tsai/Pipeline"

source "${REPO_ROOT}/config/paths.sh"

# NOTE: SLURM log redirection is handled by submit_pipeline.sh via
# command-line --output/--error flags.
# Stage 3 is a single integration job (not a SLURM array) because it
# must load and concatenate all samples simultaneously.
mkdir -p "${TSAI_PROCESSING_LOGS}"

export HDF5_USE_FILE_LOCKING=FALSE

# Temporarily relax nounset for conda activation (some activate.d scripts
# reference unset variables like ADDR2LINE).
set +u
source "${CONDA_INIT_SCRIPT}"
conda activate "${BATCHCORR_ENV}"
set -u
if [[ -z "${CONDA_PREFIX:-}" ]]; then
    echo "ERROR: Failed to activate conda environment: ${BATCHCORR_ENV}"
    exit 1
fi

INPUT_DIR="${TSAI_DOUBLET_REMOVED}"
OUTPUT_DIR="${TSAI_PROCESSING_OUTPUTS}/03_Integrated_projid"
METADATA_CSV="${TSAI_METADATA_CSV}"
MARKERS_RDS="${TSAI_MARKERS_RDS}"

python "${SCRIPT_DIR}/03_integration_annotation.py" \
    --input-dir "${INPUT_DIR}" \
    --output-dir "${OUTPUT_DIR}" \
    --metadata-csv "${METADATA_CSV}" \
    --markers-rds "${MARKERS_RDS}" \
    --harmony-batch-key projid
