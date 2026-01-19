#!/usr/bin/env python3
"""
Generate Cell Ranger batch scripts for all Tsai patients.

This script reads All_ROSMAP_FASTQs.csv and generates:
1. Per-patient Cell Ranger batch scripts
2. Per-patient CellBender batch scripts
3. Batch assignment file (which patients are in which batch)
4. Master patient list with all metadata

Usage:
    python generate_batch_scripts.py [--dry-run]
"""

import argparse
import os
import pandas as pd
from pathlib import Path


# Configuration - these match cellranger_config.sh
BATCH_SIZE = 30
REPO_ROOT = "/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/ROSMAP-SingleNucleusRNAseq"
PIPELINE_DIR = f"{REPO_ROOT}/Preprocessing/Tsai/02_Cellranger_Counts"
INPUT_CSV = f"{REPO_ROOT}/Data/Tsai/All_ROSMAP_FASTQs.csv"

# Output directories
BATCH_SCRIPTS_DIR = f"{PIPELINE_DIR}/Batch_Scripts"
TRACKING_DIR = f"{PIPELINE_DIR}/Tracking"
LOGS_OUT = f"{PIPELINE_DIR}/Logs/Outs"
LOGS_ERR = f"{PIPELINE_DIR}/Logs/Errs"

# Scratch paths
SCRATCH_ROOT = "/home/mabdel03/orcd/scratch"
CELLRANGER_OUTPUT = f"{SCRATCH_ROOT}/Tsai/Cellranger_Counts"
CELLBENDER_OUTPUT = f"{SCRATCH_ROOT}/Tsai/Cellbender_Output"

# Permanent storage
FINAL_OUTPUT = "/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/Tsai/Preprocessed_Counts"

# Cell Ranger settings
CELLRANGER_PATH = "/orcd/data/lhtsai/001/om2/mabdel03/apps/yard/cellranger-8.0.0"
TRANSCRIPTOME = "/orcd/data/lhtsai/001/om2/mabdel03/yard/references/human/refdata-gex-GRCh38-2020-A"

# Conda
CONDA_INIT = "/orcd/data/lhtsai/001/om2/mabdel03/miniforge3/etc/profile.d/conda.sh"
CELLBENDER_ENV = "/orcd/data/lhtsai/001/om2/mabdel03/conda_envs/Cellbender_env"


def generate_cellranger_script(projid, library_id, fastq_dirs, batch_num):
    """Generate a Cell Ranger batch script for a single patient."""
    
    fastqs_arg = ",".join(fastq_dirs)
    output_dir = f"{CELLRANGER_OUTPUT}/{projid}"
    
    script = f"""#!/bin/bash
#SBATCH --job-name=cr_{projid}
#SBATCH --partition=mit_preemptable
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output={LOGS_OUT}/cellranger_{projid}_%j.out
#SBATCH --error={LOGS_ERR}/cellranger_{projid}_%j.err
#SBATCH --mail-user=mabdel03@mit.edu
#SBATCH --mail-type=FAIL

# Cell Ranger count for patient {projid} (Batch {batch_num})
# Library ID: {library_id}
# FASTQ directories: {len(fastq_dirs)}

set -e

# Add Cell Ranger to PATH
export PATH={CELLRANGER_PATH}:$PATH

# Create output directory
mkdir -p {output_dir}

# Run Cell Ranger count
cellranger count \\
    --create-bam=false \\
    --include-introns=true \\
    --nosecondary \\
    --r1-length=26 \\
    --id={library_id} \\
    --transcriptome={TRANSCRIPTOME} \\
    --sample={library_id} \\
    --fastqs={fastqs_arg} \\
    --output-dir={output_dir}

# Mark completion
echo "{projid}" >> {TRACKING_DIR}/cellranger_completed.txt

echo "Cell Ranger completed for {projid}"
"""
    return script


def generate_cellbender_script(projid, batch_num):
    """Generate a CellBender batch script for a single patient."""
    
    input_h5 = f"{CELLRANGER_OUTPUT}/{projid}/outs/raw_feature_bc_matrix.h5"
    output_dir = f"{CELLBENDER_OUTPUT}/{projid}"
    output_h5 = f"{output_dir}/cellbender_output.h5"
    final_dir = f"{FINAL_OUTPUT}/{projid}"
    
    script = f"""#!/bin/bash
#SBATCH --job-name=cb_{projid}
#SBATCH --partition=mit_normal_gpu
#SBATCH --gres=gpu:1
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --output={LOGS_OUT}/cellbender_{projid}_%j.out
#SBATCH --error={LOGS_ERR}/cellbender_{projid}_%j.err
#SBATCH --mail-user=mabdel03@mit.edu
#SBATCH --mail-type=FAIL

# CellBender for patient {projid} (Batch {batch_num})

set -e

# Initialize conda
source {CONDA_INIT}
conda activate {CELLBENDER_ENV}

# Create output directory
mkdir -p {output_dir}

# Check input exists
if [[ ! -f "{input_h5}" ]]; then
    echo "ERROR: Input file not found: {input_h5}"
    echo "{projid}" >> {TRACKING_DIR}/cellbender_failed.txt
    exit 1
fi

# Run CellBender
cellbender remove-background \\
    --input {input_h5} \\
    --output {output_h5} \\
    --expected-cells 5000 \\
    --total-droplets-included 20000 \\
    --fpr 0.01 \\
    --epochs 150 \\
    --cuda

# Copy to permanent storage
mkdir -p {final_dir}
cp {output_h5} {final_dir}/
cp {output_dir}/cellbender_output_filtered.h5 {final_dir}/ 2>/dev/null || true

# Mark completion
echo "{projid}" >> {TRACKING_DIR}/cellbender_completed.txt

echo "CellBender completed for {projid}"
"""
    return script


def main():
    parser = argparse.ArgumentParser(description="Generate Cell Ranger batch scripts")
    parser.add_argument("--dry-run", action="store_true", help="Print summary without creating files")
    args = parser.parse_args()
    
    print("=" * 60)
    print("Generating Cell Ranger + CellBender Batch Scripts")
    print("=" * 60)
    
    # Read the FASTQ CSV
    print(f"\nReading: {INPUT_CSV}")
    df = pd.read_csv(INPUT_CSV)
    
    print(f"Total FASTQ records: {len(df)}")
    print(f"Unique patients (projids): {df['projid'].nunique()}")
    print(f"Unique Library_IDs: {df['Library_ID'].nunique()}")
    
    # Group by projid to get patient-level information
    patient_data = []
    
    for projid in df['projid'].unique():
        patient_df = df[df['projid'] == projid]
        
        # Get unique source directories
        source_dirs = patient_df['source_dir'].unique().tolist()
        
        # Get the primary Library_ID (first one)
        library_id = patient_df['Library_ID'].iloc[0]
        
        # Count FASTQ files
        n_fastqs = len(patient_df)
        
        patient_data.append({
            'projid': projid,
            'library_id': library_id,
            'source_dirs': source_dirs,
            'n_source_dirs': len(source_dirs),
            'n_fastqs': n_fastqs
        })
    
    # Sort by projid for consistent ordering
    patient_data = sorted(patient_data, key=lambda x: x['projid'])
    
    # Assign batch numbers
    for i, patient in enumerate(patient_data):
        patient['batch'] = (i // BATCH_SIZE) + 1
        patient['index'] = i + 1
    
    # Summary
    n_patients = len(patient_data)
    n_batches = (n_patients + BATCH_SIZE - 1) // BATCH_SIZE
    
    print(f"\nPatient Summary:")
    print(f"  Total patients: {n_patients}")
    print(f"  Batch size: {BATCH_SIZE}")
    print(f"  Number of batches: {n_batches}")
    
    for batch_num in range(1, n_batches + 1):
        batch_patients = [p for p in patient_data if p['batch'] == batch_num]
        print(f"  Batch {batch_num}: {len(batch_patients)} patients")
    
    if args.dry_run:
        print("\n[DRY RUN] Skipping file creation")
        return
    
    # Create directories
    os.makedirs(BATCH_SCRIPTS_DIR, exist_ok=True)
    os.makedirs(TRACKING_DIR, exist_ok=True)
    os.makedirs(LOGS_OUT, exist_ok=True)
    os.makedirs(LOGS_ERR, exist_ok=True)
    
    # Create batch subdirectories
    for batch_num in range(1, n_batches + 1):
        os.makedirs(f"{BATCH_SCRIPTS_DIR}/batch_{batch_num}/cellranger", exist_ok=True)
        os.makedirs(f"{BATCH_SCRIPTS_DIR}/batch_{batch_num}/cellbender", exist_ok=True)
    
    # Generate scripts for each patient
    print(f"\nGenerating scripts...")
    
    for patient in patient_data:
        projid = patient['projid']
        library_id = patient['library_id']
        source_dirs = patient['source_dirs']
        batch_num = patient['batch']
        
        # Cell Ranger script
        cr_script = generate_cellranger_script(projid, library_id, source_dirs, batch_num)
        cr_path = f"{BATCH_SCRIPTS_DIR}/batch_{batch_num}/cellranger/{projid}_cellranger.sh"
        with open(cr_path, 'w') as f:
            f.write(cr_script)
        os.chmod(cr_path, 0o755)
        
        # CellBender script
        cb_script = generate_cellbender_script(projid, batch_num)
        cb_path = f"{BATCH_SCRIPTS_DIR}/batch_{batch_num}/cellbender/{projid}_cellbender.sh"
        with open(cb_path, 'w') as f:
            f.write(cb_script)
        os.chmod(cb_path, 0o755)
    
    # Save patient metadata
    metadata_df = pd.DataFrame(patient_data)
    metadata_df['source_dirs'] = metadata_df['source_dirs'].apply(lambda x: '|'.join(x))
    metadata_df.to_csv(f"{TRACKING_DIR}/patient_metadata.csv", index=False)
    
    # Save batch assignments
    batch_df = metadata_df[['projid', 'library_id', 'batch', 'index']]
    batch_df.to_csv(f"{TRACKING_DIR}/batch_assignments.csv", index=False)
    
    # Initialize tracking files
    open(f"{TRACKING_DIR}/cellranger_completed.txt", 'a').close()
    open(f"{TRACKING_DIR}/cellbender_completed.txt", 'a').close()
    open(f"{TRACKING_DIR}/cellranger_failed.txt", 'a').close()
    open(f"{TRACKING_DIR}/cellbender_failed.txt", 'a').close()
    
    print(f"\nGenerated scripts:")
    print(f"  Cell Ranger scripts: {n_patients}")
    print(f"  CellBender scripts: {n_patients}")
    print(f"  Batch assignments: {TRACKING_DIR}/batch_assignments.csv")
    print(f"  Patient metadata: {TRACKING_DIR}/patient_metadata.csv")
    
    print("\n" + "=" * 60)
    print("Script generation complete!")
    print("=" * 60)
    print(f"\nNext steps:")
    print(f"  1. Review generated scripts in {BATCH_SCRIPTS_DIR}/")
    print(f"  2. Run batch 1: ./Scripts/run_batch.sh 1")
    print(f"  3. After completion, run remaining batches 2-{n_batches}")


if __name__ == "__main__":
    main()

