#!/usr/bin/env python3
"""
Generate Cell Ranger batch scripts for Tsai patients from local FASTQs.

Expected FASTQ layout:
    /om/scratch/Mon/mabdel03/Tsai_Data/FASTQs/<projid>/<Library_ID>/*.fastq.gz

Usage:
    python generate_batch_scripts.py [--dry-run]
"""

import argparse
import csv
import os
import shutil
from pathlib import Path


# Configuration - override via environment variables when needed
BATCH_SIZE = int(os.environ.get("TSAI_BATCH_SIZE", "30"))
REPO_ROOT = os.environ.get(
    "REPO_ROOT", "/om/scratch/Mon/mabdel03/ROSMAP-SingleNucleusRNAseq"
)
PIPELINE_DIR = os.environ.get(
    "PIPELINE_DIR", f"{REPO_ROOT}/Preprocessing/Tsai/02_Cellranger_Counts"
)

# Input FASTQs
FASTQS_ROOT = Path(
    os.environ.get("TSAI_FASTQS_DIR", "/om/scratch/Mon/mabdel03/Tsai_Data/FASTQs")
)

# Output directories
BATCH_SCRIPTS_DIR = f"{PIPELINE_DIR}/Batch_Scripts"
TRACKING_DIR = f"{PIPELINE_DIR}/Tracking"
LOGS_OUT = f"{PIPELINE_DIR}/Logs/Outs"
LOGS_ERR = f"{PIPELINE_DIR}/Logs/Errs"

# Scratch paths
CELLRANGER_OUTPUT = os.environ.get(
    "TSAI_CELLRANGER_OUTPUT", "/om/scratch/Mon/mabdel03/Tsai_Data/Cellranger_Outputs"
)
CELLBENDER_OUTPUT = os.environ.get(
    "TSAI_CELLBENDER_SCRATCH", "/om/scratch/Mon/mabdel03/Tsai/Cellbender_Output"
)

# Permanent storage
FINAL_OUTPUT = os.environ.get(
    "TSAI_PREPROCESSED",
    "/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/Data/Tsai/Preprocessed_Counts",
)

# Cell Ranger settings
CELLRANGER_PATH = os.environ.get(
    "CELLRANGER_PATH", "/om2/user/mabdel03/apps/yard/cellranger-8.0.0"
)
TRANSCRIPTOME = os.environ.get(
    "CELLRANGER_REF",
    "/om2/user/mabdel03/yard/references/human/refdata-gex-GRCh38-2020-A",
)

# Conda
CONDA_INIT = os.environ.get(
    "CONDA_INIT_SCRIPT",
    "/orcd/data/lhtsai/001/om2/mabdel03/miniforge3/etc/profile.d/conda.sh",
)
CELLBENDER_ENV = os.environ.get(
    "CELLBENDER_ENV", "/orcd/data/lhtsai/001/om2/mabdel03/conda_envs/Cellbender_env"
)


def generate_cellranger_script(projid, library_ids, fastq_dirs, batch_num):
    """Generate a Cell Ranger batch script for a single patient."""
    
    fastqs_arg = ",".join(fastq_dirs)
    output_dir = f"{CELLRANGER_OUTPUT}/{projid}"
    sample_arg = ",".join(library_ids)
    
    script = f"""#!/bin/bash
#SBATCH --job-name=cr_{projid}
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output={LOGS_OUT}/cellranger_{projid}_%j.out
#SBATCH --error={LOGS_ERR}/cellranger_{projid}_%j.err
#SBATCH --mail-user=mabdel03@mit.edu
#SBATCH --mail-type=FAIL

# Cell Ranger count for patient {projid} (Batch {batch_num})
# Library IDs: {",".join(library_ids)}
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
    --id={projid} \\
    --transcriptome={TRANSCRIPTOME} \\
    --sample={sample_arg} \\
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


def discover_patients(fastqs_root):
    """Discover patients and FASTQ directories from local layout."""
    patient_dirs = sorted([p for p in fastqs_root.iterdir() if p.is_dir()], key=lambda p: p.name)
    patient_data = []
    skipped = []

    for patient_dir in patient_dirs:
        library_dirs = sorted([d for d in patient_dir.iterdir() if d.is_dir()], key=lambda d: d.name)
        fastq_dirs = []
        library_ids = []

        for lib_dir in library_dirs:
            if any(lib_dir.glob("*.fastq.gz")):
                fastq_dirs.append(str(lib_dir))
                library_ids.append(lib_dir.name)
            else:
                run_dirs = [d for d in lib_dir.iterdir() if d.is_dir()]
                for run_dir in run_dirs:
                    if any(run_dir.glob("*.fastq.gz")):
                        fastq_dirs.append(str(run_dir))
                        library_ids.append(lib_dir.name)

        if not fastq_dirs:
            if any(patient_dir.glob("*.fastq.gz")):
                fastq_dirs.append(str(patient_dir))
                library_ids.append(patient_dir.name)

        if not fastq_dirs:
            skipped.append(patient_dir.name)
            continue

        n_fastqs = sum(len(list(Path(d).glob("*.fastq.gz"))) for d in fastq_dirs)

        patient_data.append(
            {
                "projid": patient_dir.name,
                "library_ids": library_ids,
                "fastq_dirs": fastq_dirs,
                "n_source_dirs": len(fastq_dirs),
                "n_fastqs": n_fastqs,
            }
        )

    return patient_data, skipped


def main():
    parser = argparse.ArgumentParser(description="Generate Cell Ranger batch scripts")
    parser.add_argument("--dry-run", action="store_true", help="Print summary without creating files")
    args = parser.parse_args()
    
    print("=" * 60)
    print("Generating Cell Ranger + CellBender Batch Scripts")
    print("=" * 60)
    
    # Discover patients from FASTQ directory layout
    print(f"\nScanning FASTQs: {FASTQS_ROOT}")
    if not FASTQS_ROOT.exists():
        raise SystemExit(f"ERROR: FASTQs root not found: {FASTQS_ROOT}")

    patient_data, skipped = discover_patients(FASTQS_ROOT)

    unique_library_ids = set()
    total_fastqs = 0
    for patient in patient_data:
        unique_library_ids.update(patient["library_ids"])
        total_fastqs += patient["n_fastqs"]

    print(f"Total FASTQ files: {total_fastqs}")
    print(f"Unique patients (projids): {len(patient_data)}")
    print(f"Unique Library_IDs: {len(unique_library_ids)}")

    if skipped:
        print(f"WARNING: Skipped {len(skipped)} patients with no FASTQs.")
    
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
    
    # Clean existing batch scripts to avoid stale jobs
    if os.path.isdir(BATCH_SCRIPTS_DIR):
        for entry in os.listdir(BATCH_SCRIPTS_DIR):
            if entry.startswith("batch_"):
                shutil.rmtree(os.path.join(BATCH_SCRIPTS_DIR, entry), ignore_errors=True)

    # Create batch subdirectories
    for batch_num in range(1, n_batches + 1):
        os.makedirs(f"{BATCH_SCRIPTS_DIR}/batch_{batch_num}/cellranger", exist_ok=True)
        os.makedirs(f"{BATCH_SCRIPTS_DIR}/batch_{batch_num}/cellbender", exist_ok=True)
    
    # Generate scripts for each patient
    print(f"\nGenerating scripts...")
    
    for patient in patient_data:
        projid = patient["projid"]
        library_ids = patient["library_ids"]
        source_dirs = patient["fastq_dirs"]
        batch_num = patient["batch"]
        
        # Cell Ranger script
        cr_script = generate_cellranger_script(projid, library_ids, source_dirs, batch_num)
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
    metadata_path = f"{TRACKING_DIR}/patient_metadata.csv"
    with open(metadata_path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "projid",
                "library_ids",
                "fastq_dirs",
                "n_source_dirs",
                "n_fastqs",
                "batch",
                "index",
            ],
        )
        writer.writeheader()
        for patient in patient_data:
            writer.writerow(
                {
                    "projid": patient["projid"],
                    "library_ids": "|".join(patient["library_ids"]),
                    "fastq_dirs": "|".join(patient["fastq_dirs"]),
                    "n_source_dirs": patient["n_source_dirs"],
                    "n_fastqs": patient["n_fastqs"],
                    "batch": patient["batch"],
                    "index": patient["index"],
                }
            )

    # Save batch assignments
    batch_path = f"{TRACKING_DIR}/batch_assignments.csv"
    with open(batch_path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["projid", "library_ids", "batch", "index"],
        )
        writer.writeheader()
        for patient in patient_data:
            writer.writerow(
                {
                    "projid": patient["projid"],
                    "library_ids": "|".join(patient["library_ids"]),
                    "batch": patient["batch"],
                    "index": patient["index"],
                }
            )
    
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

