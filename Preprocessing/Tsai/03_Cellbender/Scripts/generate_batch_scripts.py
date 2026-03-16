#!/usr/bin/env python3
"""
Generate CellBender batch scripts for Tsai patients from CellRanger outputs.

Uses batch assignments from:
  Preprocessing/Tsai/02_Cellranger_Counts/Tracking/batch_assignments.csv
"""

import argparse
import csv
import os
import shutil
from pathlib import Path


def env_bool(value, default=False):
    if value is None:
        return default
    return value.strip().lower() in {"1", "true", "yes", "y"}


# Compute repo/workspace roots from script location (depth 4 from repo root)
_script_path = Path(__file__).resolve()
_repo_root = _script_path.parents[3]
_workspace_root = _repo_root.parent

# Defaults (override with environment variables)
REPO_ROOT = os.environ.get(
    "REPO_ROOT", str(_repo_root)
)
PIPELINE_DIR = os.environ.get(
    "PIPELINE_DIR", f"{REPO_ROOT}/Preprocessing/Tsai/03_Cellbender"
)

BATCH_SCRIPTS_DIR = os.environ.get(
    "BATCH_SCRIPTS_DIR", f"{PIPELINE_DIR}/Batch_Scripts"
)
TRACKING_DIR = os.environ.get("TRACKING_DIR", f"{PIPELINE_DIR}/Tracking")
LOGS_OUT = os.environ.get("LOGS_OUT", f"{PIPELINE_DIR}/Logs/Outs")
LOGS_ERR = os.environ.get("LOGS_ERR", f"{PIPELINE_DIR}/Logs/Errs")

CELLRANGER_OUTPUT = os.environ.get(
    "CELLRANGER_OUTPUT",
    os.environ.get("TSAI_CELLRANGER_OUTPUT", str(_workspace_root / "Tsai_Data" / "Cellranger_Outputs")),
)
CELLBENDER_OUTPUT = os.environ.get(
    "CELLBENDER_OUTPUT",
    os.environ.get("TSAI_CELLBENDER_SCRATCH", str(_workspace_root / "Tsai" / "Cellbender_Output")),
)
FINAL_OUTPUT = os.environ.get(
    "FINAL_OUTPUT",
    os.environ.get(
        "TSAI_PREPROCESSED",
        str(_workspace_root / "Tsai_Data" / "Cellbender_Outputs"),
    ),
)

BATCH_ASSIGNMENTS_CSV = os.environ.get(
    "BATCH_ASSIGNMENTS_CSV",
    f"{REPO_ROOT}/Preprocessing/Tsai/02_Cellranger_Counts/Tracking/batch_assignments.csv",
)
CELLRANGER_COMPLETED = os.environ.get(
    "CELLRANGER_COMPLETED",
    f"{REPO_ROOT}/Preprocessing/Tsai/02_Cellranger_Counts/Tracking/cellranger_completed.txt",
)

# CellBender parameters
CB_FPR = os.environ.get("CB_FPR", "0")
CB_CUDA = env_bool(os.environ.get("CB_CUDA", "true"), default=True)

# SLURM settings
CB_SLURM_TIME = os.environ.get("CB_SLURM_TIME", "47:00:00")
CB_SLURM_CPUS = os.environ.get("CB_SLURM_CPUS", "32")
CB_SLURM_MEM = os.environ.get("CB_SLURM_MEM", "128G")
CB_SLURM_GPU = os.environ.get("CB_SLURM_GPU", "1")
SLURM_MAIL_USER = os.environ.get("SLURM_MAIL_USER", "")

CONDA_INIT_SCRIPT = os.environ.get(
    "CONDA_INIT_SCRIPT",
    os.path.join(os.environ.get("HOME", ""), "miniforge3/etc/profile.d/conda.sh"),
)
CELLBENDER_CONDA = os.environ.get(
    "CELLBENDER_CONDA",
    os.environ.get("CELLBENDER_ENV", os.path.join(os.environ.get("HOME", ""), "conda_envs/Cellbender_env")),
)


def load_batch_assignments(csv_path):
    assignments = []
    with open(csv_path, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            assignments.append(
                {
                    "projid": row["projid"],
                    "library_ids": row.get("library_ids", ""),
                    "batch": int(row["batch"]),
                    "index": int(row["index"]),
                }
            )
    return assignments


def load_completed(path):
    if not os.path.exists(path):
        return set()
    with open(path) as handle:
        return {line.strip() for line in handle if line.strip()}


def generate_cellbender_script(projid, batch_num):
    input_h5 = f"{CELLRANGER_OUTPUT}/{projid}/outs/raw_feature_bc_matrix.h5"
    output_dir = f"{FINAL_OUTPUT}/{projid}"
    output_h5 = f"{output_dir}/processed_feature_bc_matrix.h5"

    cuda_flag = "--cuda " if CB_CUDA else ""

    script = f"""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={CB_SLURM_CPUS}
#SBATCH --time={CB_SLURM_TIME}
#SBATCH --output={LOGS_OUT}/slurm-%j.out
#SBATCH --error={LOGS_ERR}/slurm-%j.err
#SBATCH --gres=gpu:{CB_SLURM_GPU}
#SBATCH --mem={CB_SLURM_MEM}
#SBATCH --mail-user={SLURM_MAIL_USER}
#SBATCH --mail-type=FAIL

# Disable HDF5 file locking to avoid shared filesystem issues
export HDF5_USE_FILE_LOCKING=FALSE

source {CONDA_INIT_SCRIPT}

conda activate {CELLBENDER_CONDA}

mkdir -p {output_dir}
cd {output_dir}

cellbender remove-background {cuda_flag}--input {input_h5} --fpr {CB_FPR} --output {output_h5}

# Mark completion
echo "{projid}" >> {TRACKING_DIR}/cellbender_completed.txt
"""
    return script


def main():
    parser = argparse.ArgumentParser(description="Generate CellBender batch scripts")
    parser.add_argument("--dry-run", action="store_true", help="Print summary without creating files")
    parser.add_argument("--batch", type=int, help="Generate scripts for a specific batch only")
    parser.add_argument(
        "--only-completed",
        action="store_true",
        help="Only generate scripts for samples with completed CellRanger runs",
    )
    args = parser.parse_args()

    print("=" * 60)
    print("Generating CellBender Batch Scripts")
    print("=" * 60)

    if not os.path.exists(BATCH_ASSIGNMENTS_CSV):
        raise SystemExit(f"ERROR: Batch assignments not found: {BATCH_ASSIGNMENTS_CSV}")

    assignments = load_batch_assignments(BATCH_ASSIGNMENTS_CSV)
    if args.batch:
        assignments = [a for a in assignments if a["batch"] == args.batch]

    if args.only_completed:
        completed = load_completed(CELLRANGER_COMPLETED)
        assignments = [a for a in assignments if a["projid"] in completed]

    if not assignments:
        raise SystemExit("No matching patients found to generate scripts.")

    batches = sorted({a["batch"] for a in assignments})

    print(f"Total patients: {len(assignments)}")
    print(f"Batches: {', '.join(str(b) for b in batches)}")
    if args.only_completed:
        print(f"Filtering: only completed CellRanger samples ({CELLRANGER_COMPLETED})")

    if args.dry_run:
        print("[DRY RUN] Skipping file creation")
        return

    os.makedirs(BATCH_SCRIPTS_DIR, exist_ok=True)
    os.makedirs(TRACKING_DIR, exist_ok=True)
    os.makedirs(LOGS_OUT, exist_ok=True)
    os.makedirs(LOGS_ERR, exist_ok=True)

    # Clean existing batch script directories to avoid stale jobs
    for entry in os.listdir(BATCH_SCRIPTS_DIR):
        if entry.startswith("batch_"):
            shutil.rmtree(os.path.join(BATCH_SCRIPTS_DIR, entry), ignore_errors=True)

    # Create batch subdirectories
    for batch_num in batches:
        os.makedirs(f"{BATCH_SCRIPTS_DIR}/batch_{batch_num}/cellbender", exist_ok=True)

    # Generate scripts
    for patient in assignments:
        projid = patient["projid"]
        batch_num = patient["batch"]
        cb_script = generate_cellbender_script(projid, batch_num)
        cb_path = f"{BATCH_SCRIPTS_DIR}/batch_{batch_num}/cellbender/{projid}_cellbender.sh"
        with open(cb_path, "w") as handle:
            handle.write(cb_script)
        os.chmod(cb_path, 0o755)

    # Initialize tracking files
    open(f"{TRACKING_DIR}/cellbender_completed.txt", "a").close()
    open(f"{TRACKING_DIR}/cellbender_failed.txt", "a").close()

    print(f"Generated CellBender scripts: {len(assignments)}")
    print(f"Scripts directory: {BATCH_SCRIPTS_DIR}")
    print("=" * 60)


if __name__ == "__main__":
    main()
