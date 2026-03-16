"""
Generate and submit Cell Ranger count batch scripts for each DeJager library.

Reads paths from environment variables set by config/paths.sh.

Usage:
    source config/paths.sh
    python Count_DeJager.py
"""

import os
import subprocess
from pathlib import Path

import pandas as pd


# Compute repo/workspace roots from script location (depth 3 from repo root)
_script_path = Path(__file__).resolve()
_repo_root = _script_path.parents[3]
_workspace_root = _repo_root.parent

# Resolve paths from environment (set by config/paths.sh)
DATA_ROOT = os.environ.get(
    "DATA_ROOT",
    "__UNCONFIGURED__set_DATA_ROOT_in_paths_local_sh",
)
DEJAGER_FASTQS = os.environ.get(
    "DEJAGER_FASTQS",
    str(_workspace_root / "DeJager_Data" / "FASTQs"),
)
DEJAGER_COUNTS = os.environ.get(
    "DEJAGER_COUNTS",
    str(_workspace_root / "DeJager_Data" / "Counts"),
)
CELLRANGER_PATH = os.environ.get(
    "CELLRANGER_PATH",
    "__UNCONFIGURED__set_CELLRANGER_PATH_in_paths_local_sh",
)
CELLRANGER_REF = os.environ.get(
    "CELLRANGER_REF",
    "__UNCONFIGURED__set_CELLRANGER_REF_in_paths_local_sh",
)
SLURM_MAIL_USER = os.environ.get("SLURM_MAIL_USER", "")

SYNAPSE_CSV = os.path.join(
    DATA_ROOT, "Data/DeJager/FASTQs_Download/FASTQ_Download_CSVs/Synapse_FASTQ_IDs.csv"
)

# Output directories
BATCH_SCRIPTS_DIR = os.path.join(
    DATA_ROOT, "Data/DeJager/Original_Counts/Counts_Batch_Scripts"
)
LOGS_OUT = os.path.join(DATA_ROOT, "Data/DeJager/Original_Counts/Counts_outs")
LOGS_ERR = os.path.join(DATA_ROOT, "Data/DeJager/Original_Counts/Counts_errors")

df = pd.read_csv(SYNAPSE_CSV)
in_root = DEJAGER_FASTQS
out_root = DEJAGER_COUNTS

# Create output directories
os.makedirs(BATCH_SCRIPTS_DIR, exist_ok=True)
os.makedirs(LOGS_OUT, exist_ok=True)
os.makedirs(LOGS_ERR, exist_ok=True)

# Make folders for each library's counts outputs
for libID in set(df["LibraryID"]):
    dir_path = os.path.join(out_root, libID)
    os.makedirs(dir_path, exist_ok=True)

# Build mail line for SBATCH header
mail_line = ""
if SLURM_MAIL_USER:
    mail_line = f"#SBATCH --mail-user={SLURM_MAIL_USER}"

# Iterate over each library and generate + submit batch scripts
for libID in set(df["LibraryID"]):
    fastqs = os.path.join(in_root, libID)
    files = os.listdir(fastqs)
    samples = [piece.split("_")[0] + "_" + piece.split("_")[1] for piece in files]
    unique_samples = list(set(samples))

    if len(unique_samples) == 2:
        sample = f"{unique_samples[0]},{unique_samples[1]}"
    else:
        sample = unique_samples[0]

    filename = libID + "_count.sh"
    filepath = os.path.join(BATCH_SCRIPTS_DIR, filename)
    out_dir = os.path.join(out_root, libID)

    output = f"""#!/bin/bash
#SBATCH -t 47:00:00
#SBATCH -n 32
#SBATCH --mem=128G
#SBATCH --output={LOGS_OUT}/slurm-%j.out
#SBATCH --error={LOGS_ERR}/slurm-%j.err
{mail_line}
#SBATCH --mail-type=FAIL
export PATH={CELLRANGER_PATH}:$PATH
cellranger count --create-bam true --include-introns true --nosecondary --r1-length 26 --id {libID} --transcriptome={CELLRANGER_REF} --sample {sample} --fastqs {fastqs} --output-dir={out_dir}
"""
    with open(filepath, "w") as f:
        f.write(output)

    sbatch_command = f"sbatch {filepath}"
    process = subprocess.Popen(sbatch_command.split(), stdout=subprocess.PIPE)
    stdout, _ = process.communicate()
    print(f"Submitted {libID}: {stdout.decode().strip()}")
