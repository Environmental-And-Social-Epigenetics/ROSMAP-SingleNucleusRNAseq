"""
Generate and submit demuxlet batch scripts for each DeJager library.

The demuxlet workflow has two stages per library:
  A) BAM filtering  — reduce BAM to reads overlapping SNPs with valid barcodes
  B) Pileup + Demuxlet — generate pileup and run demuxlet for cell assignment

Reads paths from environment variables set by config/paths.sh.

Usage:
    source config/paths.sh
    python Demuxlet_DeJager.py --all            # generate scripts only
    python Demuxlet_DeJager.py --all --submit   # generate and submit
    python Demuxlet_DeJager.py --bam-only       # BAM filtering scripts only
    python Demuxlet_DeJager.py --demux-only     # pileup+demuxlet scripts only
"""

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path


# ---------------------------------------------------------------------------
# Resolve paths from environment (set by config/paths.sh)
# ---------------------------------------------------------------------------

# Compute repo/workspace roots from script location (depth 4 from repo root)
_script_path = Path(__file__).resolve()
_repo_root = _script_path.parents[3]
_dejager_root = _repo_root / "Data" / "Transcriptomics" / "DeJager"

DEJAGER_WGS_DIR = os.environ.get("DEJAGER_WGS_DIR", "__UNCONFIGURED__set_DEJAGER_WGS_DIR")
DEJAGER_COUNTS = os.environ.get(
    "DEJAGER_COUNTS",
    os.environ.get("DEJAGER_CELLRANGER", str(_dejager_root / "Cellranger_Output")),
)
DEJAGER_PREPROCESSED = os.environ.get(
    "DEJAGER_PREPROCESSED",
    os.environ.get("DEJAGER_CELLBENDER", str(_dejager_root / "Cellbender_Output")),
)
DEJAGER_DEMUX_VCF = os.environ.get(
    "DEJAGER_DEMUX_VCF",
    os.path.join(DEJAGER_WGS_DIR, "snp_fixedconcatenated_liftedROSMAP.vcf.gz"),
)
DEMUXAFY_SIF = os.environ.get(
    "DEMUXAFY_SIF", os.path.join(DEJAGER_WGS_DIR, "Demuxafy.sif")
)
DEJAGER_PATIENT_IDS_DIR = os.environ.get(
    "DEJAGER_PATIENT_IDS_DIR", os.path.join(DEJAGER_WGS_DIR, "individ")
)
CONDA_INIT_SCRIPT = os.environ.get(
    "CONDA_INIT_SCRIPT",
    os.path.join(os.environ.get("HOME", ""), "miniforge3/etc/profile.d/conda.sh"),
)
BCFTOOLS_ENV = os.environ.get(
    "BCFTOOLS_ENV",
    os.path.join(os.environ.get("HOME", ""), "conda_envs/bcftools_env"),
)
SINGULARITY_MODULE = os.environ.get("SINGULARITY_MODULE", "singularity/3.10.4")
SLURM_MAIL_USER = os.environ.get("SLURM_MAIL_USER", "")

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

# Output directories for generated scripts and logs
GENERATED_SCRIPTS_DIR = os.path.join(DEJAGER_WGS_DIR, "generated_scripts")
LOGS_DIR = os.path.join(DEJAGER_WGS_DIR, "logs")


# ---------------------------------------------------------------------------
# SLURM script templates
# ---------------------------------------------------------------------------

def _mail_line():
    if SLURM_MAIL_USER:
        return f"#SBATCH --mail-user={SLURM_MAIL_USER}"
    return ""


BAM_TEMPLATE = """#!/bin/bash
#SBATCH -n 45
#SBATCH -t 3:00:00
#SBATCH --mem=350G
#SBATCH --partition=pi_lhtsai
#SBATCH -o {logs_dir}/bam-%j.out
#SBATCH -e {logs_dir}/bam-%j.err
{mail_line}
#SBATCH --mail-type=FAIL

set -euo pipefail

source {conda_init}

cd {wgs_dir}

conda activate {bcftools_env}
cd {helper_tools_dir}

mkdir -p {wgs_dir}/{lib_id} || exit 1

./filter_bam_file_for_popscle_dsc_pileup.sh \\
    {counts_dir}/{lib_id}/outs/possorted_genome_bam.bam \\
    {preprocessed_dir}/{lib_id}/processed_feature_bc_matrix_cell_barcodes.csv \\
    {vcf} \\
    {wgs_dir}/{lib_id}/BAMOutput1.bam
"""

DEMUX_TEMPLATE = """#!/bin/bash
#SBATCH -n 10
#SBATCH -t 2-00:00:00
#SBATCH --mem=350G
#SBATCH --partition=pi_lhtsai
#SBATCH -o {logs_dir}/demux-%j.out
#SBATCH -e {logs_dir}/demux-%j.err
{mail_line}
#SBATCH --mail-type=FAIL

set -euo pipefail

source /etc/profile.d/modules.sh

cd {wgs_dir}

module load {singularity_module}

unset SINGULARITY_VERIFY_CHECKS

apptainer exec \\
    --pwd {wgs_dir} \\
    --bind {wgs_dir}:/mnt \\
    {sif} \\
    /opt/popscle/bin/popscle dsc-pileup \\
        --sam /mnt/{lib_id}/BAMOutput1.bam \\
        --vcf /mnt/{vcf_basename} \\
        --group-list /mnt/processed_feature_bc_matrix_cell_barcodes_{lib_id}.csv \\
        --out /mnt/{lib_id}/plpDemux1 \\
        --sm-list /mnt/individ/individPat{lib_id}.txt

apptainer exec \\
    --pwd {wgs_dir} \\
    --bind {wgs_dir}:/mnt \\
    {sif} \\
    /opt/popscle/bin/popscle demuxlet \\
        --plp /mnt/{lib_id}/plpDemux1 \\
        --vcf /mnt/{vcf_basename} \\
        --field "PL" \\
        --group-list /mnt/processed_feature_bc_matrix_cell_barcodes_{lib_id}.csv \\
        --sm-list /mnt/individ/individPat{lib_id}.txt \\
        --out /mnt/{lib_id}/demux1 \\
        --alpha 0.05 \\
        --min-mac 1 \\
        --doublet-prior 0.1
"""


# ---------------------------------------------------------------------------
# Library discovery
# ---------------------------------------------------------------------------

def discover_libraries(force=False):
    """Find libraries that have patient ID files and Cell Ranger BAMs."""
    libraries = []
    if not os.path.isdir(DEJAGER_PATIENT_IDS_DIR):
        print(f"ERROR: Patient IDs directory not found: {DEJAGER_PATIENT_IDS_DIR}")
        sys.exit(1)

    # Get libraries from individPat*.txt files
    for fname in sorted(os.listdir(DEJAGER_PATIENT_IDS_DIR)):
        if fname.startswith("individPat") and fname.endswith(".txt"):
            lib_id = fname[len("individPat"):-len(".txt")]
            if "alone" in lib_id:
                continue
            if not force:
                demux_best = os.path.join(DEJAGER_WGS_DIR, lib_id, "demux1.best")
                if os.path.exists(demux_best):
                    continue
            libraries.append(lib_id)
    return libraries


def require_configured(names):
    """Fail early when cluster-specific Demuxlet paths were not configured."""
    missing = []
    for name in names:
        value = globals()[name]
        if not value or "__UNCONFIGURED__" in str(value):
            missing.append(name)
    if missing:
        print("ERROR: DeJager demuxlet paths are not configured.")
        print("Source config/paths.sh and set these variables in config/paths.local.sh:")
        for name in missing:
            print(f"  - {name}")
        sys.exit(1)


# ---------------------------------------------------------------------------
# Script generation
# ---------------------------------------------------------------------------

def generate_bam_scripts(libraries):
    """Generate BAM filtering SLURM scripts."""
    os.makedirs(GENERATED_SCRIPTS_DIR, exist_ok=True)
    os.makedirs(LOGS_DIR, exist_ok=True)
    scripts = []
    for lib_id in libraries:
        content = BAM_TEMPLATE.format(
            logs_dir=LOGS_DIR,
            mail_line=_mail_line(),
            conda_init=CONDA_INIT_SCRIPT,
            wgs_dir=DEJAGER_WGS_DIR,
            bcftools_env=BCFTOOLS_ENV,
            helper_tools_dir=os.path.join(SCRIPT_DIR, "popscle_helper_tools"),
            counts_dir=DEJAGER_COUNTS,
            preprocessed_dir=DEJAGER_PREPROCESSED,
            vcf=DEJAGER_DEMUX_VCF,
            lib_id=lib_id,
        )
        filepath = os.path.join(GENERATED_SCRIPTS_DIR, f"scriptBAM_{lib_id}.sh")
        with open(filepath, "w") as f:
            f.write(content)
        scripts.append((lib_id, filepath))
    return scripts


def generate_demux_scripts(libraries):
    """Generate pileup + demuxlet SLURM scripts."""
    os.makedirs(GENERATED_SCRIPTS_DIR, exist_ok=True)
    os.makedirs(LOGS_DIR, exist_ok=True)
    vcf_basename = os.path.basename(DEJAGER_DEMUX_VCF)
    scripts = []
    for lib_id in libraries:
        content = DEMUX_TEMPLATE.format(
            logs_dir=LOGS_DIR,
            mail_line=_mail_line(),
            wgs_dir=DEJAGER_WGS_DIR,
            singularity_module=SINGULARITY_MODULE,
            sif=DEMUXAFY_SIF,
            vcf_basename=vcf_basename,
            lib_id=lib_id,
        )
        filepath = os.path.join(GENERATED_SCRIPTS_DIR, f"scriptDemux_{lib_id}.sh")
        with open(filepath, "w") as f:
            f.write(content)
        scripts.append((lib_id, filepath))
    return scripts


def submit_scripts(scripts, label=""):
    """Submit scripts to SLURM via sbatch."""
    for lib_id, filepath in scripts:
        try:
            result = subprocess.run(
                ["sbatch", filepath],
                capture_output=True, text=True, check=True,
            )
            print(f"  [{label}] {lib_id}: {result.stdout.strip()}")
        except subprocess.CalledProcessError as e:
            print(f"  [{label}] {lib_id}: FAILED - {e.stderr.strip()}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def _staged_barcode_path(lib_id):
    """Destination path for a library's barcode file inside DEJAGER_WGS_DIR."""
    return os.path.join(
        DEJAGER_WGS_DIR,
        f"processed_feature_bc_matrix_cell_barcodes_{lib_id}.csv",
    )


def stage_barcodes(libraries):
    """Copy each library's CellBender barcode CSV into DEJAGER_WGS_DIR with the
    library-suffixed name the demux scripts expect. Idempotent: skips files that
    are already staged. Replaces the manual copy loop previously in the README."""
    staged, missing = 0, []
    for lib_id in libraries:
        src = os.path.join(
            DEJAGER_PREPROCESSED, lib_id, "processed_feature_bc_matrix_cell_barcodes.csv"
        )
        dst = _staged_barcode_path(lib_id)
        if os.path.exists(dst):
            continue
        if not os.path.exists(src):
            missing.append((lib_id, src))
            continue
        shutil.copy2(src, dst)
        staged += 1
    print(f"  Staged {staged} barcode file(s) into {DEJAGER_WGS_DIR}")
    if missing:
        print(f"  WARNING: {len(missing)} libraries missing CellBender barcode files:")
        for lib_id, src in missing[:10]:
            print(f"    {lib_id}: {src}")
    return staged, missing


def check_barcodes_staged(libraries):
    """Preflight: confirm every library's barcode file is staged in DEJAGER_WGS_DIR.
    Returns the list of libraries whose barcode file is missing."""
    return [lib for lib in libraries if not os.path.exists(_staged_barcode_path(lib))]


def main():
    parser = argparse.ArgumentParser(
        description="Generate and optionally submit demuxlet batch scripts."
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--all", action="store_true",
                      help="Generate both BAM and demuxlet scripts")
    mode.add_argument("--bam-only", action="store_true",
                      help="Generate BAM filtering scripts only")
    mode.add_argument("--demux-only", action="store_true",
                      help="Generate pileup+demuxlet scripts only")
    mode.add_argument("--stage-barcodes", action="store_true",
                      help="Copy CellBender barcode CSVs into DEJAGER_WGS_DIR (run before --demux-only/--all)")

    parser.add_argument("--submit", action="store_true",
                        help="Submit generated scripts to SLURM")
    parser.add_argument("--force", action="store_true",
                        help="Regenerate scripts even if demux1.best already exists")
    args = parser.parse_args()

    required = ["DEJAGER_WGS_DIR", "DEJAGER_DEMUX_VCF", "DEMUXAFY_SIF", "DEJAGER_PATIENT_IDS_DIR"]
    if args.all or args.bam_only or args.stage_barcodes:
        required.extend(["DEJAGER_COUNTS", "DEJAGER_PREPROCESSED", "CONDA_INIT_SCRIPT", "BCFTOOLS_ENV"])
    require_configured(required)

    libraries = discover_libraries(force=args.force)
    if not libraries:
        print("No libraries to process (all have demux1.best or no patient IDs found).")
        print("Use --force to regenerate scripts for completed libraries.")
        return

    print(f"Found {len(libraries)} libraries to process.")

    if args.stage_barcodes:
        print("\nStaging CellBender barcode files into DEJAGER_WGS_DIR ...")
        stage_barcodes(libraries)
        return

    # Preflight: the pileup/demuxlet step reads the staged barcode files. Fail
    # loudly (with the exact fix) instead of letting popscle error on missing input.
    if args.all or args.demux_only:
        unstaged = check_barcodes_staged(libraries)
        if unstaged:
            print(f"\nERROR: {len(unstaged)} libraries are missing staged barcode files in "
                  f"{DEJAGER_WGS_DIR}.")
            print("       Run barcode staging first:")
            print("         python Demuxlet_DeJager.py --stage-barcodes")
            print(f"       Missing (first 10): {', '.join(unstaged[:10])}")
            sys.exit(1)

    if args.all or args.bam_only:
        print(f"\nGenerating BAM filtering scripts...")
        bam_scripts = generate_bam_scripts(libraries)
        print(f"  Generated {len(bam_scripts)} scripts in {GENERATED_SCRIPTS_DIR}")
        if args.submit:
            print("Submitting BAM filtering jobs...")
            submit_scripts(bam_scripts, label="BAM")

    if args.all or args.demux_only:
        print(f"\nGenerating pileup+demuxlet scripts...")
        demux_scripts = generate_demux_scripts(libraries)
        print(f"  Generated {len(demux_scripts)} scripts in {GENERATED_SCRIPTS_DIR}")
        if args.submit:
            print("Submitting pileup+demuxlet jobs...")
            submit_scripts(demux_scripts, label="Demux")

    if not args.submit:
        print("\nScripts generated but NOT submitted. Use --submit to submit to SLURM.")


if __name__ == "__main__":
    main()
