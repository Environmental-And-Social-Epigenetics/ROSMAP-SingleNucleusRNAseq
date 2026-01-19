#!/usr/bin/env python3
"""
Generate Globus batch transfer file for FASTQ files.

This script reads All_ROSMAP_FASTQs.csv and generates a Globus batch file
that maps each FASTQ from its source location on Engaging to the destination
on Openmind, organized by projid/Library_ID/run_id.

Note: Some samples have FASTQs with identical filenames from different 
sequencing runs (flow cells). This script creates unique subdirectories
per source directory to preserve all files.

Usage:
    python generate_globus_batch.py [--dry-run]

Output:
    Batch_Files/globus_batch.txt - Globus batch transfer file
"""

import argparse
import csv
import os
import re
import hashlib
from pathlib import Path
from datetime import datetime
from collections import defaultdict


# =============================================================================
# Configuration
# =============================================================================

# Input CSV with FASTQ locations
REPO_ROOT = "/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/ROSMAP-SingleNucleusRNAseq"
INPUT_CSV = f"{REPO_ROOT}/Data/Tsai/All_ROSMAP_FASTQs.csv"

# Script directory (for relative paths to output)
SCRIPT_DIR = Path(__file__).parent.parent
BATCH_FILE = SCRIPT_DIR / "Batch_Files" / "globus_batch.txt"
LOG_DIR = SCRIPT_DIR / "Logs"

# Destination base path on Openmind
DEST_BASE = "/om/scratch/Mon/mabdel03/Tsai_Data/FASTQs"

# Globus endpoint IDs
# Note: Use MIT ORCD Engaging Collection, not the older mithpc#engaging
ENGAGING_ENDPOINT = "ec54b570-cac5-47f7-b2a1-100c2078686f"  # MIT ORCD Engaging Collection
OPENMIND_ENDPOINT = "cbc6f8da-d37e-11eb-bde9-5111456017d9"


def extract_run_id(source_dir):
    """
    Extract a unique run identifier from the source directory path.
    
    Examples:
        /nfs/.../190930Tsa/10x-4182G/mkfastq/... -> "10x-4182G"
        /nfs/.../220311Tsa -> "220311Tsa"
        /nfs/.../SM-GZQS1/HKMG5DMXX/D19-5948 -> "HKMG5DMXX"
    """
    # Try to extract flow cell ID (common patterns)
    # Pattern 1: 10x-XXXXL or 10x-XXXXF or 10x-XXXXG or 10x-XXXXH
    match = re.search(r'(10x-\d+[A-Z])', source_dir)
    if match:
        return match.group(1)
    
    # Pattern 2: Flow cell ID like HKMG5DMXX, HL2FCDMXX
    match = re.search(r'/([A-Z0-9]{9,10})/', source_dir)
    if match:
        return match.group(1)
    
    # Pattern 3: Date-based directory like 220311Tsa
    match = re.search(r'/(\d{6}Tsa[A-Za-z]?)(?:/|$)', source_dir)
    if match:
        return match.group(1)
    
    # Fallback: use hash of source_dir for uniqueness
    return hashlib.md5(source_dir.encode()).hexdigest()[:8]


def generate_batch_file(dry_run=False):
    """Generate the Globus batch transfer file."""
    
    print("=" * 70)
    print("Globus Batch File Generator for Tsai FASTQ Transfer")
    print("=" * 70)
    print()
    
    # Verify input CSV exists
    if not os.path.exists(INPUT_CSV):
        print(f"ERROR: Input CSV not found: {INPUT_CSV}")
        return False
    
    print(f"Input CSV: {INPUT_CSV}")
    print(f"Destination base: {DEST_BASE}")
    print(f"Output batch file: {BATCH_FILE}")
    print()
    
    # First pass: identify source directories per projid/library_id
    # to determine if we need run subdirectories
    projid_lib_sources = defaultdict(set)
    
    with open(INPUT_CSV, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            key = (row['projid'], row['Library_ID'])
            projid_lib_sources[key].add(row['source_dir'])
    
    # Identify which projid/library combinations have multiple sources
    multi_source = {k: v for k, v in projid_lib_sources.items() if len(v) > 1}
    print(f"Patients with multiple sequencing runs: {len(multi_source)}")
    
    # Read CSV and generate batch lines
    batch_lines = []
    projid_set = set()
    library_id_set = set()
    run_id_set = set()
    total_size = 0
    
    # Track source_dir -> run_id mapping for consistency
    source_to_run = {}
    
    with open(INPUT_CSV, 'r') as f:
        reader = csv.DictReader(f)
        
        for row in reader:
            projid = row['projid']
            library_id = row['Library_ID']
            fastq_filename = row['fastq_filename']
            source_path = row['full_path']
            source_dir = row['source_dir']
            file_size = int(row['file_size'])
            
            key = (projid, library_id)
            
            # Determine if we need a run subdirectory
            if key in multi_source:
                # Multiple sources - add run_id subdirectory
                if source_dir not in source_to_run:
                    run_id = extract_run_id(source_dir)
                    # Handle potential run_id collisions within same projid/lib
                    base_run_id = run_id
                    counter = 1
                    while any(source_to_run.get(s) == run_id 
                             for s in projid_lib_sources[key] 
                             if s in source_to_run and s != source_dir):
                        run_id = f"{base_run_id}_{counter}"
                        counter += 1
                    source_to_run[source_dir] = run_id
                else:
                    run_id = source_to_run[source_dir]
                
                # Build path with run subdirectory
                dest_path = f"{DEST_BASE}/{projid}/{library_id}/{run_id}/{fastq_filename}"
                run_id_set.add(run_id)
            else:
                # Single source - no run subdirectory needed
                dest_path = f"{DEST_BASE}/{projid}/{library_id}/{fastq_filename}"
            
            # Globus batch format: source_path dest_path
            batch_lines.append(f"{source_path} {dest_path}")
            
            projid_set.add(projid)
            library_id_set.add(library_id)
            total_size += file_size
    
    # Summary statistics
    print("Transfer Summary:")
    print(f"  Total FASTQ files: {len(batch_lines)}")
    print(f"  Unique patients (projids): {len(projid_set)}")
    print(f"  Unique Library_IDs: {len(library_id_set)}")
    print(f"  Unique run IDs (for multi-source samples): {len(run_id_set)}")
    print(f"  Total data size: {total_size / 1e12:.2f} TB ({total_size / 1e9:.1f} GB)")
    print()
    
    # Show sample entries
    print("Sample batch entries (first 5):")
    for line in batch_lines[:5]:
        parts = line.split(' ')
        print(f"  Source: {parts[0][:60]}...")
        print(f"  Dest:   {parts[1]}")
        print()
    
    if dry_run:
        print("[DRY RUN] Skipping file creation")
        return True
    
    # Ensure output directory exists
    BATCH_FILE.parent.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    
    # Write batch file
    with open(BATCH_FILE, 'w') as f:
        for line in batch_lines:
            f.write(line + '\n')
    
    print(f"Batch file written: {BATCH_FILE}")
    print(f"  Lines: {len(batch_lines)}")
    
    # Write summary log
    log_file = LOG_DIR / f"batch_generation_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    with open(log_file, 'w') as f:
        f.write("Globus Batch File Generation Log\n")
        f.write("=" * 50 + "\n")
        f.write(f"Timestamp: {datetime.now().isoformat()}\n")
        f.write(f"Input CSV: {INPUT_CSV}\n")
        f.write(f"Output batch file: {BATCH_FILE}\n")
        f.write(f"Destination base: {DEST_BASE}\n")
        f.write("\n")
        f.write("Statistics:\n")
        f.write(f"  Total FASTQ files: {len(batch_lines)}\n")
        f.write(f"  Unique projids: {len(projid_set)}\n")
        f.write(f"  Unique Library_IDs: {len(library_id_set)}\n")
        f.write(f"  Unique run IDs: {len(run_id_set)}\n")
        f.write(f"  Multi-source samples: {len(multi_source)}\n")
        f.write(f"  Total data size: {total_size / 1e12:.2f} TB\n")
        f.write("\n")
        f.write("Globus Endpoints:\n")
        f.write(f"  Source (Engaging): {ENGAGING_ENDPOINT}\n")
        f.write(f"  Destination (Openmind): {OPENMIND_ENDPOINT}\n")
    
    print(f"Log file written: {log_file}")
    
    print()
    print("=" * 70)
    print("Batch file generation complete!")
    print("=" * 70)
    print()
    print("Next steps:")
    print("  1. Review the batch file: less Batch_Files/globus_batch.txt")
    print("  2. Submit the transfer: ./Scripts/submit_globus_transfer.sh")
    
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Generate Globus batch transfer file for FASTQ files"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print summary without creating files"
    )
    args = parser.parse_args()
    
    success = generate_batch_file(dry_run=args.dry_run)
    return 0 if success else 1


if __name__ == "__main__":
    exit(main())
