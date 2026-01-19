#!/usr/bin/env python3
"""
Verify Globus transfer completion for FASTQ files.

This script compares the expected files (from All_ROSMAP_FASTQs.csv) against
the actual files on Openmind after the Globus transfer completes.

Note: This script uses the same path logic as generate_globus_batch.py,
including run_id subdirectories for multi-source samples.

Usage:
    # On Openmind after transfer:
    python verify_transfer.py
    
    # Generate expected file list only (for pre-transfer):
    python verify_transfer.py --generate-manifest
    
    # Check specific directory:
    python verify_transfer.py --dest-dir /path/to/fastqs

Output:
    - Verification report with missing/extra files
    - Summary statistics
"""

import argparse
import csv
import os
import re
import sys
import hashlib
from pathlib import Path
from datetime import datetime
from collections import defaultdict


# =============================================================================
# Configuration
# =============================================================================

# Input CSV with expected FASTQ locations
REPO_ROOT = "/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis/ROSMAP-SingleNucleusRNAseq"
INPUT_CSV = f"{REPO_ROOT}/Data/Tsai/All_ROSMAP_FASTQs.csv"

# Default destination base path on Openmind
DEFAULT_DEST_BASE = "/om/scratch/Mon/mabdel03/Tsai_Data/FASTQs"

# Script directory (for relative paths to output)
SCRIPT_DIR = Path(__file__).parent.parent
LOG_DIR = SCRIPT_DIR / "Logs"


def extract_run_id(source_dir):
    """
    Extract a unique run identifier from the source directory path.
    Must match the logic in generate_globus_batch.py exactly.
    """
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


def load_expected_files(csv_path, dest_base):
    """Load expected files from CSV and build expected destination paths."""
    
    # First pass: identify source directories per projid/library_id
    projid_lib_sources = defaultdict(set)
    
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            key = (row['projid'], row['Library_ID'])
            projid_lib_sources[key].add(row['source_dir'])
    
    # Identify which projid/library combinations have multiple sources
    multi_source = {k: v for k, v in projid_lib_sources.items() if len(v) > 1}
    
    expected_files = {}  # dest_path -> source_info
    projid_files = defaultdict(list)  # projid -> list of expected files
    source_to_run = {}  # Track source_dir -> run_id mapping
    
    with open(csv_path, 'r') as f:
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
                    # Handle potential run_id collisions
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
                
                dest_path = f"{dest_base}/{projid}/{library_id}/{run_id}/{fastq_filename}"
            else:
                # Single source - no run subdirectory
                dest_path = f"{dest_base}/{projid}/{library_id}/{fastq_filename}"
            
            expected_files[dest_path] = {
                'projid': projid,
                'library_id': library_id,
                'filename': fastq_filename,
                'source_path': source_path,
                'expected_size': file_size
            }
            
            projid_files[projid].append(dest_path)
    
    return expected_files, projid_files


def scan_actual_files(dest_base):
    """Scan the destination directory for actual files.
    
    Handles both structures:
    - projid/Library_ID/filename.fastq.gz (single source)
    - projid/Library_ID/run_id/filename.fastq.gz (multi source)
    """
    
    actual_files = {}  # path -> size
    
    dest_path = Path(dest_base)
    if not dest_path.exists():
        return actual_files
    
    # Walk through all projid directories
    for projid_dir in dest_path.iterdir():
        if not projid_dir.is_dir():
            continue
        
        # Walk through Library_ID subdirectories
        for lib_dir in projid_dir.iterdir():
            if not lib_dir.is_dir():
                continue
            
            # Check contents - could be FASTQ files or run_id subdirectories
            for item in lib_dir.iterdir():
                if item.is_file() and item.suffix == '.gz':
                    # Direct FASTQ file (single source case)
                    actual_files[str(item)] = item.stat().st_size
                elif item.is_dir():
                    # run_id subdirectory (multi source case)
                    for fastq_file in item.iterdir():
                        if fastq_file.is_file() and fastq_file.suffix == '.gz':
                            actual_files[str(fastq_file)] = fastq_file.stat().st_size
    
    return actual_files


def verify_transfer(dest_base, csv_path, verbose=False):
    """Verify the transfer by comparing expected vs actual files."""
    
    print("=" * 70)
    print("Globus Transfer Verification")
    print("=" * 70)
    print()
    
    print(f"Expected files from: {csv_path}")
    print(f"Checking destination: {dest_base}")
    print()
    
    # Load expected files
    print("Loading expected files...")
    expected_files, projid_files = load_expected_files(csv_path, dest_base)
    print(f"  Expected files: {len(expected_files)}")
    print(f"  Unique projids: {len(projid_files)}")
    print()
    
    # Scan actual files
    print("Scanning destination directory...")
    actual_files = scan_actual_files(dest_base)
    print(f"  Found files: {len(actual_files)}")
    print()
    
    # Compare
    expected_set = set(expected_files.keys())
    actual_set = set(actual_files.keys())
    
    missing_files = expected_set - actual_set
    extra_files = actual_set - expected_set
    present_files = expected_set & actual_set
    
    # Check sizes for present files
    size_mismatches = []
    for path in present_files:
        expected_size = expected_files[path]['expected_size']
        actual_size = actual_files[path]
        if expected_size != actual_size:
            size_mismatches.append({
                'path': path,
                'expected': expected_size,
                'actual': actual_size
            })
    
    # Results
    print("=" * 70)
    print("Verification Results")
    print("=" * 70)
    print()
    print(f"Expected files:     {len(expected_files)}")
    print(f"Present files:      {len(present_files)}")
    print(f"Missing files:      {len(missing_files)}")
    print(f"Extra files:        {len(extra_files)}")
    print(f"Size mismatches:    {len(size_mismatches)}")
    print()
    
    # Success check
    if len(missing_files) == 0 and len(size_mismatches) == 0:
        print("SUCCESS: All expected files are present with correct sizes!")
        success = True
    else:
        print("WARNING: Transfer verification found issues")
        success = False
    
    print()
    
    # Report missing files by projid
    if missing_files:
        print("Missing Files by Patient:")
        print("-" * 40)
        
        missing_by_projid = defaultdict(list)
        for path in missing_files:
            info = expected_files[path]
            missing_by_projid[info['projid']].append(info['filename'])
        
        # Show summary
        for projid in sorted(missing_by_projid.keys()):
            files = missing_by_projid[projid]
            print(f"  {projid}: {len(files)} files missing")
            if verbose:
                for f in files[:5]:
                    print(f"    - {f}")
                if len(files) > 5:
                    print(f"    ... and {len(files) - 5} more")
        print()
    
    # Report size mismatches
    if size_mismatches:
        print("Size Mismatches:")
        print("-" * 40)
        for mismatch in size_mismatches[:10]:
            print(f"  {mismatch['path']}")
            print(f"    Expected: {mismatch['expected']}, Actual: {mismatch['actual']}")
        if len(size_mismatches) > 10:
            print(f"  ... and {len(size_mismatches) - 10} more")
        print()
    
    # Write verification report
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    report_file = LOG_DIR / f"verification_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
    
    with open(report_file, 'w') as f:
        f.write("Globus Transfer Verification Report\n")
        f.write("=" * 50 + "\n")
        f.write(f"Timestamp: {datetime.now().isoformat()}\n")
        f.write(f"Destination: {dest_base}\n")
        f.write(f"CSV source: {csv_path}\n")
        f.write("\n")
        f.write("Summary:\n")
        f.write(f"  Expected files: {len(expected_files)}\n")
        f.write(f"  Present files: {len(present_files)}\n")
        f.write(f"  Missing files: {len(missing_files)}\n")
        f.write(f"  Extra files: {len(extra_files)}\n")
        f.write(f"  Size mismatches: {len(size_mismatches)}\n")
        f.write(f"  Status: {'SUCCESS' if success else 'ISSUES FOUND'}\n")
        f.write("\n")
        
        if missing_files:
            f.write("Missing Files:\n")
            f.write("-" * 50 + "\n")
            for path in sorted(missing_files):
                info = expected_files[path]
                f.write(f"{info['projid']}/{info['library_id']}/{info['filename']}\n")
            f.write("\n")
        
        if size_mismatches:
            f.write("Size Mismatches:\n")
            f.write("-" * 50 + "\n")
            for mismatch in size_mismatches:
                f.write(f"{mismatch['path']}\n")
                f.write(f"  Expected: {mismatch['expected']}, Actual: {mismatch['actual']}\n")
            f.write("\n")
    
    print(f"Verification report saved: {report_file}")
    
    return success


def generate_manifest(dest_base, csv_path):
    """Generate a manifest of expected files for pre-transfer verification."""
    
    print("Generating expected file manifest...")
    
    expected_files, projid_files = load_expected_files(csv_path, dest_base)
    
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    manifest_file = LOG_DIR / "expected_files_manifest.txt"
    
    with open(manifest_file, 'w') as f:
        f.write("# Expected FASTQ files after Globus transfer\n")
        f.write(f"# Generated: {datetime.now().isoformat()}\n")
        f.write(f"# Total files: {len(expected_files)}\n")
        f.write(f"# Destination base: {dest_base}\n")
        f.write("#\n")
        f.write("# Format: expected_size destination_path\n")
        f.write("#\n")
        
        for path in sorted(expected_files.keys()):
            info = expected_files[path]
            f.write(f"{info['expected_size']} {path}\n")
    
    print(f"Manifest saved: {manifest_file}")
    print(f"  Total files: {len(expected_files)}")
    print(f"  Unique projids: {len(projid_files)}")
    
    # Also generate per-projid summary
    summary_file = LOG_DIR / "expected_files_by_projid.txt"
    with open(summary_file, 'w') as f:
        f.write("# Expected files per projid\n")
        f.write(f"# Generated: {datetime.now().isoformat()}\n")
        f.write("#\n")
        f.write("# Format: projid file_count\n")
        f.write("#\n")
        
        for projid in sorted(projid_files.keys()):
            f.write(f"{projid} {len(projid_files[projid])}\n")
    
    print(f"Per-projid summary saved: {summary_file}")


def main():
    parser = argparse.ArgumentParser(
        description="Verify Globus transfer completion for FASTQ files"
    )
    parser.add_argument(
        "--dest-dir",
        default=DEFAULT_DEST_BASE,
        help=f"Destination directory to check (default: {DEFAULT_DEST_BASE})"
    )
    parser.add_argument(
        "--csv",
        default=INPUT_CSV,
        help=f"Input CSV with expected files (default: {INPUT_CSV})"
    )
    parser.add_argument(
        "--generate-manifest",
        action="store_true",
        help="Only generate manifest of expected files, don't verify"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Show detailed output"
    )
    
    args = parser.parse_args()
    
    # Check CSV exists
    if not os.path.exists(args.csv):
        print(f"ERROR: CSV file not found: {args.csv}")
        return 1
    
    if args.generate_manifest:
        generate_manifest(args.dest_dir, args.csv)
        return 0
    
    success = verify_transfer(args.dest_dir, args.csv, verbose=args.verbose)
    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
