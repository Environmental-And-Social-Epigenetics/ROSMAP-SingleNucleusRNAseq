#!/usr/bin/env python3
"""
Validate the organized FASTQ directory structure.

This script:
1. Verifies all expected symlinks exist and are not broken
2. Compares file counts against the master CSV
3. Outputs a detailed validation report
"""

import argparse
import pandas as pd
import os
import sys
from pathlib import Path
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        description="Validate organized FASTQ directory structure"
    )
    parser.add_argument(
        "--fastq-csv",
        required=True,
        help="Path to All_ROSMAP_FASTQs.csv"
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Path to FASTQs_By_Patient directory"
    )
    parser.add_argument(
        "--report",
        required=True,
        help="Path to output validation report"
    )
    return parser.parse_args()


def count_symlinks(directory):
    """Count all symlinks in a directory recursively."""
    count = 0
    broken = 0
    for root, dirs, files in os.walk(directory):
        for f in files:
            filepath = os.path.join(root, f)
            if os.path.islink(filepath):
                count += 1
                if not os.path.exists(filepath):
                    broken += 1
    return count, broken


def main():
    args = parse_args()
    
    report_lines = []
    
    def log(msg):
        print(msg)
        report_lines.append(msg)
    
    log("=" * 60)
    log("FASTQ Organization Validation Report")
    log("=" * 60)
    log(f"Date: {pd.Timestamp.now()}")
    log(f"Master CSV: {args.fastq_csv}")
    log(f"Output directory: {args.output_dir}")
    log("")
    
    # Load the master FASTQ CSV
    log("Loading master FASTQ CSV...")
    df = pd.read_csv(args.fastq_csv)
    
    expected_files = len(df)
    expected_projids = df['projid'].nunique()
    expected_library_ids = df['Library_ID'].nunique()
    
    log(f"Expected FASTQ files: {expected_files}")
    log(f"Expected projids: {expected_projids}")
    log(f"Expected Library_IDs: {expected_library_ids}")
    log("")
    
    # Check output directory exists
    if not os.path.isdir(args.output_dir):
        log(f"ERROR: Output directory does not exist: {args.output_dir}")
        log("Validation FAILED")
        with open(args.report, 'w') as f:
            f.write('\n'.join(report_lines))
        return 1
    
    # Count patient directories
    patient_dirs = [d for d in os.listdir(args.output_dir) 
                    if os.path.isdir(os.path.join(args.output_dir, d))]
    log(f"Patient directories found: {len(patient_dirs)}")
    
    # Count symlinks
    log("Counting symlinks (this may take a moment)...")
    total_links, broken_links = count_symlinks(args.output_dir)
    
    log(f"Total symlinks created: {total_links}")
    log(f"Broken symlinks: {broken_links}")
    log("")
    
    # Per-projid validation
    log("-" * 60)
    log("Per-Patient Validation")
    log("-" * 60)
    
    issues = []
    projid_stats = defaultdict(dict)
    
    for projid in df['projid'].unique():
        projid_str = str(projid)
        projid_dir = os.path.join(args.output_dir, projid_str)
        
        expected_for_projid = len(df[df['projid'] == projid])
        
        if not os.path.isdir(projid_dir):
            issues.append(f"Missing directory for projid {projid}")
            projid_stats[projid]['status'] = 'MISSING'
            projid_stats[projid]['expected'] = expected_for_projid
            projid_stats[projid]['found'] = 0
            continue
        
        actual_links, actual_broken = count_symlinks(projid_dir)
        
        projid_stats[projid]['expected'] = expected_for_projid
        projid_stats[projid]['found'] = actual_links
        projid_stats[projid]['broken'] = actual_broken
        
        if actual_links < expected_for_projid:
            issues.append(f"projid {projid}: Expected {expected_for_projid}, found {actual_links}")
            projid_stats[projid]['status'] = 'INCOMPLETE'
        elif actual_broken > 0:
            issues.append(f"projid {projid}: {actual_broken} broken symlinks")
            projid_stats[projid]['status'] = 'BROKEN_LINKS'
        else:
            projid_stats[projid]['status'] = 'OK'
    
    # Summary statistics
    ok_count = sum(1 for p, s in projid_stats.items() if s.get('status') == 'OK')
    incomplete_count = sum(1 for p, s in projid_stats.items() if s.get('status') == 'INCOMPLETE')
    missing_count = sum(1 for p, s in projid_stats.items() if s.get('status') == 'MISSING')
    broken_count = sum(1 for p, s in projid_stats.items() if s.get('status') == 'BROKEN_LINKS')
    
    log(f"Patients OK: {ok_count}")
    log(f"Patients incomplete: {incomplete_count}")
    log(f"Patients missing: {missing_count}")
    log(f"Patients with broken links: {broken_count}")
    log("")
    
    # List issues (first 50)
    if issues:
        log("-" * 60)
        log("Issues Found (first 50):")
        log("-" * 60)
        for issue in issues[:50]:
            log(f"  - {issue}")
        if len(issues) > 50:
            log(f"  ... and {len(issues) - 50} more issues")
        log("")
    
    # Final verdict
    log("=" * 60)
    if total_links >= expected_files * 0.95 and broken_links == 0:
        log("VALIDATION: PASSED")
        log(f"Coverage: {total_links}/{expected_files} ({100*total_links/expected_files:.1f}%)")
    elif total_links >= expected_files * 0.80:
        log("VALIDATION: PASSED WITH WARNINGS")
        log(f"Coverage: {total_links}/{expected_files} ({100*total_links/expected_files:.1f}%)")
        if broken_links > 0:
            log(f"Warning: {broken_links} broken symlinks detected")
    else:
        log("VALIDATION: FAILED")
        log(f"Coverage: {total_links}/{expected_files} ({100*total_links/expected_files:.1f}%)")
    log("=" * 60)
    
    # Write report
    with open(args.report, 'w') as f:
        f.write('\n'.join(report_lines))
    
    log(f"\nReport saved to: {args.report}")
    
    return 0 if (total_links >= expected_files * 0.80 and broken_links == 0) else 1


if __name__ == "__main__":
    sys.exit(main())

