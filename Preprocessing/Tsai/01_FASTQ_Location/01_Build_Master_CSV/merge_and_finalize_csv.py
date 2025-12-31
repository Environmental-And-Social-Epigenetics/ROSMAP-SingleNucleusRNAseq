#!/usr/bin/env python3
"""
Merge and finalize the master FASTQ CSV.

This script:
1. Reads the combined worker outputs
2. Validates all file paths exist
3. Deduplicates entries
4. Outputs the final All_ROSMAP_FASTQs.csv with statistics
"""

import argparse
import pandas as pd
import os
import sys
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge and finalize FASTQ CSV from worker outputs"
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to combined worker output CSV"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to output the final CSV"
    )
    parser.add_argument(
        "--master-csv",
        required=True,
        help="Path to original Tsai_To_Openmind.csv for comparison"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    
    print(f"Reading combined output from: {args.input}")
    
    # Read the combined CSV
    df = pd.read_csv(args.input)
    
    print(f"Total raw records: {len(df)}")
    
    # Remove duplicates (same full_path)
    df_dedup = df.drop_duplicates(subset=['full_path'], keep='first')
    print(f"After deduplication: {len(df_dedup)}")
    
    # Validate file existence
    print("Validating file existence...")
    valid_mask = df_dedup['full_path'].apply(lambda x: os.path.isfile(x) if pd.notna(x) else False)
    df_valid = df_dedup[valid_mask].copy()
    df_invalid = df_dedup[~valid_mask].copy()
    
    print(f"Valid files: {len(df_valid)}")
    print(f"Invalid/missing files: {len(df_invalid)}")
    
    # Sort by projid, then Library_ID, then filename
    df_valid = df_valid.sort_values(['projid', 'Library_ID', 'fastq_filename'])
    
    # Calculate statistics
    n_projids = df_valid['projid'].nunique()
    n_library_ids = df_valid['Library_ID'].nunique()
    total_size_gb = df_valid['file_size'].sum() / (1024**3)
    
    print("\n" + "="*50)
    print("SUMMARY STATISTICS")
    print("="*50)
    print(f"Total FASTQ files indexed: {len(df_valid)}")
    print(f"Unique projids: {n_projids}")
    print(f"Unique Library_IDs: {n_library_ids}")
    print(f"Total file size: {total_size_gb:.2f} GB")
    print("="*50 + "\n")
    
    # Compare with master CSV
    print("Comparing with master CSV...")
    master_df = pd.read_csv(args.master_csv)
    master_projids = set(master_df['projid'].astype(str).unique())
    found_projids = set(df_valid['projid'].astype(str).unique())
    
    missing_projids = master_projids - found_projids
    if missing_projids:
        print(f"WARNING: {len(missing_projids)} projids from master CSV have no FASTQ files:")
        for pid in sorted(missing_projids)[:10]:
            print(f"  - {pid}")
        if len(missing_projids) > 10:
            print(f"  ... and {len(missing_projids) - 10} more")
    else:
        print("All projids from master CSV have at least one FASTQ file.")
    
    # Save the final CSV
    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    df_valid.to_csv(args.output, index=False)
    print(f"\nFinal CSV saved to: {args.output}")
    
    # Save invalid entries for debugging
    if len(df_invalid) > 0:
        invalid_path = args.output.replace('.csv', '_missing.csv')
        df_invalid.to_csv(invalid_path, index=False)
        print(f"Missing files log saved to: {invalid_path}")
    
    # Print sample of the output
    print("\nSample of output:")
    print(df_valid.head(10).to_string())
    
    return 0


if __name__ == "__main__":
    sys.exit(main())

