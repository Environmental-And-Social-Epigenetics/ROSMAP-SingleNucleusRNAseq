"""
Aggregate demuxlet results into a single cell-to-patient assignment CSV.

Iterates over all library directories in the WGS output directory, reads each
demux1.best file, extracts the cell barcode and best patient assignment, and
writes a combined CSV.

Usage:
    source config/paths.sh
    python postprocess_assignments.py
    python postprocess_assignments.py --wgs-dir /path/to/WGS --output assignments.csv
"""

import argparse
import csv
import os
import sys

import pandas as pd


def aggregate_assignments(wgs_dir, output_path):
    """Read demux1.best from each library and write combined CSV."""
    rows = []
    libraries_processed = 0
    libraries_skipped = 0

    for lib_id in sorted(os.listdir(wgs_dir)):
        lib_path = os.path.join(wgs_dir, lib_id)
        if not os.path.isdir(lib_path):
            continue
        if "alone" in lib_id:
            libraries_skipped += 1
            continue

        demux_file = os.path.join(lib_path, "demux1.best")
        if not os.path.isfile(demux_file):
            continue

        df = pd.read_csv(demux_file, sep="\t")
        # Extract first patient ID from comma-separated BEST.GUESS
        df["BEST.GUESS"] = df["BEST.GUESS"].str.split(",").str[0]

        for _, row in df.iterrows():
            rows.append({
                "Cell Barcode": row["BARCODE"],
                "Assigned Patient": row["BEST.GUESS"],
                "Library": lib_id,
            })
        libraries_processed += 1

    if not rows:
        print("No demux1.best files found. Run demuxlet first.")
        sys.exit(1)

    # Write output
    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["Cell Barcode", "Assigned Patient", "Library"])
        writer.writeheader()
        writer.writerows(rows)

    # Summary
    df_out = pd.DataFrame(rows)
    print(f"Libraries processed: {libraries_processed}")
    print(f"Libraries skipped ('alone'): {libraries_skipped}")
    print(f"Total cells assigned: {len(rows)}")
    print(f"Unique patients: {df_out['Assigned Patient'].nunique()}")
    print(f"\nCells per library:")
    for lib, count in df_out.groupby("Library").size().items():
        print(f"  {lib}: {count}")
    print(f"\nOutput written to: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Aggregate demuxlet results into cell-to-patient assignment CSV."
    )
    parser.add_argument(
        "--wgs-dir",
        default=os.environ.get("DEJAGER_WGS_DIR", "/om/scratch/Mon/shared_folder/WGS"),
        help="Directory containing per-library demuxlet outputs (default: $DEJAGER_WGS_DIR)",
    )
    parser.add_argument(
        "--output", "-o",
        default="cell_to_patient_assignments.csv",
        help="Output CSV path (default: cell_to_patient_assignments.csv)",
    )
    args = parser.parse_args()

    if not os.path.isdir(args.wgs_dir):
        print(f"ERROR: WGS directory not found: {args.wgs_dir}")
        sys.exit(1)

    aggregate_assignments(args.wgs_dir, args.output)


if __name__ == "__main__":
    main()
