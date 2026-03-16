#!/usr/bin/env python3
"""
Derive real sequencing batch groups from FASTQ headers.

For each sample, this script:
  1. Locates the R1 FASTQ file(s) across all FASTQ directories
  2. Reads the first line to extract the Illumina flowcell ID
     (format: @instrument:run:flowcellID:lane:tile:x:y)
  3. Groups samples by their set of flowcell IDs
  4. Assigns a human-readable derived_batch label

Output: Resources/derived_batches.csv with columns:
    projid, flowcell_ids, instrument, derived_batch

Usage:
    python derive_batches.py [--metadata-csv PATH] [--output PATH]
"""

import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path

try:
    from typing import Dict, List, Optional, Set, Tuple
except ImportError:
    pass


def default_paths():
    # type: () -> Dict[str, Path]
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parents[2]
    return {
        "metadata_csv": repo_root
        / "Preprocessing"
        / "Tsai"
        / "02_Cellranger_Counts"
        / "Tracking"
        / "patient_metadata.csv",
        "output": script_dir / "Resources" / "derived_batches.csv",
    }


def parse_args():
    # type: () -> argparse.Namespace
    paths = default_paths()
    parser = argparse.ArgumentParser(description="Derive sequencing batch groups from FASTQ headers.")
    parser.add_argument("--metadata-csv", type=Path, default=paths["metadata_csv"])
    parser.add_argument("--output", type=Path, default=paths["output"])
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print summary without writing output file.",
    )
    return parser.parse_args()


def extract_flowcell_from_fastq(fastq_path):
    # type: (Path) -> Optional[Tuple[str, str]]
    """Read the first line of a FASTQ file and extract (instrument, flowcell_id).

    Illumina header format: @instrument:run:flowcellID:lane:tile:x:y
    """
    try:
        with gzip.open(str(fastq_path), "rt") as gz:
            header = gz.readline().strip()
        parts = header.split(":")
        if len(parts) < 3:
            return None
        instrument = parts[0].lstrip("@")
        flowcell_id = parts[2]
        return instrument, flowcell_id
    except Exception as exc:
        print("  [WARN] Cannot read {}: {}".format(fastq_path, exc))
        return None


def find_r1_fastqs(fastq_dir):
    # type: (Path) -> List[Path]
    """Find R1 FASTQ files in a directory."""
    return sorted(fastq_dir.glob("*_R1_*.fastq.gz"))


def extract_sample_flowcells(projid, fastq_dirs):
    # type: (str, List[str]) -> Tuple[Set[str], Set[str]]
    """Extract all unique (instrument, flowcell) pairs for a sample.

    Returns (flowcell_set, instrument_set).
    """
    flowcells = set()  # type: Set[str]
    instruments = set()  # type: Set[str]

    for fd in fastq_dirs:
        fd_path = Path(fd)
        if not fd_path.exists():
            continue

        r1_files = find_r1_fastqs(fd_path)
        if not r1_files:
            # Try one level deeper (some dirs have run subdirectories)
            for sub in sorted(fd_path.iterdir()):
                if sub.is_dir():
                    r1_files.extend(find_r1_fastqs(sub))

        if not r1_files:
            continue

        # Only need the first R1 per directory (all share the same flowcell)
        result = extract_flowcell_from_fastq(r1_files[0])
        if result is not None:
            instrument, flowcell_id = result
            flowcells.add(flowcell_id)
            instruments.add(instrument)

    return flowcells, instruments


def assign_batch_labels(sample_flowcells):
    # type: (Dict[str, Tuple[Set[str], Set[str]]]) -> Dict[str, Dict[str, str]]
    """Group samples by their flowcell signature and assign batch labels.

    Returns dict[projid] -> {flowcell_ids, instrument, derived_batch}.
    """
    # Build flowcell signature -> list of projids
    sig_to_samples = defaultdict(list)  # type: Dict[Tuple[str, ...], List[str]]
    for projid, (flowcells, _instruments) in sample_flowcells.items():
        sig = tuple(sorted(flowcells))
        sig_to_samples[sig].append(projid)

    # For samples with complex signatures (3+ flowcells, i.e. re-sequenced),
    # find the best matching 2-flowcell group to merge them into.
    two_fc_sigs = {sig for sig in sig_to_samples if len(sig) == 2}
    one_fc_sigs = {sig for sig in sig_to_samples if len(sig) == 1}
    complex_sigs = {sig for sig in sig_to_samples if len(sig) > 2}

    # Map complex signatures to the 2-flowcell group they overlap with most
    reassignments = {}  # type: Dict[str, Tuple[str, ...]]
    for sig in complex_sigs:
        fc_set = set(sig)
        best_match = None
        best_overlap = 0
        for two_sig in two_fc_sigs:
            overlap = len(fc_set & set(two_sig))
            if overlap > best_overlap:
                best_overlap = overlap
                best_match = two_sig
        if best_match is None:
            for one_sig in one_fc_sigs:
                if set(one_sig) & fc_set:
                    best_match = one_sig
                    break
        if best_match is not None:
            for projid in sig_to_samples[sig]:
                reassignments[projid] = best_match

    # Apply reassignments
    final_groups = defaultdict(list)  # type: Dict[Tuple[str, ...], List[str]]
    for sig, projids in sig_to_samples.items():
        if len(sig) <= 2:
            final_groups[sig].extend(projids)
    for projid, target_sig in reassignments.items():
        final_groups[target_sig].append(projid)
    # Any remaining unmerged complex samples go into their own group
    for sig in complex_sigs:
        for projid in sig_to_samples[sig]:
            if projid not in reassignments:
                final_groups[sig].append(projid)

    # Assign labels: FC_{flowcell1} or FC_{flowcell1}+{flowcell2}
    results = {}  # type: Dict[str, Dict[str, str]]
    for sig in sorted(final_groups, key=lambda s: (-len(final_groups[s]), s)):
        label = "FC_" + "+".join(sig)
        for projid in final_groups[sig]:
            flowcells, instruments = sample_flowcells[projid]
            results[projid] = {
                "flowcell_ids": "+".join(sorted(flowcells)),
                "instrument": "+".join(sorted(instruments)),
                "derived_batch": label,
            }

    return results


def main():
    # type: () -> None
    args = parse_args()

    print("Reading metadata: {}".format(args.metadata_csv))
    with open(str(args.metadata_csv)) as f:
        rows = list(csv.DictReader(f))
    print("  {} samples found".format(len(rows)))

    # Extract flowcell IDs for each sample
    print("Extracting flowcell IDs from FASTQ headers...")
    sample_flowcells = {}  # type: Dict[str, Tuple[Set[str], Set[str]]]
    failed = []  # type: List[str]

    for i, row in enumerate(rows):
        projid = row["projid"]
        fastq_dirs = row["fastq_dirs"].split("|")

        flowcells, instruments = extract_sample_flowcells(projid, fastq_dirs)
        if not flowcells:
            print("  [WARN] {}: no flowcell extracted — skipping".format(projid))
            failed.append(projid)
            continue

        sample_flowcells[projid] = (flowcells, instruments)

        if (i + 1) % 50 == 0:
            print("  processed {}/{} samples...".format(i + 1, len(rows)))

    print("\nExtraction complete: {} succeeded, {} failed".format(
        len(sample_flowcells), len(failed)))

    # Assign batch labels
    results = assign_batch_labels(sample_flowcells)

    # Print summary
    batch_counts = defaultdict(int)
    for info in results.values():
        batch_counts[info["derived_batch"]] += 1

    print("\nDerived {} batch groups:".format(len(batch_counts)))
    for batch, count in sorted(batch_counts.items(), key=lambda x: -x[1]):
        print("  {}: {} samples".format(batch, count))

    if args.dry_run:
        print("\n[DRY RUN] Skipping file output")
        return

    # Write output
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["projid", "flowcell_ids", "instrument", "derived_batch"]
    with open(str(args.output), "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for projid in sorted(results):
            writer.writerow(dict(projid=projid, **results[projid]))

    print("\nWrote {} rows to {}".format(len(results), args.output))

    if failed:
        print("\n[WARN] {} samples had no extractable flowcell ID:".format(len(failed)))
        for projid in failed:
            print("  {}".format(projid))


if __name__ == "__main__":
    main()
