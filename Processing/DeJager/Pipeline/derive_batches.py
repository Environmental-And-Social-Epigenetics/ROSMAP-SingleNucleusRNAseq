#!/usr/bin/env python3
"""
Derive real sequencing batch groups from Cell Ranger BAM headers and FASTQ headers.

For each DeJager library, this script:
  1. Reads the BAM @RG (read group) header from the Cell Ranger output
  2. Extracts the Illumina flowcell ID from the PU (Platform Unit) tag
  3. Infers flowcells for A/B sibling libraries that lack BAMs
  4. Optionally downloads one R1 FASTQ per missing B-number from Synapse
     to fill in remaining gaps (--fill-from-synapse)
  5. Groups libraries by their set of flowcell IDs
  6. Assigns a human-readable derived_batch label

Output: Resources/derived_batches.csv with columns:
    library_id, flowcell_ids, instrument, derived_batch, source

Usage:
    # BAM-only (with pool fallbacks for missing B-numbers):
    python derive_batches.py

    # Full coverage via Synapse FASTQ download:
    python derive_batches.py --fill-from-synapse
"""

import argparse
import csv
import gzip
import os
import re
import struct
import tempfile
from collections import defaultdict
from pathlib import Path

try:
    from typing import Dict, List, Optional, Set, Tuple
except ImportError:
    pass


# Pattern matching DeJager library naming convention: YYMMDD-BXX[-suffix]
_LIBRARY_RE = re.compile(r"^\d{6}-B\d+")
# Extract the B-number (pool prefix) from a library name
_BNUMBER_RE = re.compile(r"^\d{6}-(B\d+)")


def default_paths():
    # type: () -> Dict[str, Path]
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parents[2]
    workspace_root = repo_root.parent
    return {
        "cellranger_dir": Path(
            os.environ.get(
                "DEJAGER_CELLRANGER",
                str(workspace_root / "DeJager_Data" / "Cellranger_Output"),
            )
        ),
        "cellbender_dir": Path(
            os.environ.get(
                "DEJAGER_PREPROCESSED",
                str(workspace_root / "DeJager_Data" / "Cellbender_Output"),
            )
        ),
        "output": script_dir / "Resources" / "derived_batches.csv",
        "synapse_csv": Path(
            os.environ.get(
                "DEJAGER_SYNAPSE_FASTQ_CSV",
                "/orcd/data/lhtsai/001/om2/mabdel03/files/ACE_Analysis"
                "/Data/DeJager/FASTQs_Download"
                "/FASTQ_Download_CSVs/Synapse_FASTQ_IDs.csv",
            )
        ),
    }


def parse_args():
    # type: () -> argparse.Namespace
    paths = default_paths()
    parser = argparse.ArgumentParser(
        description="Derive sequencing batch groups from Cell Ranger BAM headers."
    )
    parser.add_argument("--cellranger-dir", type=Path, default=paths["cellranger_dir"])
    parser.add_argument("--cellbender-dir", type=Path, default=paths["cellbender_dir"])
    parser.add_argument("--output", type=Path, default=paths["output"])
    parser.add_argument(
        "--fill-from-synapse",
        action="store_true",
        help="Download one R1 FASTQ per missing B-number from Synapse to "
        "fill in flowcell IDs that cannot be determined from BAMs.",
    )
    parser.add_argument(
        "--synapse-csv",
        type=Path,
        default=paths["synapse_csv"],
        help="CSV mapping SynID -> FileName -> LibraryID for Synapse FASTQ files.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print summary without writing output file.",
    )
    return parser.parse_args()


def extract_bnumber(library_id):
    # type: (str) -> Optional[str]
    """Extract pool prefix (e.g., 'B12') from a library name."""
    match = _BNUMBER_RE.match(library_id)
    return match.group(1) if match else None


def discover_cellbender_libraries(cellbender_dir):
    # type: (Path) -> List[str]
    """Return sorted list of library IDs from the CellBender output directory."""
    libraries = []
    for entry in sorted(cellbender_dir.iterdir()):
        if entry.is_dir() and _LIBRARY_RE.match(entry.name):
            libraries.append(entry.name)
    return libraries


def discover_cellranger_bams(cellranger_dir):
    # type: (Path) -> Dict[str, Path]
    """Return dict mapping library_id -> BAM path for available Cell Ranger outputs."""
    bams = {}
    if not cellranger_dir.exists():
        return bams
    for entry in sorted(cellranger_dir.iterdir()):
        if not entry.is_dir() or not _LIBRARY_RE.match(entry.name):
            continue
        bam_path = entry / "outs" / "possorted_genome_bam.bam"
        if bam_path.exists():
            bams[entry.name] = bam_path
    return bams


def extract_flowcells_from_bam(bam_path):
    # type: (Path) -> Tuple[Set[str], Set[str]]
    """Extract flowcell IDs and instrument IDs from BAM @RG headers.

    Reads the binary BAM header directly (no pysam dependency).
    BAM format: magic(4) + header_len(4) + header_text(header_len) + ...
    """
    flowcells = set()  # type: Set[str]
    instruments = set()  # type: Set[str]

    try:
        with open(str(bam_path), "rb") as f:
            data = gzip.GzipFile(fileobj=f).read(200000)

        if data[:4] != b"BAM\1":
            print("  [WARN] Not a valid BAM file: {}".format(bam_path))
            return flowcells, instruments

        header_len = struct.unpack("<I", data[4:8])[0]
        header_text = data[8 : 8 + header_len].decode("ascii", errors="replace")

        for line in header_text.split("\n"):
            if not line.startswith("@RG"):
                continue
            for field in line.split("\t"):
                if field.startswith("PU:"):
                    # PU format: SAMPLE:run:lane:FLOWCELL:LANE
                    parts = field[3:].split(":")
                    if len(parts) >= 5:
                        flowcells.add(parts[3])
                    elif len(parts) >= 2:
                        # Fallback: FLOWCELL:LANE
                        flowcells.add(parts[0])
                elif field.startswith("PM:"):
                    instruments.add(field[3:])

        # Also try extracting instrument from ID field if PM not present
        if not instruments:
            for line in header_text.split("\n"):
                if not line.startswith("@RG"):
                    continue
                for field in line.split("\t"):
                    if field.startswith("ID:"):
                        # ID format: SAMPLE:run:lane:FLOWCELL:LANE
                        parts = field[3:].split(":")
                        if len(parts) >= 5:
                            # Instrument is not directly in ID; skip
                            pass

    except Exception as exc:
        print("  [WARN] Cannot read {}: {}".format(bam_path, exc))

    return flowcells, instruments


def extract_flowcell_from_fastq(fastq_path):
    # type: (Path) -> Optional[Tuple[str, str]]
    """Read the first line of a gzipped FASTQ and extract (instrument, flowcell_id).

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


def find_synapse_r1_for_missing(synapse_csv, missing_bnumbers):
    # type: (Path, Set[str]) -> Dict[str, List[Tuple[str, str, str]]]
    """Find R1 FASTQ SynIDs for missing B-numbers from the Synapse CSV.

    Downloads one R1 per unique sequencing center (Broad vs NYGC) per B-number,
    since different centers use different flowcells.

    Returns dict[bnumber] -> [(syn_id, library_id, filename), ...].
    """
    # Collect all R1 files per B-number, keyed by sequencing center
    bnumber_center_r1 = defaultdict(
        dict
    )  # type: Dict[str, Dict[str, Tuple[str, str, str]]]
    with open(str(synapse_csv)) as f:
        reader = csv.DictReader(f)
        for row in reader:
            lib = row["LibraryID"]
            fname = row["FileName"]
            m = _BNUMBER_RE.match(lib)
            if not m:
                continue
            bn = m.group(1)
            if bn not in missing_bnumbers or "_R1_" not in fname:
                continue
            # Determine sequencing center from filename (e.g., _Broad_, _NYGC_)
            if "_Broad_" in fname:
                center = "Broad"
            elif "_NYGC_" in fname:
                center = "NYGC"
            else:
                center = "unknown"
            # Keep only one R1 per center per B-number
            if center not in bnumber_center_r1[bn]:
                bnumber_center_r1[bn][center] = (row["SynID"], lib, fname)

    # Flatten: one entry per center per B-number
    result = {}  # type: Dict[str, List[Tuple[str, str, str]]]
    for bn, centers in bnumber_center_r1.items():
        result[bn] = list(centers.values())
    return result


def download_and_extract_flowcells(bnumber_to_r1s):
    # type: (Dict[str, List[Tuple[str, str, str]]]) -> Dict[str, Tuple[Set[str], Set[str]]]
    """Download R1 FASTQs per B-number from Synapse, extract flowcells, delete.

    Downloads one R1 per sequencing center per B-number to capture all flowcells.

    Returns dict[bnumber] -> (flowcell_set, instrument_set).
    """
    try:
        import synapseclient
    except ImportError:
        raise SystemExit(
            "synapseclient is required for --fill-from-synapse.\n"
            "Activate the Synapse conda environment first."
        )

    syn = synapseclient.Synapse()
    token = os.environ.get("SYNAPSE_AUTH_TOKEN", "")
    if token:
        syn.login(authToken=token, silent=True)
    else:
        syn.login(silent=True)
    print("  Synapse login successful")

    results = {}  # type: Dict[str, Tuple[Set[str], Set[str]]]

    with tempfile.TemporaryDirectory(prefix="dejager_fc_") as tmp_dir:
        for bn in sorted(bnumber_to_r1s):
            r1_list = bnumber_to_r1s[bn]
            flowcells = set()  # type: Set[str]
            instruments = set()  # type: Set[str]
            for syn_id, lib_id, fname in r1_list:
                print("  {} ({}): downloading {}...".format(bn, lib_id, syn_id))
                try:
                    entity = syn.get(syn_id, downloadLocation=tmp_dir)
                    fastq_path = Path(entity.path)
                    result = extract_flowcell_from_fastq(fastq_path)
                    # Delete immediately to save space
                    fastq_path.unlink(missing_ok=True)
                    if result is not None:
                        instrument, flowcell_id = result
                        flowcells.add(flowcell_id)
                        instruments.add(instrument)
                        print("    -> flowcell={}, instrument={}".format(
                            flowcell_id, instrument))
                    else:
                        print("    [WARN] Could not parse FASTQ header")
                except Exception as exc:
                    print("    [WARN] Download failed: {}".format(exc))
            if flowcells:
                results[bn] = (flowcells, instruments)

    return results


def assign_batch_labels(library_flowcells):
    # type: (Dict[str, Tuple[Set[str], Set[str], str]]) -> Dict[str, Dict[str, str]]
    """Group libraries by flowcell signature and assign batch labels.

    Returns dict[library_id] -> {flowcell_ids, instrument, derived_batch, source}.
    """
    # Separate flowcell-based and fallback libraries
    fc_libraries = {}  # type: Dict[str, Tuple[Set[str], Set[str], str]]
    fallback_libraries = {}  # type: Dict[str, str]

    for lib_id, (flowcells, instruments, source) in library_flowcells.items():
        if source == "fallback":
            bnumber = extract_bnumber(lib_id)
            fallback_libraries[lib_id] = bnumber or lib_id
        else:
            fc_libraries[lib_id] = (flowcells, instruments, source)

    # Group flowcell-based libraries by signature
    sig_to_libs = defaultdict(list)  # type: Dict[Tuple[str, ...], List[str]]
    for lib_id, (flowcells, _instruments, _source) in fc_libraries.items():
        sig = tuple(sorted(flowcells))
        sig_to_libs[sig].append(lib_id)

    # For libraries with complex signatures (3+ flowcells),
    # find the best matching 2-flowcell group to merge them into.
    two_fc_sigs = {sig for sig in sig_to_libs if len(sig) == 2}
    one_fc_sigs = {sig for sig in sig_to_libs if len(sig) == 1}
    complex_sigs = {sig for sig in sig_to_libs if len(sig) > 2}

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
            for lib_id in sig_to_libs[sig]:
                reassignments[lib_id] = best_match

    # Build final flowcell groups
    final_groups = defaultdict(list)  # type: Dict[Tuple[str, ...], List[str]]
    for sig, lib_ids in sig_to_libs.items():
        if len(sig) <= 2:
            final_groups[sig].extend(lib_ids)
    for lib_id, target_sig in reassignments.items():
        final_groups[target_sig].append(lib_id)
    for sig in complex_sigs:
        for lib_id in sig_to_libs[sig]:
            if lib_id not in reassignments:
                final_groups[sig].append(lib_id)

    # Assign flowcell-based labels
    results = {}  # type: Dict[str, Dict[str, str]]
    for sig in sorted(final_groups, key=lambda s: (-len(final_groups[s]), s)):
        label = "FC_" + "+".join(sig)
        for lib_id in final_groups[sig]:
            flowcells, instruments, source = fc_libraries[lib_id]
            results[lib_id] = {
                "flowcell_ids": "+".join(sorted(flowcells)),
                "instrument": "+".join(sorted(instruments)),
                "derived_batch": label,
                "source": source,
            }

    # Assign pool-based fallback labels
    for lib_id, bnumber in fallback_libraries.items():
        results[lib_id] = {
            "flowcell_ids": "",
            "instrument": "",
            "derived_batch": "POOL_{}".format(bnumber),
            "source": "fallback",
        }

    return results


def main():
    # type: () -> None
    args = parse_args()

    print("Discovering libraries...")
    all_libraries = discover_cellbender_libraries(args.cellbender_dir)
    print("  {} CellBender libraries found".format(len(all_libraries)))

    bam_map = discover_cellranger_bams(args.cellranger_dir)
    print("  {} Cell Ranger BAMs found".format(len(bam_map)))

    # Extract flowcells from available BAMs
    print("\nExtracting flowcell IDs from BAM headers...")
    bam_flowcells = {}  # type: Dict[str, Tuple[Set[str], Set[str]]]
    for i, (lib_id, bam_path) in enumerate(sorted(bam_map.items())):
        flowcells, instruments = extract_flowcells_from_bam(bam_path)
        if flowcells:
            bam_flowcells[lib_id] = (flowcells, instruments)
            print("  {}: flowcells={}".format(lib_id, sorted(flowcells)))
        else:
            print("  [WARN] {}: no flowcell extracted".format(lib_id))

    # Build B-number -> flowcell mapping for inference
    bnumber_flowcells = defaultdict(
        lambda: (set(), set())
    )  # type: Dict[str, Tuple[Set[str], Set[str]]]
    for lib_id, (flowcells, instruments) in bam_flowcells.items():
        bn = extract_bnumber(lib_id)
        if bn:
            existing_fc, existing_inst = bnumber_flowcells[bn]
            existing_fc.update(flowcells)
            existing_inst.update(instruments)
            bnumber_flowcells[bn] = (existing_fc, existing_inst)

    # Assign flowcells to all libraries
    print("\nAssigning flowcell IDs to all libraries...")
    library_flowcells = (
        {}
    )  # type: Dict[str, Tuple[Set[str], Set[str], str]]  # lib -> (fc, inst, source)

    n_bam = 0
    n_inferred = 0
    n_fallback = 0

    for lib_id in all_libraries:
        if lib_id in bam_flowcells:
            fc, inst = bam_flowcells[lib_id]
            library_flowcells[lib_id] = (fc, inst, "bam")
            n_bam += 1
        else:
            bn = extract_bnumber(lib_id)
            if bn and bn in bnumber_flowcells:
                fc, inst = bnumber_flowcells[bn]
                library_flowcells[lib_id] = (set(fc), set(inst), "inferred")
                n_inferred += 1
                print(
                    "  {} -> inferred from {} (flowcells={})".format(
                        lib_id, bn, sorted(fc)
                    )
                )
            else:
                library_flowcells[lib_id] = (set(), set(), "fallback")
                n_fallback += 1
                print("  {} -> fallback (POOL_{})".format(lib_id, bn or lib_id))

    print(
        "\nSources: {} from BAM, {} inferred, {} fallback".format(
            n_bam, n_inferred, n_fallback
        )
    )

    # Fill in fallback libraries from Synapse FASTQ headers if requested.
    if args.fill_from_synapse and n_fallback > 0:
        missing_bnumbers = set()
        for lib_id, (_fc, _inst, source) in library_flowcells.items():
            if source == "fallback":
                bn = extract_bnumber(lib_id)
                if bn:
                    missing_bnumbers.add(bn)

        print("\nFilling {} missing B-numbers from Synapse FASTQ headers...".format(
            len(missing_bnumbers)))

        if not args.synapse_csv.exists():
            raise SystemExit(
                "Synapse CSV not found: {}\n"
                "Set DEJAGER_SYNAPSE_FASTQ_CSV or pass --synapse-csv".format(
                    args.synapse_csv))

        bnumber_to_r1 = find_synapse_r1_for_missing(args.synapse_csv, missing_bnumbers)
        not_found = missing_bnumbers - set(bnumber_to_r1)
        if not_found:
            print("  [WARN] No R1 FASTQ found in Synapse CSV for: {}".format(
                sorted(not_found)))

        synapse_flowcells = download_and_extract_flowcells(bnumber_to_r1)

        # Merge Synapse results into bnumber_flowcells and reassign
        for bn, (fc, inst) in synapse_flowcells.items():
            bnumber_flowcells[bn] = (fc, inst)

        # Reassign all fallback libraries that now have flowcell data
        n_synapse = 0
        for lib_id in list(library_flowcells):
            _fc, _inst, source = library_flowcells[lib_id]
            if source != "fallback":
                continue
            bn = extract_bnumber(lib_id)
            if bn and bn in synapse_flowcells:
                fc, inst = synapse_flowcells[bn]
                library_flowcells[lib_id] = (set(fc), set(inst), "synapse")
                n_synapse += 1
                n_fallback -= 1

        print("\nUpdated sources: {} from BAM, {} inferred, {} synapse, {} fallback".format(
            n_bam, n_inferred, n_synapse, n_fallback))

    # Assign batch labels
    results = assign_batch_labels(library_flowcells)

    # Print summary
    batch_counts = defaultdict(int)  # type: Dict[str, int]
    for info in results.values():
        batch_counts[info["derived_batch"]] += 1

    print("\nDerived {} batch groups:".format(len(batch_counts)))
    for batch, count in sorted(batch_counts.items(), key=lambda x: -x[1]):
        print("  {}: {} libraries".format(batch, count))

    if args.dry_run:
        print("\n[DRY RUN] Skipping file output")
        return

    # Write output
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["library_id", "flowcell_ids", "instrument", "derived_batch", "source"]
    with open(str(args.output), "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for lib_id in sorted(results):
            writer.writerow(dict(library_id=lib_id, **results[lib_id]))

    print("\nWrote {} rows to {}".format(len(results), args.output))


if __name__ == "__main__":
    main()
