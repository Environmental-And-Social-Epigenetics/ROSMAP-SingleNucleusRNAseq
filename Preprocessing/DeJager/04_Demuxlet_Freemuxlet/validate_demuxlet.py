"""
Validate demuxlet or freemuxlet results against ground truth annotations.

Computes Adjusted Rand Index (ARI), Normalized Mutual Information (NMI),
and generates a confusion matrix heatmap comparing demuxlet assignments
to known cell annotations.

Usage:
    python validate_demuxlet.py \
        --demux-file output/demux1.best \
        --annotation-file cell-annotation.csv \
        --library-id 191121-B6 \
        --output-dir results/

    # For freemuxlet output (gzipped):
    python validate_demuxlet.py \
        --demux-file output/freemux.clust1.samples.gz \
        --annotation-file cell-annotation.csv \
        --library-id 190409-B5-B \
        --output-dir results/ \
        --freemuxlet
"""

import argparse
import gzip
import os
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.metrics import adjusted_rand_score, confusion_matrix, normalized_mutual_info_score


def load_demux_output(filepath, is_freemuxlet=False):
    """Load demuxlet or freemuxlet output file."""
    if is_freemuxlet and filepath.endswith(".gz"):
        decompressed = filepath.rstrip(".gz")
        with gzip.open(filepath, "rb") as f_in:
            with open(decompressed, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        filepath = decompressed

    df = pd.read_csv(filepath, sep="\t")
    df.rename(columns={"BARCODE": "barcode"}, inplace=True)

    if not is_freemuxlet:
        # Demuxlet: extract first patient from comma-separated BEST.GUESS
        df["BEST.GUESS"] = df["BEST.GUESS"].str.split(",").str[0]

    return df


def load_annotations(filepath, library_id):
    """Load and filter cell annotations for a specific library."""
    df = pd.read_csv(filepath)
    prefix = f"{library_id}_"
    df = df[df["barcode"].str.startswith(prefix)].copy()
    # Remove library prefix from barcodes to match demuxlet output
    df["barcode"] = df["barcode"].str[len(prefix):]
    return df


def validate(demux_df, annotation_df, output_dir, library_id, is_freemuxlet=False):
    """Compute metrics and generate confusion matrix heatmap."""
    merged = pd.merge(demux_df, annotation_df, on="barcode")
    merged = merged.dropna(subset=["BEST.GUESS", "individualID"])

    if merged.empty:
        print("ERROR: No overlapping barcodes found between demux output and annotations.")
        return

    # Compute metrics
    ari = adjusted_rand_score(merged["BEST.GUESS"], merged["individualID"])
    nmi = normalized_mutual_info_score(merged["BEST.GUESS"], merged["individualID"])
    mismatches = (merged["BEST.GUESS"] != merged["individualID"]).sum()

    print(f"Library: {library_id}")
    print(f"Overlapping cells: {len(merged)}")
    print(f"Adjusted Rand Index: {ari:.4f}")
    print(f"Normalized Mutual Information: {nmi:.4f}")
    print(f"Mismatched assignments: {mismatches} / {len(merged)}")

    # Confusion matrix heatmap
    labels_true = merged["individualID"]
    labels_pred = merged["BEST.GUESS"]
    all_labels = sorted(set(labels_true) | set(labels_pred))

    cm = confusion_matrix(labels_true, labels_pred, labels=all_labels)
    plt.figure(figsize=(max(8, len(all_labels)), max(6, len(all_labels) * 0.8)))
    sns.heatmap(cm, annot=True, fmt="d", xticklabels=all_labels, yticklabels=all_labels)
    plt.xlabel("Demuxlet Assignment" if not is_freemuxlet else "Freemuxlet Cluster")
    plt.ylabel("Ground Truth (Annotation)")
    plt.title(f"Confusion Matrix - {library_id}")
    plt.tight_layout()

    os.makedirs(output_dir, exist_ok=True)
    method = "freemuxlet" if is_freemuxlet else "demuxlet"
    plot_path = os.path.join(output_dir, f"confusion_matrix_{method}_{library_id}.png")
    plt.savefig(plot_path, dpi=150)
    plt.close()
    print(f"Heatmap saved to: {plot_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Validate demuxlet/freemuxlet results against ground truth annotations."
    )
    parser.add_argument("--demux-file", required=True,
                        help="Path to demux1.best or freemux.clust1.samples.gz")
    parser.add_argument("--annotation-file", required=True,
                        help="Path to cell-annotation.csv with ground truth")
    parser.add_argument("--library-id", required=True,
                        help="Library ID (e.g., 191121-B6)")
    parser.add_argument("--output-dir", default=".",
                        help="Directory for output plots (default: current directory)")
    parser.add_argument("--freemuxlet", action="store_true",
                        help="Input is freemuxlet output (not demuxlet)")
    args = parser.parse_args()

    demux_df = load_demux_output(args.demux_file, is_freemuxlet=args.freemuxlet)
    annotation_df = load_annotations(args.annotation_file, args.library_id)

    if annotation_df.empty:
        print(f"ERROR: No annotations found for library '{args.library_id}'")
        return

    validate(demux_df, annotation_df, args.output_dir, args.library_id,
             is_freemuxlet=args.freemuxlet)


if __name__ == "__main__":
    main()
