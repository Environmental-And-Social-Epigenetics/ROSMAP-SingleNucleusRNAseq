#!/usr/bin/env python3
"""
Generate side-by-side UMAPs of v1 (ORA) vs v2 (specificity) cell-type labels.

Reads obs/cell_type and obsm/X_umap directly from each annotated h5ad via
h5py (avoids the ~30 min full anndata load of the 88 GB object). Subsamples
deterministically for legible rendering.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
sys.stdout.reconfigure(line_buffering=True)

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

CELLTYPE_ORDER = [
    "Ast", "Endo", "Ex-L2/3", "Ex-L4", "Ex-L4/5", "Ex-L5", "Ex-L5/6",
    "Ex-L5/6-CC", "Ex-NRGN", "In-PV (Basket)", "In-PV (Chandelier)",
    "In-Rosehip", "In-SST", "In-VIP", "Mic", "OPC", "Oli",
]
PAL = plt.get_cmap("tab20")(np.linspace(0, 1, 20))
COLORS = {ct: PAL[i % 20] for i, ct in enumerate(CELLTYPE_ORDER)}


def read_obs_celltype_and_umap(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with h5py.File(path, "r") as f:
        ct_grp = f["obs/cell_type"]
        codes = ct_grp["codes"][:]
        raw_cats = ct_grp["categories"][:]
        cats = [c.decode("utf-8") if isinstance(c, bytes) else str(c) for c in raw_cats]
        cell_type = np.array(cats)[codes]
        umap = f["obsm/X_umap"][:]
    return cell_type, umap


def scatter(ax, umap: np.ndarray, ct: np.ndarray, title: str, n_sub: int, seed: int) -> None:
    rng = np.random.default_rng(seed)
    n_total = umap.shape[0]
    keep = (rng.choice(n_total, size=n_sub, replace=False)
            if n_sub < n_total else np.arange(n_total))
    sub_ct = ct[keep]
    sub_umap = umap[keep]
    perm = rng.permutation(len(keep))
    sub_ct = sub_ct[perm]; sub_umap = sub_umap[perm]
    colors = np.array([COLORS.get(c, (0.5, 0.5, 0.5, 1.0)) for c in sub_ct])
    ax.scatter(sub_umap[:, 0], sub_umap[:, 1], s=1.0, c=colors,
               linewidths=0, rasterized=True, alpha=0.55)
    ax.set_title(f"{title}\n(showing {len(keep):,} of {n_total:,} cells)", fontsize=11)
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--v1", type=Path, required=True)
    p.add_argument("--v2", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--n-sub", type=int, default=300_000)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    print(f"reading v1: {args.v1}", flush=True)
    ct1, umap1 = read_obs_celltype_and_umap(args.v1)
    print(f"  {len(ct1):,} cells; categories: {sorted(set(ct1))}", flush=True)
    print(f"reading v2: {args.v2}", flush=True)
    ct2, umap2 = read_obs_celltype_and_umap(args.v2)
    print(f"  {len(ct2):,} cells; categories: {sorted(set(ct2))}", flush=True)
    if umap1.shape != umap2.shape or not np.allclose(umap1, umap2, equal_nan=True):
        diff = float(np.nanmax(np.abs(umap1 - umap2))) if umap1.shape == umap2.shape else float("nan")
        print(f"  WARN: UMAP basis differs between v1/v2 (max abs diff={diff:.4e})", flush=True)
    else:
        print("  UMAP basis matches between v1 and v2.", flush=True)

    # Side-by-side panel.
    fig, axes = plt.subplots(1, 2, figsize=(20, 10))
    scatter(axes[0], umap1, ct1, "v1 — ORA top-1 per cluster (Mohammadi)",
            args.n_sub, args.seed)
    scatter(axes[1], umap2, ct2, "v2 — specificity-trimmed AUCell + cluster majority",
            args.n_sub, args.seed)
    handles = [mpatches.Patch(color=COLORS[c], label=c) for c in CELLTYPE_ORDER]
    fig.legend(handles=handles, loc="center right", bbox_to_anchor=(1.0, 0.5),
               fontsize=9, frameon=False, title="cell_type", title_fontsize=10)
    fig.suptitle("UMAP — v1 (ORA) vs v2 (specificity) cell-type labels", fontsize=13)
    fig.tight_layout(rect=[0, 0, 0.86, 0.96])
    fig.savefig(args.out_dir / "umap_v1_vs_v2.pdf", bbox_inches="tight")
    fig.savefig(args.out_dir / "umap_v1_vs_v2.png", bbox_inches="tight", dpi=180)
    plt.close(fig)
    print(f"wrote {args.out_dir/'umap_v1_vs_v2.{pdf,png}'}", flush=True)

    for tag, ct, umap, title in [
        ("v1", ct1, umap1, "v1 — ORA top-1 per cluster (Mohammadi)"),
        ("v2", ct2, umap2, "v2 — specificity-trimmed AUCell + cluster majority"),
    ]:
        fig, ax = plt.subplots(figsize=(12, 11))
        scatter(ax, umap, ct, title, args.n_sub, args.seed)
        handles = [mpatches.Patch(color=COLORS[c], label=c) for c in CELLTYPE_ORDER]
        ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1.01, 0.5),
                  fontsize=9, frameon=False, title="cell_type", title_fontsize=10)
        fig.tight_layout()
        fig.savefig(args.out_dir / f"umap_{tag}.pdf", bbox_inches="tight")
        fig.savefig(args.out_dir / f"umap_{tag}.png", bbox_inches="tight", dpi=180)
        plt.close(fig)
    print("Done.", flush=True)


if __name__ == "__main__":
    main()
