#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Cluster sensitivity analysis: recurrent-locus filtering.

Purpose
-------
Test whether short-range clustering (within +/- bp) persists after removing loci that recur
across multiple samples (potential locus-level artifacts / high-callability sites).

Input
-----
--events   Path to events_long.csv(.gz) with columns:
           sample_id, rad_type, contig, pos
--outdir   Output directory

Key parameters
--------------
--cluster_bp   Clustering distance (default: 10)
--min_k_list   Comma-separated thresholds for recurrence filtering.
               A locus observed in >= k samples will be removed.
               (default: "2,3,5,10")

Outputs
-------
tables/
  recurrent_locus_counts.csv
  cluster_fraction_sensitivity.csv

figures/
  cluster_fraction_sensitivity.pdf
  cluster_fraction_sensitivity.png (300 dpi)

Notes
-----
- No hard-coded paths.
- Figures are saved as PDF + PNG (300 dpi).
- No numeric labels are embedded in the figure (values are in CSV).
"""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# -----------------------------
# Helpers
# -----------------------------

def normalize_chrom(x) -> str:
    s = str(x).strip()
    # Prefer mapping CM027680.1 -> 1..10 if present
    cm_map = {
        "CM027680.1": "1", "CM027681.1": "2", "CM027682.1": "3", "CM027683.1": "4",
        "CM027684.1": "5", "CM027685.1": "6", "CM027686.1": "7", "CM027687.1": "8",
        "CM027688.1": "9", "CM027689.1": "10"
    }
    if s in cm_map:
        return cm_map[s]
    m = re.search(r"(\d+)", s)
    if m:
        v = int(m.group(1))
        if 1 <= v <= 10:
            return str(v)
    return s


def read_events(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Events file not found: {path}")
    if str(path).endswith(".gz"):
        df = pd.read_csv(path, compression="gzip")
    else:
        df = pd.read_csv(path)

    needed = {"sample_id", "rad_type", "contig", "pos"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"Events missing required columns: {sorted(missing)}")

    df = df.copy()
    df["sample_id"] = pd.to_numeric(df["sample_id"], errors="coerce").astype(int)
    df["rad_type"] = df["rad_type"].astype(str).str.lower()
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce").astype(int)
    df["chrom"] = df["contig"].apply(normalize_chrom)
    # keep main chromosomes only if recognizable
    df = df[df["chrom"].isin([str(i) for i in range(1, 11)])].copy()
    return df


def compute_cluster_fraction(df: pd.DataFrame, cluster_bp: int) -> pd.DataFrame:
    ev = df.sort_values(["sample_id", "chrom", "pos"]).copy()
    ev["prev_pos"] = ev.groupby(["sample_id", "chrom"])["pos"].shift(1)
    ev["next_pos"] = ev.groupby(["sample_id", "chrom"])["pos"].shift(-1)
    ev["d_prev"] = ev["pos"] - ev["prev_pos"]
    ev["d_next"] = ev["next_pos"] - ev["pos"]
    ev["clustered"] = (ev["d_prev"] <= cluster_bp) | (ev["d_next"] <= cluster_bp)

    out = (ev.groupby(["sample_id", "rad_type"])
             .agg(total_events=("pos", "count"),
                  clustered_events=("clustered", "sum"))
             .reset_index())
    out["cluster_fraction"] = out["clustered_events"] / out["total_events"]
    return out


# -----------------------------
# Main
# -----------------------------

def main():
    ap = argparse.ArgumentParser(description="Cluster sensitivity to recurrent-locus filtering.")
    ap.add_argument("--events", required=True, type=str, help="Path to events_long.csv(.gz)")
    ap.add_argument("--outdir", required=True, type=str, help="Output directory")
    ap.add_argument("--cluster_bp", default=10, type=int, help="Clustering distance (bp)")
    ap.add_argument("--min_k_list", default="2,3,5,10", type=str,
                    help="Comma-separated recurrence thresholds (remove loci seen in >=k samples)")
    args = ap.parse_args()

    events_path = Path(args.events)
    outdir = Path(args.outdir)
    tabdir = outdir / "tables"
    figdir = outdir / "figures"
    tabdir.mkdir(parents=True, exist_ok=True)
    figdir.mkdir(parents=True, exist_ok=True)

    sns.set_style("whitegrid")

    df = read_events(events_path)

    # Define locus key at chrom:pos
    df["locus"] = df["chrom"].astype(str) + ":" + df["pos"].astype(str)

    # Recurrence: number of unique samples per locus
    locus_sample_counts = (df.groupby("locus")["sample_id"]
                             .nunique()
                             .reset_index(name="n_samples_with_locus"))
    locus_sample_counts.to_csv(tabdir / "recurrent_locus_counts.csv", index=False)

    # Prepare thresholds
    min_k: List[int] = []
    for tok in args.min_k_list.split(","):
        tok = tok.strip()
        if tok:
            min_k.append(int(tok))
    min_k = sorted(set(min_k))

    # Baseline (no filtering) = k = inf (represented as 0)
    results = []
    base_cf = compute_cluster_fraction(df, args.cluster_bp)
    base_cf["filter_k"] = 0  # 0 means no recurrent-locus filtering
    results.append(base_cf)

    # Apply recurrence filtering
    # Remove loci with recurrence >= k
    for k in min_k:
        bad = set(locus_sample_counts.loc[locus_sample_counts["n_samples_with_locus"] >= k, "locus"])
        sub = df[~df["locus"].isin(bad)].copy()
        cf = compute_cluster_fraction(sub, args.cluster_bp)
        cf["filter_k"] = k
        results.append(cf)

    sens = pd.concat(results, ignore_index=True)

    # Summarize by rad_type and filter_k
    summ = (sens.groupby(["rad_type", "filter_k"])
                 .agg(n_samples=("sample_id", "nunique"),
                      mean_cluster_fraction=("cluster_fraction", "mean"),
                      median_cluster_fraction=("cluster_fraction", "median"))
                 .reset_index())
    summ.to_csv(tabdir / "cluster_fraction_sensitivity.csv", index=False)

    # Plot: mean cluster fraction by filter_k
    # Use lineplot with points; no numeric annotation.
    fig, ax = plt.subplots(figsize=(6.8, 4.2))

    plot_df = summ.copy()
    # Make filter_k labels: 0 -> "None", else ">=k"
    plot_df["filter_label"] = plot_df["filter_k"].apply(lambda x: "None" if x == 0 else f">={x} samples")

    # Ensure order
    order = ["None"] + [f">={k} samples" for k in min_k]
    plot_df["filter_label"] = pd.Categorical(plot_df["filter_label"], categories=order, ordered=True)

    sns.lineplot(data=plot_df, x="filter_label", y="mean_cluster_fraction",
                 hue="rad_type", style="rad_type", markers=True, dashes=False, ax=ax)

    ax.set_title("Cluster fraction sensitivity to recurrent-locus filtering")
    ax.set_xlabel("Removed loci recurring in")
    ax.set_ylabel(f"Mean cluster fraction (±{args.cluster_bp} bp)")
    ax.tick_params(axis="x", rotation=25)
    fig.tight_layout()

    fig.savefig(figdir / "cluster_fraction_sensitivity.pdf", format="pdf")
    fig.savefig(figdir / "cluster_fraction_sensitivity.png", format="png", dpi=300)

    print("Done.")
    print(str(outdir))


if __name__ == "__main__":
    main()
