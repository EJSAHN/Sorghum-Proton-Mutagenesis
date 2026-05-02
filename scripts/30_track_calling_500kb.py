# -*- coding: utf-8 -*-
"""
30_track_calling_500kb.py

Detect mutation "tracks" as high-density windows in induced SNVs.
Legacy definition:
- Bin induced SNVs into 500 kb windows per chromosome.
- Determine a global track threshold from gamma samples: 95th percentile of non-zero window counts.
- For each sample: track windows are those with count >= threshold.

Outputs:
  outputs_radclock/tables/track_summary_per_sample.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd

from pipeline_utils import resolve_project_paths


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", default=None)
    ap.add_argument("--window", type=int, default=500_000, help="Window size in bp (default 500,000).")
    ap.add_argument("--gamma-percentile", type=float, default=95.0, help="Percentile for gamma-derived threshold (default 95).")
    args = ap.parse_args()

    paths = resolve_project_paths(args.project_root)
    events_path = paths.outputs_radclock / "tables" / "events_long.csv.gz"
    if not events_path.exists():
        raise FileNotFoundError(f"events_long not found: {events_path}. Run RadClock step first.")

    df = pd.read_csv(events_path, compression="gzip")
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    df = df.dropna(subset=["pos"]).copy()
    df["pos"] = df["pos"].astype(int)

    # Normalize contigs (legacy)
    df["contig"] = df["contig"].astype(str).str.replace("Chr", "", regex=False).str.replace("chr", "", regex=False)

    window = int(args.window)

    # ------------------------------------------------------------------
    # 1) Determine threshold from gamma window counts (non-zero windows only)
    # ------------------------------------------------------------------
    gamma_df = df[df["rad_type"] == "gamma"].copy()
    gamma_window_counts: list[int] = []

    for (sid, contig), sub in gamma_df.groupby(["sample_id", "contig"]):
        counts = sub.groupby(sub["pos"] // window).size()  # non-zero windows only
        gamma_window_counts.extend(counts.values.tolist())

    if len(gamma_window_counts) == 0:
        raise RuntimeError("No gamma window counts found. Check rad_type labels and events_long content.")

    threshold = float(np.percentile(gamma_window_counts, args.gamma_percentile))
    print(f"[INFO] Track threshold (gamma {args.gamma_percentile}th percentile) = {threshold:.2f} mutations/window")

    # ------------------------------------------------------------------
    # 2) Per-sample track summary
    # ------------------------------------------------------------------
    results = []
    for sid, sample_df in df.groupby("sample_id"):
        rad_type = str(sample_df["rad_type"].iloc[0])
        dose = int(sample_df["dose_Gy"].iloc[0]) if "dose_Gy" in sample_df.columns and pd.notna(sample_df["dose_Gy"].iloc[0]) else None

        total_mut = int(len(sample_df))
        track_mut = 0
        n_tracks = 0
        track_densities: list[int] = []

        for contig, chrom_df in sample_df.groupby("contig"):
            counts = chrom_df.groupby(chrom_df["pos"] // window).size()
            high = counts[counts >= threshold]
            if len(high) == 0:
                continue
            track_mut += int(high.sum())
            n_tracks += int(high.shape[0])
            track_densities.extend(high.values.tolist())

        track_fraction = (track_mut / total_mut) if total_mut > 0 else 0.0
        mean_track_density = float(np.mean(track_densities)) if track_densities else 0.0
        max_track_density = float(np.max(track_densities)) if track_densities else 0.0

        results.append({
            "sample_id": int(sid),
            "rad_type": rad_type,
            "dose_Gy": dose,
            "total_mutations": total_mut,
            "track_mutations": int(track_mut),
            "track_fraction": float(track_fraction),
            "n_tracks": int(n_tracks),
            "mean_track_density": mean_track_density,
            "max_track_density": max_track_density,
        })

    out_df = pd.DataFrame(results).sort_values(["rad_type", "dose_Gy", "sample_id"], na_position="last")
    out_path = paths.outputs_radclock / "tables" / "track_summary_per_sample.csv"
    out_df.to_csv(out_path, index=False)
    print(f"[OK] Saved: {out_path} ({len(out_df)} samples)")


if __name__ == "__main__":
    main()
