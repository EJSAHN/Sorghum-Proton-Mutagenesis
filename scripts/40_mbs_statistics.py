# -*- coding: utf-8 -*-
"""
40_mbs_statistics.py

Compute Multi-Base Substitution (MBS) statistics.
Legacy definition: count adjacent SNVs (distance == 1 bp) within each sample+chromosome.

Outputs:
  outputs_radclock/tables/mbs_statistics.csv
"""
from __future__ import annotations

import argparse
import pandas as pd

from pipeline_utils import resolve_project_paths


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", default=None)
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

    # Sort and compute adjacent differences
    df = df.sort_values(["sample_id", "contig", "pos"]).reset_index(drop=True)
    df["pos_diff"] = df.groupby(["sample_id", "contig"])["pos"].diff()

    # Adjacent SNVs (distance 1)
    df["is_mbs"] = df["pos_diff"] == 1

    mbs_stats = (
        df.groupby(["sample_id", "rad_type", "dose_Gy"])
          .agg(mbs_count=("is_mbs", "sum"), total_snv=("pos", "count"))
          .reset_index()
    )
    mbs_stats["mbs_rate"] = mbs_stats["mbs_count"] / mbs_stats["total_snv"]

    out_path = paths.outputs_radclock / "tables" / "mbs_statistics.csv"
    mbs_stats.to_csv(out_path, index=False)
    print(f"[OK] Saved: {out_path} ({len(mbs_stats)} samples)")


if __name__ == "__main__":
    main()
