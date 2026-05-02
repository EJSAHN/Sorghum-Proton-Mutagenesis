# -*- coding: utf-8 -*-
"""
20_spatial_stats_1mb.py  (Paper/RAW reproduction: histogram range=(0,max_pos))

This reproduces the original RAW logic used in sorghum_spatial_stats.py:
- chr1-10 only
- per chromosome:
    n_bins = int(max_pos / 1Mb) + 1
    hist, _ = np.histogram(pos, bins=n_bins, range=(0, max_pos))
    densities = hist / n_samples
- concatenate densities across chromosomes
- metrics:
    Mean_Density
    Max_Spike_Density
    CV (Dispersion) = std/mean
    Gini_Coefficient
    Hotspots_Count = bins where density > mean + 3*std

Inputs:
  outputs_radclock/tables/events_long.csv.gz

Outputs:
  outputs_radclock/tables/spatial_statistics.csv
  outputs_radclock/tables/ks_test_results.csv
"""

from __future__ import annotations

import argparse
import re
from typing import Dict, List

import numpy as np
import pandas as pd
from scipy.stats import ks_2samp

from pipeline_utils import resolve_project_paths


BIN_SIZE = 1_000_000
CHR10 = [str(i) for i in range(1, 11)]

_CONTIG_MAP = {
    "CM027680.1": "1", "CM027681.1": "2", "CM027682.1": "3", "CM027683.1": "4",
    "CM027684.1": "5", "CM027685.1": "6", "CM027686.1": "7", "CM027687.1": "8",
    "CM027688.1": "9", "CM027689.1": "10",
    "Chr01": "1", "Chr02": "2", "Chr03": "3", "Chr04": "4", "Chr05": "5",
    "Chr06": "6", "Chr07": "7", "Chr08": "8", "Chr09": "9", "Chr10": "10",
    "chr01": "1", "chr02": "2", "chr03": "3", "chr04": "4", "chr05": "5",
    "chr06": "6", "chr07": "7", "chr08": "8", "chr09": "9", "chr10": "10",
    "1": "1", "2": "2", "3": "3", "4": "4", "5": "5",
    "6": "6", "7": "7", "8": "8", "9": "9", "10": "10",
}

def norm_chrom(x: object) -> str:
    s = str(x).strip()
    if not s:
        return ""
    if s in _CONTIG_MAP:
        return _CONTIG_MAP[s]
    s0 = s.lower().strip().replace("chromosome", "").replace("chr", "").strip()
    s0 = s0.replace("_", " ").replace("-", " ").strip()
    m = re.match(r"^0*([1-9]|10)\b", s0)
    if m:
        return str(int(m.group(1)))
    return s


def gini_coefficient(x: np.ndarray) -> float:
    if x.size == 0:
        return float("nan")
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return float("nan")
    if np.any(x < 0):
        x = x - np.min(x)
    s = float(x.sum())
    if s == 0.0:
        return 0.0
    x_sorted = np.sort(x)
    n = x_sorted.size
    idx = np.arange(1, n + 1, dtype=float)
    return float(np.sum((2 * idx - n - 1) * x_sorted) / (n * s))


def build_raw_hist_densities(rad_df: pd.DataFrame, n_samples: int) -> np.ndarray:
    """
    RAW logic:
      n_bins = int(max_pos/bin) + 1
      histogram range = (0, max_pos)  <-- this is the key difference
    """
    dens_all: List[float] = []

    for chrom in CHR10:
        chrom_df = rad_df[rad_df["chrom_norm"] == chrom]
        if chrom_df.empty:
            continue

        max_pos = int(chrom_df["pos"].max())
        if max_pos <= 0:
            continue

        n_bins = int(max_pos / BIN_SIZE) + 1
        # IMPORTANT: range=(0, max_pos) (RAW behavior)
        hist, _ = np.histogram(chrom_df["pos"].values, bins=n_bins, range=(0, max_pos))

        dens = hist.astype(float) / float(max(1, n_samples))
        dens_all.extend(dens.tolist())

    return np.array(dens_all, dtype=float)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", default=None)
    args = ap.parse_args()

    paths = resolve_project_paths(args.project_root)

    events_path = paths.outputs_radclock / "tables" / "events_long.csv.gz"
    if not events_path.exists():
        raise FileNotFoundError(f"events_long not found: {events_path}")

    print(f"[INFO] Loading events: {events_path}")
    df = pd.read_csv(events_path, compression="gzip")

    df["rad_type"] = df["rad_type"].astype(str).str.lower().str.strip()
    df["chrom_norm"] = df["contig"].apply(norm_chrom).astype(str)
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    df = df.dropna(subset=["pos"]).copy()
    df["pos"] = df["pos"].astype(int)

    before = len(df)
    df = df[df["chrom_norm"].isin(CHR10)].copy()
    after = len(df)
    print(f"[INFO] Filter chr1-10: {before:,} -> {after:,}")
    if df.empty:
        raise RuntimeError("No events left after chr1-10 filtering.")

    out_rows = []
    densities_by_rad: Dict[str, np.ndarray] = {}

    # Paper row order: Gamma then Proton
    for rad in ["gamma", "proton"]:
        sub = df[df["rad_type"] == rad].copy()
        if sub.empty:
            continue

        n_samples = int(sub["sample_id"].nunique())
        print(f"[INFO] [{rad}] samples={n_samples:,} events={len(sub):,}")

        dens = build_raw_hist_densities(sub, n_samples=n_samples)
        if dens.size == 0:
            raise RuntimeError(f"[{rad}] No densities computed.")

        densities_by_rad[rad] = dens

        mean_d = float(np.mean(dens))
        std_d = float(np.std(dens, ddof=0))
        max_d = float(np.max(dens))
        cv = float(std_d / mean_d) if mean_d > 0 else float("nan")
        gini = float(gini_coefficient(dens))
        hotspots = int(np.sum(dens > (mean_d + 3.0 * std_d)))

        out_rows.append({
            "Radiation": "Gamma" if rad == "gamma" else "Proton",
            "Mean_Density": mean_d,
            "Max_Spike_Density": max_d,
            "CV (Dispersion)": cv,
            "Gini_Coefficient": gini,
            "Hotspots_Count": hotspots,
        })

        print(f"[INFO] [{rad}] Mean={mean_d:.6f} Max={max_d:.6f} CV={cv:.6f} Gini={gini:.6f} Hotspots={hotspots}")

    out_df = pd.DataFrame(out_rows, columns=[
        "Radiation",
        "Mean_Density",
        "Max_Spike_Density",
        "CV (Dispersion)",
        "Gini_Coefficient",
        "Hotspots_Count",
    ])

    out_dir = paths.outputs_radclock / "tables"
    out_dir.mkdir(parents=True, exist_ok=True)

    out_csv = out_dir / "spatial_statistics.csv"
    out_df.to_csv(out_csv, index=False, float_format="%.6f")
    print(f"[OK] Saved: {out_csv}")

    # KS test (optional, keeps your pipeline stable)
    ks_rows = []
    if "gamma" in densities_by_rad and "proton" in densities_by_rad:
        stat, pval = ks_2samp(densities_by_rad["gamma"], densities_by_rad["proton"])
        ks_rows.append({"Test": "KS_2samp(density): gamma vs proton", "KS_stat": float(stat), "P_value": float(pval)})

    ks_df = pd.DataFrame(ks_rows)
    ks_csv = out_dir / "ks_test_results.csv"
    ks_df.to_csv(ks_csv, index=False, float_format="%.6g")
    print(f"[OK] Saved: {ks_csv}")


if __name__ == "__main__":
    main()
