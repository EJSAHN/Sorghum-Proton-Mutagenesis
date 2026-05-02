# -*- coding: utf-8 -*-
"""
80_te_enrichment.py

TE enrichment analysis (Mask-based logic, tables-only).

Methodology:
- Parse repeatmasked GFF3 -> TE intervals per chromosome (chr1-10).
- Build boolean TE mask per chromosome (length = max(TE_end, max_event_pos) + 1000).
- TE fraction per chromosome = sum(mask)/len(mask), then averaged across chromosomes (unweighted mean).
- Mark each induced mutation as in_TE by mask lookup.
- Expected TE mutations = avg_te_fraction * total_mutations per radiation type.
- Binomial test (two-sided).

Inputs:
  - outputs_radclock/tables/events_long.csv.gz
  - inputs/SbicolorRio_468_v2.1.repeatmasked_assembly_v2.0.gff3.gz

Output:
  - outputs_functional_shielding/tables/te_mobilization_stats.csv
"""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
from scipy.stats import binomtest

from pipeline_utils import find_first_existing, resolve_project_paths


CHR10 = [str(i) for i in range(1, 11)]

# GenBank CM contigs -> chromosome numbers (with/without .1)
CM_MAP = {
    "CM027680": "1", "CM027681": "2", "CM027682": "3", "CM027683": "4",
    "CM027684": "5", "CM027685": "6", "CM027686": "7", "CM027687": "8",
    "CM027688": "9", "CM027689": "10",
    "CM027680.1": "1", "CM027681.1": "2", "CM027682.1": "3", "CM027683.1": "4",
    "CM027684.1": "5", "CM027685.1": "6", "CM027686.1": "7", "CM027687.1": "8",
    "CM027688.1": "9", "CM027689.1": "10",
}

# Phytozome-style
PHYTO_MAP = {
    "Chr01": "1", "Chr02": "2", "Chr03": "3", "Chr04": "4", "Chr05": "5",
    "Chr06": "6", "Chr07": "7", "Chr08": "8", "Chr09": "9", "Chr10": "10",
    "chr01": "1", "chr02": "2", "chr03": "3", "chr04": "4", "chr05": "5",
    "chr06": "6", "chr07": "7", "chr08": "8", "chr09": "9", "chr10": "10",
}

def normalize_chrom(name: object) -> str:
    """
    Normalize seqid/contig to '1'..'10' when possible.

    Handles:
      - CM027680(.1) with/without suffix
      - Chr01/chr01
      - chr1/Chr1/Chromosome 1/1
      - chromosome_9 / chromosome-9 / chromosome 9  -> 9

    IMPORTANT: do NOT extract arbitrary large digits (e.g., 27680) as chromosome.
    """
    s = str(name).strip()
    if not s:
        return ""

    # keep first token if seqid has extra suffix separators
    s0 = re.split(r"[ \t|;]", s)[0]
    s_base = s0.split(".")[0]

    if s0 in CM_MAP:
        return CM_MAP[s0]
    if s_base in CM_MAP:
        return CM_MAP[s_base]
    if s0 in PHYTO_MAP:
        return PHYTO_MAP[s0]

    # normalize patterns like chromosome_9 / chromosome-9
    s2 = s0.lower().replace("chromosome", "").replace("chr", "").strip()
    s2 = s2.replace("_", " ").replace("-", " ").strip()

    m = re.match(r"^0*([1-9]|10)\b", s2)
    if m:
        return str(int(m.group(1)))

    return s0


def format_p_value(p: float) -> str:
    if p < 1e-100:
        return "<1e-100"
    return f"{p:.3e}"


def parse_te_gff(gff_path: Path, debug: bool = True) -> Dict[str, List[Tuple[int, int]]]:
    """
    Return dict chrom -> list of (start,end) TE intervals for chrom in 1..10.
    We do not depend on feature type; repeatmasked gff already encodes repeats.
    """
    opener = gzip.open if str(gff_path).endswith(".gz") else open

    # --- DEBUG: show seqid samples ---
    if debug:
        raw_seqids = []
        with opener(gff_path, "rt", encoding="utf-8", errors="replace") as fdbg:
            for line in fdbg:
                if not line or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 1:
                    continue
                raw_seqids.append(parts[0])
                if len(raw_seqids) >= 120:
                    break
        uniq = []
        for x in raw_seqids:
            if x not in uniq:
                uniq.append(x)
            if len(uniq) >= 20:
                break
        print("[DEBUG] repeatmasked seqid samples (raw -> normalized):")
        for x in uniq:
            print("   ", x, "->", normalize_chrom(x))

    intervals: Dict[str, List[Tuple[int, int]]] = {c: [] for c in CHR10}

    with opener(gff_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            seqid = parts[0]
            start = parts[3]
            end = parts[4]

            chrom = normalize_chrom(seqid)
            if chrom not in CHR10:
                continue

            try:
                s = int(start)
                e = int(end)
            except ValueError:
                continue
            if s > e:
                s, e = e, s

            intervals[chrom].append((s, e))

    intervals = {c: v for c, v in intervals.items() if len(v) > 0}
    return intervals


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", default=None)
    args = ap.parse_args()

    paths = resolve_project_paths(args.project_root)

    events_path = paths.outputs_radclock / "tables" / "events_long.csv.gz"
    if not events_path.exists():
        raise FileNotFoundError(f"events_long not found: {events_path}")

    repeat_gff = find_first_existing(
        [paths.inputs, paths.root],
        [
            "SbicolorRio_468_v2.1.repeatmasked_assembly_v2.0.gff3.gz",
            "SbicolorRio_468_v2.1.repeatmasked_assembly_v2.0.gff3",
            "repeatmasked.gff3.gz",
            "repeats.gff3.gz",
        ],
    )

    print(f"[INFO] Loading events: {events_path}")
    events = pd.read_csv(events_path, compression="gzip")
    if not {"rad_type", "contig", "pos"}.issubset(set(events.columns)):
        raise ValueError("events_long missing required columns: rad_type, contig, pos")

    events["contig_norm"] = events["contig"].apply(normalize_chrom)
    events["pos"] = pd.to_numeric(events["pos"], errors="coerce")
    events = events.dropna(subset=["pos"]).copy()
    events["pos"] = events["pos"].astype(int)

    # keep chr1-10 only
    events = events[events["contig_norm"].isin(CHR10)].copy()
    if events.empty:
        raise RuntimeError("No events remain after filtering to chr1-10. Check contig normalization.")

    print(f"[INFO] Parsing repeatmasked GFF: {repeat_gff}")
    te_intervals = parse_te_gff(Path(repeat_gff), debug=False)
    if not te_intervals:
        raise RuntimeError("Parsed 0 TE intervals on chr1-10. Seqid mapping failing (see DEBUG lines above).")

    print(f"[INFO] TE chromosomes found: {sorted(te_intervals.keys())}")

    # Build per-chrom masks and TE fractions; mark TE hits
    events["in_TE"] = False
    te_fraction_dict: Dict[str, float] = {}
    total_te_hits = 0

    for chrom in sorted(te_intervals.keys(), key=lambda x: int(x)):
        iv = te_intervals[chrom]
        sub_ev = events[events["contig_norm"] == chrom]

        max_te_end = max(e for _s, e in iv)
        max_ev_pos = int(sub_ev["pos"].max()) if not sub_ev.empty else 0
        L = max(max_te_end, max_ev_pos) + 1000
        if L <= 1000:
            continue

        mask = np.zeros(L + 1, dtype=bool)
        for s, e in iv:
            if s < 1:
                s = 1
            if e >= len(mask):
                e = len(mask) - 1
            if e >= s:
                mask[s:e+1] = True

        te_len = int(mask.sum())
        te_fraction_dict[chrom] = te_len / float(len(mask))

        if not sub_ev.empty:
            idx = sub_ev.index.to_numpy()
            pos = sub_ev["pos"].to_numpy()
            valid = pos < len(mask)
            hits = mask[pos[valid]]
            events.loc[idx[valid], "in_TE"] = hits
            total_te_hits += int(hits.sum())

    if not te_fraction_dict:
        raise RuntimeError("No TE fractions computed after mask building. Check TE intervals/mask length.")
    avg_te_frac = float(np.mean(list(te_fraction_dict.values())))
    print(f"[INFO] Avg TE fraction (unweighted mean across chroms) = {avg_te_frac:.4f}")
    print(f"[INFO] Total induced mutations inside TE = {total_te_hits:,}")

    if total_te_hits == 0:
        raise RuntimeError("Total TE hits is 0. Likely coordinate/seqid mismatch in repeatmasked GFF vs events.")

    # Stats per radiation type
    rows = []
    for rad, sub in events.groupby("rad_type"):
        total = int(len(sub))
        obs = int(sub["in_TE"].sum())
        exp = float(total * avg_te_frac)
        enrich = (obs / exp) if exp > 0 else float("nan")
        pval = binomtest(k=obs, n=total, p=avg_te_frac, alternative="two-sided").pvalue

        rows.append({
            "Radiation": str(rad),
            "Observed_TE": obs,
            "Expected_TE": exp,
            "TE_Enrichment": float(enrich),
            "P_value": format_p_value(float(pval)),
        })

    out_df = pd.DataFrame(rows).sort_values("Radiation")
    out_path = paths.outputs_functional / "tables" / "te_mobilization_stats.csv"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_path, index=False)
    print(f"[OK] Saved: {out_path}")


if __name__ == "__main__":
    main()
