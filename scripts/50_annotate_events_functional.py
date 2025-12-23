# -*- coding: utf-8 -*-
"""
50_annotate_events_functional.py

Functional annotation (Promoter / CDS / Intron / Intergenic) for induced SNVs.
Array-based genomic map (legacy-proven): CDS > Intron(gene body) > Promoter > Intergenic
Promoter: [gene_start-2000, gene_start-1], clipped at >=1, regardless of strand.

Inputs:
  - outputs_radclock/tables/events_long.csv.gz
  - inputs/SbicolorRio_468_v2.1.gene.gff3.gz (or common variants)

Outputs:
  - outputs_functional_shielding/tables/events_annotated_v2.csv.gz
  - outputs_functional_shielding/tables/annotation_category_counts.csv
"""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from pipeline_utils import find_first_existing, resolve_project_paths

# CM contigs (events_long/Table S4) -> chromosome number strings
CONTIG_MAP = {
    "CM027680.1": "1", "CM027681.1": "2", "CM027682.1": "3", "CM027683.1": "4",
    "CM027684.1": "5", "CM027685.1": "6", "CM027686.1": "7", "CM027687.1": "8",
    "CM027688.1": "9", "CM027689.1": "10"
}

CHR10 = [str(i) for i in range(1, 11)]
PROMOTER_SIZE = 2000

# map codes
CODE_TO_CAT = {0: "Intergenic", 1: "Promoter", 2: "Intron", 3: "CDS"}

def normalize_chrom(name: object) -> str:
    s = str(name).strip()
    if not s:
        return ""
    if s in CONTIG_MAP:
        return CONTIG_MAP[s]
    # Chr01 / chr1 / Chromosome 1 / 1
    m = re.search(r"(\d+)", s)
    if m:
        n = int(m.group(1))
        if 1 <= n <= 10:
            return str(n)
    return s

def parse_gff(gff_path: Path) -> pd.DataFrame:
    """
    Parse gene.gff3.gz and keep only chr1-10 and only 'gene' and 'CDS'.
    """
    feats = []
    with gzip.open(gff_path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom_raw = parts[0]
            ftype = parts[2]
            try:
                start = int(parts[3])
                end = int(parts[4])
            except ValueError:
                continue

            chrom = normalize_chrom(chrom_raw)
            if chrom not in CHR10:
                continue

            # keep only gene/CDS (exactly as the proven v2 logic)
            if ftype not in ("gene", "CDS"):
                continue

            if start > end:
                start, end = end, start

            feats.append({"chrom": chrom, "type": ftype, "start": start, "end": end})

    return pd.DataFrame(feats)

def build_genomic_map(gff_df: pd.DataFrame) -> Dict[str, np.ndarray]:
    """
    Build chromosome arrays:
      0 Intergenic, 1 Promoter, 2 Intron(gene body), 3 CDS
    """
    genomic_map: Dict[str, np.ndarray] = {}
    for chrom in sorted(gff_df["chrom"].unique(), key=lambda x: int(x) if x.isdigit() else 999):
        sub = gff_df[gff_df["chrom"] == chrom]
        if sub.empty:
            continue

        max_pos = int(sub["end"].max()) + PROMOTER_SIZE + 10000
        arr = np.zeros(max_pos + 1, dtype=np.int8)

        # 1) gene body -> Intron placeholder (2), promoters -> 1 (only where arr==0)
        genes = sub[sub["type"] == "gene"]
        for _, r in genes.iterrows():
            gs = int(r["start"]); ge = int(r["end"])
            arr[gs:ge+1] = 2

            ps = max(1, gs - PROMOTER_SIZE)
            pe = gs - 1
            if pe >= ps:
                region = arr[ps:pe+1]
                arr[ps:pe+1] = np.where(region == 0, 1, region)

        # 2) CDS highest priority (3)
        cdss = sub[sub["type"] == "CDS"]
        for _, r in cdss.iterrows():
            cs = int(r["start"]); ce = int(r["end"])
            arr[cs:ce+1] = 3

        genomic_map[chrom] = arr

    return genomic_map

def annotate_events(events: pd.DataFrame, genomic_map: Dict[str, np.ndarray]) -> pd.DataFrame:
    """
    Fast-ish annotation by chromosome grouping.
    """
    # default all Intergenic
    cat = np.full(len(events), "Intergenic", dtype=object)

    # group by chrom
    for chrom, idx in events.groupby("contig_norm").groups.items():
        if chrom not in genomic_map:
            continue
        arr = genomic_map[chrom]
        pos = events.loc[idx, "pos"].to_numpy(dtype=np.int64, copy=False)
        valid = (pos >= 1) & (pos < len(arr))
        if valid.any():
            codes = arr[pos[valid]]
            cat_idx = np.array(list(idx))[valid]
            cat[cat_idx] = [CODE_TO_CAT[int(c)] for c in codes]

    out = events.copy()
    out["Category"] = cat
    return out

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", default=None)
    args = ap.parse_args()

    paths = resolve_project_paths(args.project_root)

    events_path = paths.outputs_radclock / "tables" / "events_long.csv.gz"
    if not events_path.exists():
        raise FileNotFoundError(f"events_long not found: {events_path}. Run RadClock first.")

    gff = find_first_existing(
        [paths.inputs, paths.root],
        [
            "SbicolorRio_468_v2.1.gene.gff3.gz",
            "SbicolorRio_468_v2.1.gene.gff3",
            "gene.gff3.gz",
            "genes.gff3.gz",
        ],
    )

    print(f"[INFO] Loading events: {events_path}")
    events = pd.read_csv(events_path, compression="gzip")
    # normalize contig and keep chr1-10 only
    events["contig_norm"] = events["contig"].apply(normalize_chrom)
    events = events[events["contig_norm"].isin(CHR10)].copy()
    events["pos"] = pd.to_numeric(events["pos"], errors="coerce")
    events = events.dropna(subset=["pos"]).copy()
    events["pos"] = events["pos"].astype(int)

    print(f"[INFO] Loading GFF3: {gff}")
    gff_df = parse_gff(Path(gff))
    if gff_df.empty:
        raise RuntimeError("GFF parsing produced 0 gene/CDS features on chr1-10. Check GFF feature types and seqids.")

    print("[INFO] Building genomic map arrays (chr1-10)...")
    genomic_map = build_genomic_map(gff_df)

    print("[INFO] Annotating events...")
    annotated = annotate_events(events, genomic_map)

    out_events = paths.outputs_functional / "tables" / "events_annotated_v2.csv.gz"
    out_events.parent.mkdir(parents=True, exist_ok=True)
    annotated.to_csv(out_events, index=False, compression="gzip")
    print(f"[OK] Saved: {out_events} ({len(annotated):,} events)")

    summary = (
        annotated.groupby(["rad_type", "Category"])
        .size()
        .reset_index(name="count")
        .sort_values(["rad_type", "Category"])
    )
    out_sum = paths.outputs_functional / "tables" / "annotation_category_counts.csv"
    summary.to_csv(out_sum, index=False)
    print(f"[OK] Saved: {out_sum}")

if __name__ == "__main__":
    main()
