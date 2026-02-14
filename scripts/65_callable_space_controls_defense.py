#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Sorghum callable-space controls for GBS-derived induced SNV analyses.

This script provides quantitative controls addressing common concerns in reduced-representation datasets:
1) Callable-space QC (per-sample callable loci and missingness) using Table S4 genotype calls
2) Feature-specific callable denominators (CDS, promoter, gene body non-CDS, intergenic, TE overlap)
3) Callable-space-corrected enrichment of induced SNVs by feature
4) Design transparency (origin x rad_type crosstab, if origin present)
5) Cluster fraction QC (local clustering within a specified bp window)

Important implementation details for Table S4:
- Genotype calls may include IUPAC ambiguity codes (e.g., R, Y, W, M, S, K) that represent
  heterozygous/ambiguous bases. These are treated as callable.
- Missing calls include N/n/empty/NaN.
- Depth columns may exist and may be encoded as strings such as "4,3|7" or "14|14".
  This script does not require depth thresholds, but it normalizes depth column names for completeness.
- Column names may contain suffixes like "23_1" / "23_1_Depth" from Excel. These are normalized to "23" / "23_Depth"
  when safe to do so.

Outputs are written to:
  <outdir>/tables/
  <outdir>/figures/

Figures are saved as:
  PDF (vector) and PNG (300 dpi)

No numeric annotations are embedded in figures; values are provided in tables.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib.pyplot as plt
import seaborn as sns


# -----------------------------
# Chromosome normalization
# -----------------------------

CONTIG_MAP: Dict[str, str] = {
    "CM027680.1": "1",
    "CM027681.1": "2",
    "CM027682.1": "3",
    "CM027683.1": "4",
    "CM027684.1": "5",
    "CM027685.1": "6",
    "CM027686.1": "7",
    "CM027687.1": "8",
    "CM027688.1": "9",
    "CM027689.1": "10",
    "Chr01": "1",
    "Chr02": "2",
    "Chr03": "3",
    "Chr04": "4",
    "Chr05": "5",
    "Chr06": "6",
    "Chr07": "7",
    "Chr08": "8",
    "Chr09": "9",
    "Chr10": "10",
    "1": "1",
    "2": "2",
    "3": "3",
    "4": "4",
    "5": "5",
    "6": "6",
    "7": "7",
    "8": "8",
    "9": "9",
    "10": "10",
}
MAIN_CHROMS = [str(i) for i in range(1, 11)]


def normalize_chrom(name: str) -> Optional[str]:
    s = str(name).strip()
    if s in CONTIG_MAP:
        return CONTIG_MAP[s]
    m = re.search(r"(\d+)", s)
    if m:
        v = int(m.group(1))
        if 1 <= v <= 10:
            return str(v)
    return None


# -----------------------------
# Genotype call parsing
# -----------------------------

# Callable includes A/C/G/T and IUPAC ambiguity codes that represent valid base sets.
CALLABLE_IUPAC = set(list("ACGTRYSWKM"))
MISSING_IUPAC = set(["N", ""])


def is_callable_genotype(x) -> bool:
    if pd.isna(x):
        return False
    s = str(x).strip().upper()
    if s in MISSING_IUPAC:
        return False
    # Sometimes missing is encoded as lowercase 'n'
    if s == "N":
        return False
    return s in CALLABLE_IUPAC


def normalize_s4_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalize Table S4 column names:
    - Convert Excel-suffixed genotype columns like "23_1" -> "23" (if "23" not already present)
    - Convert depth columns like "23_1_Depth" -> "23_Depth"
    """
    cols = list(df.columns)
    new_cols = cols[:]

    # Build a set of existing exact names for collision checks.
    existing = set(str(c) for c in cols)

    for i, c in enumerate(cols):
        cs = str(c).strip()

        # Normalize genotype numeric-with-suffix: "23_1" -> "23"
        m = re.fullmatch(r"(\d+)_\d+", cs)
        if m:
            base = m.group(1)
            if base not in existing:
                new_cols[i] = base
                continue

        # Normalize depth numeric-with-suffix: "23_1_Depth" -> "23_Depth"
        m = re.fullmatch(r"(\d+)_\d+_Depth", cs, flags=re.IGNORECASE)
        if m:
            base = m.group(1)
            target = f"{base}_Depth"
            if target not in existing:
                new_cols[i] = target
                continue

        # Normalize "23 Depth" / "23-Depth" -> "23_Depth"
        m = re.fullmatch(r"(\d+)[ \-]Depth", cs, flags=re.IGNORECASE)
        if m:
            base = m.group(1)
            target = f"{base}_Depth"
            if target not in existing:
                new_cols[i] = target
                continue

    df = df.copy()
    df.columns = new_cols
    return df


def parse_depth_total(x) -> float:
    """
    Parse depth encoded as strings such as:
    - "14|14" -> 14
    - "4,3|7" -> 7 (treat right side of '|' as total depth if present)
    - "0|1" -> 1
    - numeric -> numeric
    Returns NaN on failure.
    """
    if pd.isna(x):
        return np.nan
    if isinstance(x, (int, float, np.integer, np.floating)):
        return float(x)
    s = str(x).strip()
    if s == "" or s.lower() == "nan":
        return np.nan

    # Prefer right-hand side after '|'
    if "|" in s:
        rhs = s.split("|")[-1].strip()
        try:
            return float(rhs)
        except Exception:
            pass

    # Otherwise sum comma-separated parts if numeric
    parts = re.split(r"[,\s]+", s)
    nums = []
    for p in parts:
        p = p.strip()
        if p == "":
            continue
        try:
            nums.append(float(p))
        except Exception:
            continue
    if nums:
        return float(np.sum(nums))
    return np.nan


# -----------------------------
# GFF interval utilities
# -----------------------------

Interval = Tuple[int, int]


def merge_intervals(intervals: List[Interval]) -> List[Interval]:
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: (x[0], x[1]))
    out: List[List[int]] = [[intervals[0][0], intervals[0][1]]]
    for s, e in intervals[1:]:
        if s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(a, b) for a, b in out]


def subtract_intervals(A: List[Interval], B: List[Interval]) -> List[Interval]:
    if not A:
        return []
    if not B:
        return A.copy()
    out: List[Interval] = []
    j = 0
    for a_start, a_end in A:
        cur = a_start
        while j < len(B) and B[j][1] < cur:
            j += 1
        k = j
        while k < len(B) and B[k][0] <= a_end:
            b_start, b_end = B[k]
            if b_start > cur:
                out.append((cur, min(a_end, b_start - 1)))
            cur = max(cur, b_end + 1)
            if cur > a_end:
                break
            k += 1
        if cur <= a_end:
            out.append((cur, a_end))
    return out


def build_arrays(intervals: List[Interval]) -> Tuple[np.ndarray, np.ndarray]:
    if not intervals:
        return np.array([], dtype=np.int64), np.array([], dtype=np.int64)
    starts = np.array([s for s, _ in intervals], dtype=np.int64)
    ends = np.array([e for _, e in intervals], dtype=np.int64)
    return starts, ends


def membership_mask(positions: np.ndarray, starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
    if starts.size == 0:
        return np.zeros(len(positions), dtype=bool)
    idx = np.searchsorted(starts, positions, side="right") - 1
    valid = idx >= 0
    idx_clip = idx.copy()
    idx_clip[~valid] = 0
    return valid & (positions <= ends[idx_clip])


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def parse_gene_gff(gff_path: Path, promoter_upstream_bp: int):
    gene: Dict[str, List[Interval]] = {}
    cds: Dict[str, List[Interval]] = {}
    promoter: Dict[str, List[Interval]] = {}

    with open_text(gff_path) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, _, ftype, start, end, _, strand, _, _ = parts
            chrom = normalize_chrom(seqid)
            if chrom is None or chrom not in MAIN_CHROMS:
                continue
            try:
                s = int(start)
                e = int(end)
            except ValueError:
                continue
            if s > e:
                s, e = e, s

            ftype_l = ftype.strip().lower()
            if ftype_l == "gene":
                gene.setdefault(chrom, []).append((s, e))
                if strand == "+":
                    ps = max(1, s - promoter_upstream_bp)
                    pe = s - 1
                elif strand == "-":
                    ps = e + 1
                    pe = e + promoter_upstream_bp
                else:
                    continue
                if pe >= ps:
                    promoter.setdefault(chrom, []).append((ps, pe))
            elif ftype_l == "cds":
                cds.setdefault(chrom, []).append((s, e))

    for d in (gene, cds, promoter):
        for k in list(d.keys()):
            d[k] = merge_intervals(d[k])

    for chrom in list(promoter.keys()):
        prom_only = subtract_intervals(promoter.get(chrom, []), gene.get(chrom, []))
        promoter[chrom] = merge_intervals(prom_only)

    return cds, gene, promoter


def parse_te_gff(repeat_gff_path: Path) -> Dict[str, List[Interval]]:
    te: Dict[str, List[Interval]] = {}
    with open_text(repeat_gff_path) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, _, _, start, end, *_ = parts
            chrom = normalize_chrom(seqid)
            if chrom is None or chrom not in MAIN_CHROMS:
                continue
            try:
                s = int(start)
                e = int(end)
            except ValueError:
                continue
            if s > e:
                s, e = e, s
            te.setdefault(chrom, []).append((s, e))

    for k in list(te.keys()):
        te[k] = merge_intervals(te[k])
    return te


# -----------------------------
# Loading inputs
# -----------------------------

def read_meta(meta_path: Path) -> pd.DataFrame:
    if not meta_path.exists():
        raise FileNotFoundError(f"Meta file not found: {meta_path}")

    if meta_path.suffix.lower() in [".xlsx", ".xls"]:
        df = pd.read_excel(meta_path, sheet_name="S1_Sample_Metadata")
    else:
        df = pd.read_csv(meta_path)

    required = ["sample_id", "rad_type", "dose_Gy", "parent_id"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Meta file missing required columns: {missing}")

    df = df.copy()
    df["sample_id"] = pd.to_numeric(df["sample_id"], errors="coerce").astype(int)
    df["rad_type"] = df["rad_type"].astype(str).str.lower()
    df["parent_id"] = pd.to_numeric(df["parent_id"], errors="coerce").astype("Int64")
    if "origin" in df.columns:
        df["origin"] = df["origin"].astype(str)
    return df


def read_events(events_path: Path) -> pd.DataFrame:
    if not events_path.exists():
        raise FileNotFoundError(f"Events file not found: {events_path}")

    if str(events_path).endswith(".gz"):
        ev = pd.read_csv(events_path, compression="gzip")
    else:
        ev = pd.read_csv(events_path)

    required = {"sample_id", "rad_type", "contig", "pos"}
    missing = required - set(ev.columns)
    if missing:
        raise ValueError(f"Events missing required columns: {sorted(missing)}")

    ev = ev.copy()
    ev["sample_id"] = pd.to_numeric(ev["sample_id"], errors="coerce").astype(int)
    ev["rad_type"] = ev["rad_type"].astype(str).str.lower()
    ev["contig_norm"] = ev["contig"].apply(lambda x: normalize_chrom(x) or "NA")
    ev = ev[ev["contig_norm"].isin(MAIN_CHROMS)].copy()
    ev["pos"] = pd.to_numeric(ev["pos"], errors="coerce").astype(int)
    return ev


def read_table_s4(s4_path: Path, sample_ids_needed: List[int], header_row: int) -> pd.DataFrame:
    if not s4_path.exists():
        raise FileNotFoundError(f"Table S4 not found: {s4_path}")

    df = pd.read_excel(s4_path, header=header_row, engine="openpyxl")
    df = df.dropna(how="all")
    df = normalize_s4_columns(df)

    rename = {}
    for c in df.columns:
        lc = str(c).strip().lower()
        if lc == "id":
            rename[c] = "contig"
        elif lc == "pos":
            rename[c] = "pos"
        elif lc == "ref":
            rename[c] = "ref"
    df = df.rename(columns=rename)

    if "contig" not in df.columns or "pos" not in df.columns:
        raise ValueError("Table S4 must contain columns 'id' (or 'contig') and 'pos'.")

    wanted = set(int(x) for x in sample_ids_needed)

    geno_cols = []
    depth_cols = []

    for c in df.columns:
        cs = str(c).strip()
        if isinstance(c, (int, np.integer)) and int(c) in wanted:
            geno_cols.append(c)
        elif cs.isdigit() and int(cs) in wanted:
            geno_cols.append(c)
        elif re.fullmatch(r"(\d+)_Depth", cs, flags=re.IGNORECASE):
            sid = int(re.fullmatch(r"(\d+)_Depth", cs, flags=re.IGNORECASE).group(1))
            if sid in wanted:
                depth_cols.append(c)

    keep_cols = ["contig", "pos"]
    if "ref" in df.columns:
        keep_cols.append("ref")
    keep_cols += geno_cols
    keep_cols += depth_cols

    df = df[keep_cols].copy()

    df["contig_norm"] = df["contig"].apply(lambda x: normalize_chrom(x) or "NA")
    df = df[df["contig_norm"].isin(MAIN_CHROMS)].copy()

    df["pos"] = pd.to_numeric(df["pos"], errors="coerce").astype("Int64")
    df = df.dropna(subset=["pos"])
    df["pos"] = df["pos"].astype(int)

    return df


# -----------------------------
# Analyses
# -----------------------------

def compute_callability_qc(s4: pd.DataFrame, meta: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    total_loci = int(len(s4))

    for sid in meta["sample_id"].unique().tolist():
        geno_col = sid if sid in s4.columns else str(sid) if str(sid) in s4.columns else None
        if geno_col is None:
            continue
        geno = s4[geno_col]
        callable_mask = geno.apply(is_callable_genotype).to_numpy()
        callable_loci = int(callable_mask.sum())
        missing_rate = float(1.0 - (callable_loci / total_loci))

        mrow = meta.loc[meta["sample_id"] == sid].iloc[0]
        rows.append({
            "sample_id": sid,
            "rad_type": str(mrow["rad_type"]),
            "dose_Gy": float(mrow["dose_Gy"]) if pd.notna(mrow["dose_Gy"]) else np.nan,
            "origin": str(mrow["origin"]) if "origin" in meta.columns else "",
            "parent_id": int(mrow["parent_id"]) if pd.notna(mrow["parent_id"]) else np.nan,
            "callable_loci": callable_loci,
            "panel_loci_chr1_10": total_loci,
            "missing_rate_panel": missing_rate,
        })

    qc = pd.DataFrame(rows)

    qc_group = (qc.groupby("rad_type")
                .agg(
                    n_samples=("sample_id", "nunique"),
                    callable_mean=("callable_loci", "mean"),
                    callable_median=("callable_loci", "median"),
                    missing_mean=("missing_rate_panel", "mean"),
                    missing_median=("missing_rate_panel", "median"),
                ).reset_index())

    return qc, qc_group


def annotate_panel_loci(s4: pd.DataFrame, gene_gff: Path, repeat_gff: Path, promoter_upstream_bp: int) -> pd.DataFrame:
    cds, gene, promoter = parse_gene_gff(gene_gff, promoter_upstream_bp)
    te = parse_te_gff(repeat_gff)

    loci = s4[["contig_norm", "pos"]].copy()
    loci["contig_norm"] = loci["contig_norm"].astype(str)
    loci["pos"] = loci["pos"].astype(int)

    loci["is_CDS"] = False
    loci["is_Promoter"] = False
    loci["is_GeneBody"] = False
    loci["is_TE"] = False

    for chrom in MAIN_CHROMS:
        idx = loci.index[loci["contig_norm"] == chrom].to_numpy()
        if idx.size == 0:
            continue
        pos = loci.loc[idx, "pos"].to_numpy(dtype=np.int64)

        cds_s, cds_e = build_arrays(cds.get(chrom, []))
        gen_s, gen_e = build_arrays(gene.get(chrom, []))
        pro_s, pro_e = build_arrays(promoter.get(chrom, []))
        te_s, te_e = build_arrays(te.get(chrom, []))

        loci.loc[idx, "is_CDS"] = membership_mask(pos, cds_s, cds_e)
        loci.loc[idx, "is_Promoter"] = membership_mask(pos, pro_s, pro_e)
        loci.loc[idx, "is_GeneBody"] = membership_mask(pos, gen_s, gen_e)
        loci.loc[idx, "is_TE"] = membership_mask(pos, te_s, te_e)

    primary = np.full(len(loci), "Intergenic", dtype=object)
    primary[loci["is_GeneBody"].to_numpy()] = "GeneBody_nonCDS"
    primary[loci["is_Promoter"].to_numpy()] = "Promoter"
    primary[loci["is_CDS"].to_numpy()] = "CDS"
    loci["primary_feature"] = primary

    return loci


def per_sample_callable_feature_fractions(s4: pd.DataFrame, loci_annot: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    feat = loci_annot["primary_feature"].to_numpy()
    te_flag = loci_annot["is_TE"].to_numpy(dtype=bool)

    out = []
    for sid in meta["sample_id"].unique().tolist():
        geno_col = sid if sid in s4.columns else str(sid) if str(sid) in s4.columns else None
        if geno_col is None:
            continue
        geno = s4[geno_col]
        callable_mask = geno.apply(is_callable_genotype).to_numpy()
        total_callable = int(callable_mask.sum())
        if total_callable == 0:
            continue

        row = meta.loc[meta["sample_id"] == sid].iloc[0]
        rec = {
            "sample_id": sid,
            "rad_type": str(row["rad_type"]),
            "dose_Gy": float(row["dose_Gy"]) if pd.notna(row["dose_Gy"]) else np.nan,
            "origin": str(row["origin"]) if "origin" in meta.columns else "",
            "parent_id": int(row["parent_id"]) if pd.notna(row["parent_id"]) else np.nan,
            "callable_total": total_callable,
            "callable_TE": int((callable_mask & te_flag).sum()),
            "callable_TE_fraction": float((callable_mask & te_flag).sum() / total_callable),
        }

        for f in ["CDS", "Promoter", "GeneBody_nonCDS", "Intergenic"]:
            rec[f"callable_{f}"] = int((callable_mask & (feat == f)).sum())
            rec[f"callable_{f}_fraction"] = float(rec[f"callable_{f}"] / total_callable)

        out.append(rec)

    return pd.DataFrame(out)


def compute_induced_feature_enrichment(events: pd.DataFrame, loci_annot: pd.DataFrame, callable_feats: pd.DataFrame) -> pd.DataFrame:
    key_loci = loci_annot.copy()
    key_loci["key"] = key_loci["contig_norm"].astype(str) + ":" + key_loci["pos"].astype(str)

    ev = events.copy()
    ev["key"] = ev["contig_norm"].astype(str) + ":" + ev["pos"].astype(int).astype(str)

    ev = ev.merge(key_loci[["key", "primary_feature", "is_TE"]], on="key", how="left")
    ev["primary_feature"] = ev["primary_feature"].fillna("Intergenic")
    ev["is_TE"] = ev["is_TE"].fillna(False)

    out_rows = []
    for rad, sub in ev.groupby("rad_type"):
        N = int(len(sub))
        if N == 0:
            continue
        cf = callable_feats[callable_feats["rad_type"] == rad]
        if cf.empty:
            continue

        exp_frac = {
            "CDS": float(cf["callable_CDS_fraction"].mean()),
            "Promoter": float(cf["callable_Promoter_fraction"].mean()),
            "GeneBody_nonCDS": float(cf["callable_GeneBody_nonCDS_fraction"].mean()),
            "Intergenic": float(cf["callable_Intergenic_fraction"].mean()),
        }
        exp_te = float(cf["callable_TE_fraction"].mean())

        obs_counts = sub["primary_feature"].value_counts().to_dict()
        obs_te = int(sub["is_TE"].sum())

        for feat, p0 in exp_frac.items():
            obs = int(obs_counts.get(feat, 0))
            exp = N * p0
            enr = (obs / N) / p0 if p0 > 0 else np.nan
            pval = stats.binomtest(obs, N, p0, alternative="two-sided").pvalue if p0 > 0 else np.nan
            out_rows.append({
                "rad_type": rad,
                "feature": feat,
                "N_events": N,
                "observed": obs,
                "expected": exp,
                "enrichment_obs_over_exp": enr,
                "p0_callable_fraction": p0,
                "pvalue_binom": pval,
            })

        obs_te_frac = obs_te / N if N > 0 else np.nan
        enr_te = obs_te_frac / exp_te if exp_te > 0 else np.nan
        pval_te = stats.binomtest(obs_te, N, exp_te, alternative="two-sided").pvalue if exp_te > 0 else np.nan
        out_rows.append({
            "rad_type": rad,
            "feature": "TE_overlap",
            "N_events": N,
            "observed": obs_te,
            "expected": N * exp_te,
            "enrichment_obs_over_exp": enr_te,
            "p0_callable_fraction": exp_te,
            "pvalue_binom": pval_te,
        })

    return pd.DataFrame(out_rows)


def compute_cluster_fraction(events: pd.DataFrame, window_bp: int) -> pd.DataFrame:
    ev = events.sort_values(["sample_id", "contig_norm", "pos"]).copy()
    ev["prev_pos"] = ev.groupby(["sample_id", "contig_norm"])["pos"].shift(1)
    ev["next_pos"] = ev.groupby(["sample_id", "contig_norm"])["pos"].shift(-1)
    ev["d_prev"] = ev["pos"] - ev["prev_pos"]
    ev["d_next"] = ev["next_pos"] - ev["pos"]
    ev["clustered"] = (ev["d_prev"] <= window_bp) | (ev["d_next"] <= window_bp)

    g = ev.groupby(["sample_id", "rad_type"]).agg(
        total_events=("pos", "count"),
        clustered_events=("clustered", "sum"),
    ).reset_index()
    g["cluster_fraction"] = g["clustered_events"] / g["total_events"]
    return g


def confounding_crosstab(meta: pd.DataFrame) -> pd.DataFrame:
    if "origin" not in meta.columns:
        return pd.DataFrame()
    ct = pd.crosstab(meta["origin"], meta["rad_type"]).reset_index().rename(columns={"origin": "origin"})
    return ct


# -----------------------------
# Plot helpers (no numbers in figures)
# -----------------------------

def save_fig(fig: plt.Figure, out_base: Path, dpi: int = 300):
    fig.savefig(str(out_base) + ".pdf", format="pdf")
    fig.savefig(str(out_base) + ".png", format="png", dpi=dpi)
    plt.close(fig)


def fig_boxstrip(df: pd.DataFrame, x: str, y: str, out_base: Path, title: str, ylabel: str):
    fig, ax = plt.subplots(figsize=(5.0, 4.2))
    sns.boxplot(data=df, x=x, y=y, ax=ax)
    sns.stripplot(data=df, x=x, y=y, ax=ax, color="black", size=3, alpha=0.5)
    ax.set_title(title)
    ax.set_xlabel(x)
    ax.set_ylabel(ylabel)
    fig.tight_layout()
    save_fig(fig, out_base)


def fig_callable_feature_fractions(callable_feats: pd.DataFrame, out_base: Path):
    parts = []
    for f in ["CDS", "Promoter", "GeneBody_nonCDS", "Intergenic"]:
        col = f"callable_{f}_fraction"
        tmp = callable_feats[["rad_type", col]].copy()
        tmp = tmp.rename(columns={col: "fraction"})
        tmp["feature"] = f
        parts.append(tmp)
    tmp_te = callable_feats[["rad_type", "callable_TE_fraction"]].copy()
    tmp_te = tmp_te.rename(columns={"callable_TE_fraction": "fraction"})
    tmp_te["feature"] = "TE_overlap"
    parts.append(tmp_te)
    ff = pd.concat(parts, ignore_index=True)

    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    sns.boxplot(data=ff, x="feature", y="fraction", hue="rad_type", ax=ax)
    ax.set_title("Callable-space feature fractions (Table S4)")
    ax.set_xlabel("Feature")
    ax.set_ylabel("Callable fraction")
    fig.tight_layout()
    save_fig(fig, out_base)


def fig_enrichment_bar(enr: pd.DataFrame, out_base: Path):
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    plot_df = enr.copy()
    plot_df = plot_df[plot_df["feature"].isin(["CDS", "Promoter", "GeneBody_nonCDS", "Intergenic", "TE_overlap"])]
    sns.barplot(data=plot_df, x="feature", y="enrichment_obs_over_exp", hue="rad_type", ax=ax)
    ax.axhline(1.0, linestyle="--", linewidth=1)
    ax.set_title("Induced SNV enrichment by feature (callable-space corrected)")
    ax.set_xlabel("Feature")
    ax.set_ylabel("Observed / expected")
    fig.tight_layout()
    save_fig(fig, out_base)


# -----------------------------
# Main
# -----------------------------

def main():
    ap = argparse.ArgumentParser(description="Callable-space controls for induced SNV analyses (GBS panel).")
    ap.add_argument("--events", required=True, type=str, help="Path to events_long.csv(.gz)")
    ap.add_argument("--meta", required=True, type=str, help="Path to metadata CSV or Excel (S1_Sample_Metadata)")
    ap.add_argument("--s4", required=True, type=str, help="Path to Table S4.xlsx")
    ap.add_argument("--gene_gff", required=True, type=str, help="Path to gene.gff3(.gz)")
    ap.add_argument("--repeat_gff", required=True, type=str, help="Path to repeatmasked gff3(.gz)")
    ap.add_argument("--outdir", required=True, type=str, help="Output directory")
    ap.add_argument("--s4_header_row", default=2, type=int, help="0-based header row index (default 2).")
    ap.add_argument("--promoter_upstream_bp", default=2000, type=int, help="Promoter upstream length (bp).")
    ap.add_argument("--cluster_bp", default=10, type=int, help="Clustering distance for QC.")
    args = ap.parse_args()

    sns.set_style("whitegrid")

    outdir = Path(args.outdir)
    out_tables = outdir / "tables"
    out_figs = outdir / "figures"
    out_tables.mkdir(parents=True, exist_ok=True)
    out_figs.mkdir(parents=True, exist_ok=True)

    events = read_events(Path(args.events))
    meta = read_meta(Path(args.meta))

    merge_cols = ["sample_id", "dose_Gy", "parent_id"]
    if "origin" in meta.columns:
        merge_cols.append("origin")
    events = events.merge(meta[merge_cols], on="sample_id", how="left")

    sids = meta["sample_id"].astype(int).unique().tolist()
    parent_ids = meta["parent_id"].dropna().astype(int).unique().tolist()
    needed = sorted(set(sids + parent_ids))

    s4 = read_table_s4(Path(args.s4), needed, args.s4_header_row)

    qc, qc_group = compute_callability_qc(s4, meta)
    qc.to_csv(out_tables / "qc_callability_per_sample.csv", index=False)
    qc_group.to_csv(out_tables / "qc_callability_by_group.csv", index=False)

    fig_boxstrip(qc, "rad_type", "callable_loci", out_figs / "fig_callability_by_group",
                 "Callable loci per sample (Table S4)", "Callable loci")
    fig_boxstrip(qc, "rad_type", "missing_rate_panel", out_figs / "fig_missingness_by_group",
                 "Missingness per sample (Table S4)", "Missing rate")

    loci_annot = annotate_panel_loci(s4, Path(args.gene_gff), Path(args.repeat_gff), args.promoter_upstream_bp)
    callable_feats = per_sample_callable_feature_fractions(s4, loci_annot, meta)
    callable_feats.to_csv(out_tables / "callable_feature_denominators_per_sample.csv", index=False)
    fig_callable_feature_fractions(callable_feats, out_figs / "fig_callable_feature_fractions")

    enr = compute_induced_feature_enrichment(events, loci_annot, callable_feats)
    enr.to_csv(out_tables / "induced_feature_enrichment_callable_corrected.csv", index=False)
    fig_enrichment_bar(enr, out_figs / "fig_enrichment_by_feature")

    ct = confounding_crosstab(meta)
    if not ct.empty:
        ct.to_csv(out_tables / "confounding_origin_by_radtype.csv", index=False)

    cf = compute_cluster_fraction(events, window_bp=args.cluster_bp)
    cf.to_csv(out_tables / "cluster_fraction_per_sample.csv", index=False)
    fig_boxstrip(cf, "rad_type", "cluster_fraction", out_figs / "fig_cluster_fraction_by_group",
                 "Fraction of events in local clusters", "Cluster fraction")

    cfg = {
        "events": str(Path(args.events)),
        "meta": str(Path(args.meta)),
        "s4": str(Path(args.s4)),
        "gene_gff": str(Path(args.gene_gff)),
        "repeat_gff": str(Path(args.repeat_gff)),
        "outdir": str(outdir),
        "s4_header_row": args.s4_header_row,
        "promoter_upstream_bp": args.promoter_upstream_bp,
        "cluster_bp": args.cluster_bp,
        "callable_definition": "ACGT + IUPAC(RYSWKM) treated as callable; N/empty treated as missing",
        "n_events_chr1_10": int(len(events)),
        "n_meta_rows": int(len(meta)),
        "n_s4_panel_loci_chr1_10": int(len(s4)),
    }
    (out_tables / "analysis_config.json").write_text(json.dumps(cfg, indent=2), encoding="utf-8")

    print("Done.")
    print(str(outdir))


if __name__ == "__main__":
    main()

