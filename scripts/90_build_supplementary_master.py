# -*- coding: utf-8 -*-
"""
90_build_supplementary_master.py

Combine final CSV tables into a single Excel workbook suitable for submission.

Default output:
  outputs/Supplementary_Data_Master_REGENERATED.xlsx
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from pipeline_utils import resolve_project_paths


def read_csv_required(path: Path, **kwargs) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Required table missing: {path}")
    return pd.read_csv(path, **kwargs)


def read_csv_optional(path: Path, **kwargs) -> pd.DataFrame | None:
    if not path.exists():
        return None
    return pd.read_csv(path, **kwargs)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", default=None)
    ap.add_argument("--out", default=None, help="Output xlsx path (default: outputs/Supplementary_Data_Master_REGENERATED.xlsx).")
    args = ap.parse_args()

    paths = resolve_project_paths(args.project_root)

    # Locate tables
    rad_tab = paths.outputs_radclock / "tables"
    fun_tab = paths.outputs_functional / "tables"

    tables = {
        "S1_Sample_Metadata": read_csv_required(rad_tab / "per_mutant_radclock.csv"),
        "S2_Spatial_Stats": read_csv_required(rad_tab / "spatial_statistics.csv"),
        "S3_Track_Analysis": read_csv_required(rad_tab / "track_summary_per_sample.csv"),
        "S4_MBS_Statistics": read_csv_required(rad_tab / "mbs_statistics.csv"),
        "S5_Shielding_Summary": read_csv_required(fun_tab / "shielding_stats_corrected.csv"),
        "S7_Evolutionary_Stats": read_csv_required(fun_tab / "evolutionary_stats.csv"),
        "S8_TE_Mobilization": read_csv_required(fun_tab / "te_mobilization_stats.csv"),
        "S9_Mutation_Spectrum_Sub6": read_csv_required(rad_tab / "sub6_fractions.csv"),
        "S10_Signature_Matrix_96": read_csv_required(rad_tab / "sub96_matrix.csv"),
    }

    dose = read_csv_optional(fun_tab / "shielding_corrected_by_dose.csv")
    if dose is not None:
        tables["S6_Shielding_Dose"] = dose

    # Output path
    out_path = Path(args.out) if args.out else (paths.outputs / "Supplementary_Data_Master_REGENERATED.xlsx")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
        # Stable order (S1..S10)
        for sheet in sorted(tables.keys(), key=lambda s: (int(s.split("_")[0][1:]), s)):
            df = tables[sheet]
            df.to_excel(writer, sheet_name=sheet, index=False)

    print(f"[OK] Saved master workbook: {out_path}")
    print("[INFO] Sheets written:")
    for k in sorted(tables.keys(), key=lambda s: (int(s.split('_')[0][1:]), s)):
        print(f"  - {k} ({tables[k].shape[0]} rows, {tables[k].shape[1]} cols)")


if __name__ == "__main__":
    main()
