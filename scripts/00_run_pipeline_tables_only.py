# -*- coding: utf-8 -*-
"""
00_run_pipeline_tables_only.py

End-to-end orchestration to regenerate all **tables** (no figures by default).

Usage (from project root):
  python scripts/00_run_pipeline_tables_only.py
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))

from pipeline_utils import resolve_project_paths


def run_step(script_path: Path, project_root: Path, extra_args: list[str] | None = None) -> None:
    cmd = [sys.executable, str(script_path), "--project-root", str(project_root)]
    if extra_args:
        cmd.extend(extra_args)
    print("\n" + "=" * 88)
    print("[STEP]", script_path.name)
    print(" ".join(cmd))
    print("=" * 88)
    subprocess.run(cmd, check=True, cwd=str(project_root))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", default=None, help="Project root (default: current working directory).")
    ap.add_argument("--max-rows", default=None, type=int, help="DEBUG: limit scanned loci in Table S4 (RadClock only).")
    ap.add_argument("--include-shielding-by-dose", action="store_true", help="Also compute shielding-by-dose table (S6).")
    args = ap.parse_args()

    paths = resolve_project_paths(args.project_root)
    project_root = paths.root

    scripts_dir = project_root / "scripts"

    # Steps (tables-only)
    steps = [
        ("10_radclock_pipeline_tables.py", ["--no-figures"] + (["--max-rows", str(args.max_rows)] if args.max_rows else [])),
        ("20_spatial_stats_1mb.py", []),
        ("30_track_calling_500kb.py", []),
        ("40_mbs_statistics.py", []),
        ("50_annotate_events_functional.py", []),
        ("60_shielding_corrected_callable_space.py", []),
        ("70_evolutionary_indifference.py", []),
        ("80_te_enrichment.py", []),
    ]

    for fname, extra in steps:
        run_step(scripts_dir / fname, project_root, extra)

    if args.include_shielding_by_dose:
        run_step(scripts_dir / "61_shielding_by_dose_optional.py", project_root, [])

    # Build master workbook at the end
    run_step(scripts_dir / "90_build_supplementary_master.py", project_root, [])

    print("\n✅ Pipeline finished.")
    print(f"   - outputs_radclock/           : {paths.outputs_radclock}")
    print(f"   - outputs_functional_shielding: {paths.outputs_functional}")
    print(f"   - outputs/                    : {paths.outputs}")


if __name__ == "__main__":
    main()
