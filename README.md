# Sorghum Mutant (Proton vs Gamma) – Reproducible Tables-Only Pipeline

This package regenerates the **tables / supplementary master workbook** for the Sorghum proton-beam vs gamma-ray mutagenesis manuscript.
**Figures are optional and disabled by default** (per your request).

## Folder layout (expected)

```
Sorghum_Mutant_Publication_2025_2/
  scripts/
    00_run_pipeline_tables_only.py
    10_radclock_pipeline_tables.py
    20_spatial_stats_1mb.py
    30_track_calling_500kb.py
    40_mbs_statistics.py
    50_annotate_events_functional.py
    60_shielding_corrected_callable_space.py
    61_shielding_by_dose_optional.py
    70_evolutionary_indifference.py
    80_te_enrichment.py
    90_build_supplementary_master.py
  inputs/
    Table S2.xlsx
    Table S4.xlsx
    SbicolorRio_468_v2.1.gene.gff3.gz
    SbicolorRio_468_v2.1.repeatmasked_assembly_v2.0.gff3.gz
  outputs_radclock/                (auto-created)
  outputs_functional_shielding/    (auto-created)
  outputs/                         (auto-created)
```

### Input filename tolerance

Some labs use underscores (`Table_S4.xlsx`) instead of spaces (`Table S4.xlsx`).
This pipeline **auto-detects** common variants.

## Environment

Minimum Python: **3.10+** (tested with 3.11)

Install dependencies (conda):

```bash
conda create -n sorghum_mut -y python=3.11
conda activate sorghum_mut
pip install -r requirements.txt
```

## Run everything (tables only)

From the project root:

```bash
python scripts/00_run_pipeline_tables_only.py
```

Outputs:
- Per-step CSV tables in `outputs_radclock/tables/` and `outputs_functional_shielding/tables/`
- Final master workbook:
  - `outputs/Supplementary_Data_Master_REGENERATED.xlsx`

## Optional: also compute “by dose” shielding table (S6)

```bash
python scripts/00_run_pipeline_tables_only.py --include-shielding-by-dose
```

## Troubleshooting

- If you see “Required input file not found”, check that your files are in `inputs/` and match one of the expected names.
- The heaviest steps are functional annotation (GFF3 parsing) and callable-space correction (scanning Table S4).

