# Sorghum Proton Mutagenesis: Tables-Only Reproducible Pipeline

This repository contains the Python scripts used to regenerate the tabular analyses for the sorghum proton-beam versus gamma-ray mutagenesis study. The workflow is designed for reduced-representation GBS data and uses explicit GBS-callable denominators for rate, enrichment, and spatial validation analyses.

The repository does not include the large input genotype or annotation files. Place the required input files in `inputs/` before running the pipeline.

## Required input files

Create an `inputs/` folder in the project root and place the following files there:

```text
inputs/
  Table S2.xlsx
  Table S4.xlsx
  SbicolorRio_468_v2.1.gene.gff3.gz
  SbicolorRio_468_v2.1.repeatmasked_assembly_v2.0.gff3.gz
```

Common filename variants with spaces or underscores are detected by the pipeline where applicable.

## Environment

Python 3.10 or later is recommended.

```bash
pip install -r requirements.txt
```

The table-only workflow requires `pandas`, `numpy`, `scipy`, `openpyxl`, and `tqdm`.

## Run the core tables-only pipeline

From the project root:

```bash
python scripts/00_run_pipeline_tables_only.py --project-root .
```

This regenerates the main event calls, spatial summaries, functional annotations, TE overlap summaries, and the supplementary master workbook.

Main outputs:

```text
outputs_radclock/tables/
outputs_functional_shielding/tables/
outputs/Supplementary_Data_Master_REGENERATED.xlsx
```

## Run additional spatial validation analyses

After the core pipeline completes, run:

```bash
python scripts/95_spatial_validation_tables.py --project-root .
```

This writes a single validation workbook:

```text
outputs_validation/Spatial_Validation_Tables.xlsx
```

The workbook includes gamma high-load comparisons, callable-density-normalized spatial summaries, callable-density-weighted spatial null analyses, track-definition sensitivity analyses, gamma dose-response regression parameters, and coding-region fraction models with radiation type terms.

## Optional robustness tables

The following scripts generate additional table-only robustness outputs and can be run after the core pipeline.

```bash
python scripts/25_cluster_sensitivity_recurrent_filter.py \
  --events outputs_radclock/tables/events_long.csv.gz \
  --outdir outputs_robustness/cluster_sensitivity
```

```bash
python scripts/65_callable_space_controls.py \
  --events outputs_radclock/tables/events_long.csv.gz \
  --meta outputs_radclock/tables/mutants_metadata.csv \
  --s4 "inputs/Table S4.xlsx" \
  --gene_gff inputs/SbicolorRio_468_v2.1.gene.gff3.gz \
  --repeat_gff inputs/SbicolorRio_468_v2.1.repeatmasked_assembly_v2.0.gff3.gz \
  --outdir outputs_robustness/callable_space_controls
```

On Windows Anaconda Prompt, use the same commands on one line or replace the line-continuation backslashes with `^`.

## Script overview

```text
scripts/00_run_pipeline_tables_only.py        Run the core table-generation workflow
scripts/10_radclock_pipeline_tables.py       Build induced event calls and mutation spectra tables
scripts/20_spatial_stats_1mb.py              Compute 1-Mb spatial statistics
scripts/25_cluster_sensitivity_recurrent_filter.py
                                             Recurrent-locus filtering sensitivity tables
scripts/30_track_calling_500kb.py            Call high-density mutation windows
scripts/40_mbs_statistics.py                 Multi-base substitution summary tables
scripts/50_annotate_events_functional.py     Annotate induced events with genomic features
scripts/60_shielding_corrected_callable_space.py
                                             Callable-space-corrected feature enrichment
scripts/61_shielding_by_dose_optional.py     Optional dose-stratified feature enrichment
scripts/65_callable_space_controls.py        Callable-space QC and denominator tables
scripts/70_cds_load_relationship.py          CDS-fraction versus mutation-load summary tables
scripts/80_te_enrichment.py                  TE-overlap enrichment summary tables
scripts/90_build_supplementary_master.py     Build regenerated supplementary master workbook
scripts/95_spatial_validation_tables.py      Additional spatial validation workbook
```

## Output policy

Generated outputs are intentionally not tracked in the repository. Re-run the commands above to regenerate outputs locally.
