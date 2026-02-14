# Sorghum Mutant (Proton vs Gamma) – Reproducible Tables-Only Pipeline

This package regenerates the **tables / supplementary master workbook** for the Sorghum proton-beam vs gamma-ray mutagenesis manuscript.  
 
## Folder layout

```text
Sorghum_Proton_Mutagenesis/
  scripts/
    00_run_pipeline_tables_only.py
    10_radclock_pipeline_tables.py
    20_spatial_stats_1mb.py
    25_cluster_sensitivity_recurrent_filter.py
    30_track_calling_500kb.py
    40_mbs_statistics.py
    50_annotate_events_functional.py
    60_shielding_corrected_callable_space.py
    61_shielding_by_dose_optional.py
    65_callable_space_controls_defense.py
    70_evolutionary_indifference.py
    80_te_enrichment.py
    90_build_supplementary_master.py
  inputs/
    Table S2.xlsx
    Table S4.xlsx
    SbicolorRio_468_v2.1.gene.gff3.gz
    SbicolorRio_468_v2.1.repeatmasked_assembly_v2.0.gff3.gz
  outputs_radclock/                
  outputs_functional_shielding/    
  outputs/                         
Input filename tolerance
Some labs use underscores (Table_S4.xlsx) instead of spaces (Table S4.xlsx).
This pipeline auto-detects common variants.

Environment
Minimum Python: 3.10+ (tested with 3.11)

Install dependencies (conda):

conda create -n sorghum_mut -y python=3.11
conda activate sorghum_mut
pip install -r requirements.txt
Run everything (tables only)
From the project root:

python scripts/00_run_pipeline_tables_only.py
Outputs:

Per-step CSV tables in outputs_radclock/tables/ and outputs_functional_shielding/tables/

Final master workbook: outputs/Supplementary_Data_Master_REGENERATED.xlsx

Optional: also compute “by dose” shielding table (S6)
python scripts/00_run_pipeline_tables_only.py --include-shielding-by-dose
Troubleshooting
If you see “Required input file not found”, check that your files are in inputs/ and match one of the expected names.

The heaviest steps are functional annotation (GFF3 parsing) and callable-space correction (scanning Table S4).

Robustness & QC (optional)
These scripts provide additional quality-control and sensitivity analyses.
They are optional and do not change the main regenerated tables/workbook produced by the core pipeline.

Step 25 (25_cluster_sensitivity_recurrent_filter.py): tests whether local clustering persists after filtering recurrent loci.

Step 65 (65_callable_space_controls_defense.py): callability QC + callable denominators + callable-space-corrected feature enrichment.

Run after the main pipeline completes:

# Step 25: Cluster sensitivity check
python scripts/25_cluster_sensitivity_recurrent_filter.py \
  --events outputs_radclock/tables/events_long.csv.gz \
  --outdir outputs_robustness/cluster_sensitivity

# Step 65: Callable-space QC + denominators + corrected enrichment
python scripts/65_callable_space_controls_defense.py \
  --events outputs_radclock/tables/events_long.csv.gz \
  --meta outputs_radclock/tables/mutants_metadata.csv \
  --s4 "inputs/Table S4.xlsx" \
  --gene_gff inputs/SbicolorRio_468_v2.1.gene.gff3.gz \
  --repeat_gff inputs/SbicolorRio_468_v2.1.repeatmasked_assembly_v2.0.gff3.gz \
  --outdir outputs_robustness/callable_space_controls
