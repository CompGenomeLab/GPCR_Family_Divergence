# Benchmark Folder Organization

## Current Structure

```
benchmark/
├── Consolidated_GPCR_mapping.tsv    # Single unified GPCR mapping file
├── Fig2cde.html                     # Combined horizontal visualization
├── FigS1_DMS_benchmark.html         # Deep Mutational Scanning benchmarks
├── FigS2_correlation_plots.html     # Correlation analysis plots
├── FigS3_evolutiıon_benchmark.html  # Evolution-based benchmarks
├── README.md                        # Main documentation
├── MAPPING_CONSOLIDATION_SUMMARY.md # Mapping consolidation details
├── FOLDER_ORGANIZATION.md           # This file
└── data/                            # All data files
    ├── ADRB2_BetaArrestin_EC100.tsv
    ├── ADRB2_Emax_Gprotein.tsv
    ├── ADRB2_potency_efficacy.tsv
    ├── ADRB2_tolerance_data_JonesEtAl.xls
    ├── ClinVar_aggregated_pathogenicity.txt
    ├── Combined_mutational_intolerance.tsv
    ├── EC100_alanine_scanning.txt
    ├── GPR68_normalized_activation_scores.tsv
    ├── GPR68_sensitivity_delta.tsv
    ├── HRH2_entropy_conservation.txt
    ├── MC4R_sensitivity.tsv
    ├── Sanders_ligand_binding_residues.txt
    ├── V2R_expression.tsv
    ├── variant_summary_classA_GPCR*.txt (3 files)
    ├── uniprot_hgcn_mapping_classA.txt
    ├── calculate_combined_significance.py
    ├── clinvar_data_manipulation.py
    └── MISSING_FILES.md
```

## File Organization

### Main Folder (benchmark/)
- **HTML Visualizations**: All interactive plots (4 files)
- **Mapping File**: `Consolidated_GPCR_mapping.tsv` - unified mapping for all receptors
- **Documentation**: README and summary files

### Data Subfolder (data/)
- **DMS Data**: Deep Mutational Scanning datasets (ADRB2, MC4R, GPR68, V2R)
- **Clinical Data**: ClinVar pathogenicity data
- **Evolution Data**: Sanders ligand binding, HRH2 conservation/entropy
- **Analysis Scripts**: Python scripts for data processing
- **Variant Data**: ClinVar variant summaries

## Path Updates

All HTML files have been updated to use `data/` prefix for data file references:
- `HRH2_entropy_conservation.txt` → `data/HRH2_entropy_conservation.txt`
- `MC4R_sensitivity.tsv` → `data/MC4R_sensitivity.tsv`
- And all other data files...

The mapping file (`Consolidated_GPCR_mapping.tsv`) remains in the main folder for easy access.

## Python Script Updates

The `calculate_combined_significance.py` script has been updated to:
- Read mapping file from `../Consolidated_GPCR_mapping.tsv`
- Read all data files from the current directory (data/)

## Moving the Benchmark Folder

To move this folder outside of `scatter_plots_and_CR_SR_calculation/`:

```powershell
# From the Supplementary Information directory
Move-Item "code_and_data\scatter_plots_and_CR_SR_calculation\benchmark" "code_and_data\benchmark"
```

Note: If OneDrive sync is active, you may need to pause sync before moving the folder.


