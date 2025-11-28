# GPCR Benchmark Visualizations

This folder contains interactive HTML visualizations for benchmarking GPCR functional data against multiple experimental datasets, including deep mutational scanning (DMS) studies, clinical variants, and evolutionary analyses.

## 📁 Folder Structure

```
benchmark/
├── *.html                           # 4 Interactive visualization files
├── Consolidated_GPCR_mapping.tsv    # Unified GPCR residue mapping
├── README.md                        # This file
├── MAPPING_CONSOLIDATION_SUMMARY.md # Technical details on mapping consolidation
├── FOLDER_ORGANIZATION.md           # Folder organization details
└── data/                            # All data files (19 files + 2 scripts)
```

## 📊 Visualization Files

### 1. **FigS1_DMS_benchmark.html**
Deep Mutational Scanning benchmark plots showing:
- ADRB2 tolerance data
- MC4R sensitivity scores  
- GPR68 pH-dependent activation
- V2R expression levels
- ADRB2 Emax (G-protein efficacy)
- ADRB2 β-arrestin recruitment
- EC100 alanine scanning data

**10 panels total** - scatter plots with entropy/conservation overlays

### 2. **FigS2_correlation_plots.html**
Correlation analysis between different experimental metrics:
- Tolerance vs. Sensitivity
- Delta vs. pH responses
- Expression correlations
- Emax vs. β-arrestin
- ClinVar pathogenicity
- EC100 potency/efficacy correlations

**20 panels total** (10 datasets × 2 correlation axes each)

### 3. **FigS3_evolutiıon_benchmark.html**
Evolution-based functional residue predictions:
- OPSD global/local toggle switches
- Sanders ligand binding residues

**2 panels** - comparing evolutionary predictions with experimental data

### 4. **Fig2cde.html**
Combined horizontal view showing:
- ClinVar pathogenicity scores
- Combined mutational intolerance
- Common activation network residues

**3 panels** - integrated view for manuscript figure

## 🗺️ Mapping File

**`Consolidated_GPCR_mapping.tsv`** - Single unified mapping file containing:
- HRH2 reference positions (residue numbers, GPCRdb IDs, regions)
- Mappings to 5 receptors: ADRB2, MC4R, GPR68, V2R, OPSD
- 361 rows covering all relevant positions

All HTML files use this single mapping file for consistency.

## 📦 Data Folder Contents

### DMS Datasets
- `ADRB2_tolerance_data_JonesEtAl.xls` - β2-adrenergic receptor tolerance
- `ADRB2_potency_efficacy.tsv` - Potency and efficacy measurements
- `ADRB2_Emax_Gprotein.tsv` - G-protein efficacy
- `ADRB2_BetaArrestin_EC100.tsv` - β-arrestin recruitment
- `MC4R_sensitivity.tsv` - Melanocortin-4 receptor sensitivity
- `GPR68_normalized_activation_scores.tsv` - GPR68 activation scores
- `V2R_expression.tsv` - Vasopressin receptor expression
- `EC100_alanine_scanning.txt` - Alanine scanning data

### Reference Data
- `HRH2_entropy_conservation.txt` - Entropy and conservation scores for HRH2 positions
- `Sanders_ligand_binding_residues.txt` - Evolutionarily predicted ligand binding residues
- `Combined_mutational_intolerance.tsv` - Integrated mutational intolerance scores

### Clinical Data
- `ClinVar_aggregated_pathogenicity.txt` - ClinVar pathogenicity aggregated by position
- `variant_summary_classA_GPCR.txt` - Class A GPCR variant summaries
- `variant_summary_classA_GPCR_missense_simplified.txt` - Simplified missense variants
- `variant_summary_classA_GPCR_missense_simplified_with_uniprot.txt` - With UniProt IDs
- `uniprot_hgcn_mapping_classA.txt` - UniProt to HGNC mapping

### Analysis Scripts
- `calculate_combined_significance.py` - Calculates combined mutational intolerance scores
- `clinvar_data_manipulation.py` - Processes ClinVar data

### Documentation
- `MISSING_FILES.md` - Notes on missing data files (OPSD_surface_expression.txt)

## 🚀 Usage

### Viewing Visualizations
1. Open any HTML file in a modern web browser (Chrome, Firefox, Edge)
2. All data is loaded locally - no internet connection required
3. Interactive features:
   - Hover over points for details
   - Click legend items to toggle datasets
   - Pan and zoom on plots

### Running Analysis Scripts

**Calculate Combined Significance:**
```bash
cd data/
python calculate_combined_significance.py
```

This will:
- Read all DMS datasets
- Map to HRH2 positions using the consolidated mapping
- Calculate combined mutational intolerance scores
- Output: `Combined_mutational_intolerance.tsv`

## 📋 Requirements

### For HTML Visualizations
- Modern web browser
- No additional dependencies (uses D3.js CDN)

### For Python Scripts
- Python 3.6+
- Required packages:
  ```
  xlrd  # For reading .xls files
  ```

## 🔧 Technical Details

### Path Configuration
All HTML files use relative paths with `data/` prefix:
```javascript
const CLASSA_LABELS = "data/HRH2_entropy_conservation.txt";
const MC4R_sensitivity_tsv_path = "data/MC4R_sensitivity.tsv";
```

The mapping file is in the parent directory:
```javascript
const CONSOLIDATED_mapping_tsv_path = "Consolidated_GPCR_mapping.tsv";
```

### Mapping File Structure
```
HRH2_region  HRH2_gpcrdb  HRH2_resNum  HRH2_AA  ADRB2_resNum  MC4R_resNum  GPR68_resNum  V2R_resNum  OPSD_resNum
N-term       2            2            A        21            23           7             
TM1          1x50         35           G        50            51           28            44           44
...
```

## 📖 References




