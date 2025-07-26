# Decoding Functional Specialization in GPCRs through Evolution-Guided Residue Profiling - Supplementary Information

This directory contains supplementary information for the manuscript "Decoding Functional Specialization in GPCRs through Evolution-Guided Residue Profiling".

## Directory Structure

```text
Supplementary Information/
├── README.md                           # This file
├── Supplementary_Table.xlsx            # Comprehensive data table with all annotations
├── Figures/                            # Manuscript figures
│   ├── Figure1.png                    # Automated ortholog identification workflow
│   ├── Figure2.png                    # Class A, olfactory, and Class T analysis
│   ├── Figure3.png                    # Class B1 and B2 analysis
│   ├── Figure4.png                    # Class F analysis
│   ├── Figure5.png                    # Class C analysis
│   └── Figure6.png                    # Subtype-specific selective residues in ligand and transducer selectivity
└── code_and_data/                     # Computational analysis files
    ├── ortholog_pipeline/             # Ortholog identification pipeline
    ├── CRs_and_SRs_on_structure/     # Structural visualizations
    └── scatter_plots_and_CR_SR_calculation/  # Scatter plot generation
```

## Supplementary Table

**Supplementary_Table.xlsx** contains comprehensive data from the GPCR conservation analysis, including:

### Data Content
- **Scatter Plot Information**: Conservation percentages and entropy values for all residues
- **Initial Annotations**: Detailed functional annotations assigned during analysis
- **Simplified Annotations**: Streamlined functional categories for visualization
- **CR/SR Information**: Classification of residues as Common (CR) or Selective (SR)
- **Class-Specific Data**: Separate sheets for each GPCR class (A, B1, B2, C, F, T, Olfactory)

### Simplified Annotation Categories
- **Ligand**: Residues involved in ligand binding
- **Transducer**: Residues involved in G protein/transducer binding
- **Known Motifs**: Conserved motifs (CWxP, NPxxY, PIF, DRY, Na+ Pocket)
- **Disulfide Bridge**: Cysteine residues forming disulfide bonds
- **Cholesterol**: Cholesterol binding site residues
- **ICL2**: Intracellular loop 2 residues (Class T specific)
- **Tethered Agonist**: Residues involved in tethered agonist binding (Class B2)
- **VFT Ligand**: Venus Flytrap domain ligand binding residues (Class C)
- **Allosteric Modulator**: Residues involved in allosteric modulation
- **WNT Binding**: Residues involved in WNT protein binding (Class F)
- **Other**: Residues without specific functional annotation

### Data Organization
Each GPCR class has its own worksheet containing:
- Residue numbers (using class-specific reference proteins)
- GPCRdb numbers for the given reference receptor
- Family-wide conservation percentages across family members
- Family-wide entropy values for sequence variability
- Initial detailed annotations
- Simplified functional categories
- CR/SR classification with thresholds

## Figures

The **Figures/** directory contains all manuscript figures:

- **Figure 1**: Automated ortholog identification and conservation analysis workflow
- **Figure 2**: Family-wide conservation and structural mapping for Class A, olfactory, and Class T receptors
- **Figure 3**: Conservation and structural organization for Class B1 and adhesion (B2) GPCRs
- **Figure 4**: Evolutionary conservation patterns in Class F GPCRs across distinct domains
- **Figure 5**: Evolutionary conservation patterns in Class C GPCRs across different domains
- **Figure 6**: Subtype-specific selective residues in ligand and transducer selectivity

## Code and Data

The **code_and_data/** directory contains all computational analysis files:

### ortholog_pipeline/
- Complete pipeline for identifying GPCR orthologs
- Python scripts for sequence processing and phylogenetic analysis
- Human protein sequences organized by GPCR class
- SLURM scripts for cluster execution

### CRs_and_SRs_on_structure/
- Structural visualizations of CR/SR mapping
- PyMOL session files for interactive 3D analysis
- AlphaFold2 models and experimental structures
- PNG images used as manuscript figures

### scatter_plots_and_CR_SR_calculation/
- Scripts for generating scatter plots
- CR/SR calculation algorithms
- Label files with simplified annotations
- Special analysis for Class F GPCRs
- **Class alignments**: Multiple sequence alignments used for receptor mapping and analysis

## Usage

### Accessing Data
1. Open **Supplementary_Table.xlsx** for comprehensive analysis data
2. Navigate to specific class worksheets for detailed information
3. Use conservation and entropy values for further analysis

### Viewing Figures
1. Open PNG files in any image viewer
2. Figures correspond to manuscript figures with same numbering

### Running Analysis
1. Navigate to **code_and_data/** for computational scripts
2. Follow README files in each subdirectory for specific instructions
3. Use provided SLURM scripts for cluster execution

## Alignment Strategy

The analysis uses representative sequences ("reps") and multiple sequence alignments to determine receptor sets and mapping:

### Representative Sequence Selection
- CD-HIT clustering identifies representative sequences from ortholog sets
- Top 5 largest clusters are selected as representatives
- Representatives provide evolutionary diversity while reducing computational complexity

### Alignment Process
1. **Initial Alignment**: Representative sequences + human sequences aligned using MAFFT
2. **Human Extraction**: Human sequences extracted from the alignment
3. **Final Analysis**: Family-wide analysis performed on human-only alignments

### Alignment Types
- **Full alignments**: Complete receptor sequences with representatives
- **Human-only alignments**: Extracted human sequences for family-wide analysis
- **Domain-specific alignments**: Specific regions (VFT domain, CRD domain, etc.)

### Impact of Alignment Choice
Different alignments produce different outcomes:
- **Receptor Set**: Different representative sequences change receptor composition
- **Mapping**: Alignment positions determine residue numbering and conservation calculations
- **Analysis Scope**: Domain-specific alignments focus on particular functional regions

### Alignment Files
The `class_alignments/` directory contains:
- Standard alignments with representatives
- Human-only extracted alignments
- Domain-specific alignments for particular regions
- Scripts for extracting human sequences from representative alignments

## Data Availability

All data and code are available at: https://github.com/CompGenomeLab/GPCR_Family_Divergence

## Contact

For questions or issues with the data or code, please contact:
- **Ogün Adebali**: oadebali@sabanciuniv.edu
