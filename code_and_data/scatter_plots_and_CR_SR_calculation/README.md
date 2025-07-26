# Scatter Plots and CR/SR Calculation - Decoding Functional Specialization in GPCRs through Evolution-Guided Residue Profiling

This directory contains scripts and data files for generating scatter plots and calculating Common Residues (CRs) and Selective Residues (SRs) across GPCR classes, supporting the manuscript "Decoding Functional Specialization in GPCRs through Evolution-Guided Residue Profiling".

## File Organization

### Scatter Plot Files
- `classA_scatter.svg` - Scatter plot for Class A GPCRs
- `classB1_scatter.svg` - Scatter plot for Class B1 GPCRs  
- `classB2_scatter.svg` - Scatter plot for Class B2 (adhesion) GPCRs
- `classC_scatter.svg` - Scatter plot for Class C GPCRs
- `classF_scatter.svg` - Scatter plot for Class F GPCRs
- `classT_scatter.svg` - Scatter plot for Class T GPCRs
- `classOlf_scatter.svg` - Scatter plot for Olfactory GPCRs

### Label Files (Simplified Annotations)
- `classA_labels.txt` - CR/SR annotations for Class A (361 residues)
- `classB1_labels.txt` - CR/SR annotations for Class B1 (595 residues)
- `classB2_labels.txt` - CR/SR annotations for Class B2 (1449 residues)
- `classC_labels.txt` - CR/SR annotations for Class C (1080 residues)
- `classF_labels.txt` - CR/SR annotations for Class F (576 residues)
- `classT_labels.txt` - CR/SR annotations for Class T (340 residues)
- `classOlf_labels.txt` - CR/SR annotations for Olfactory (352 residues)

### Analysis Scripts
- `family_wide_conservation_entropy_calculate.py` - Calculates family-wide conservation and entropy metrics
- `scatter_plot_generate.py` - Generates scatter plots with CR/SR thresholds
- `calculate_aa_counts_and_percentage.py` - Calculates amino acid conservation percentages

### Class Alignments Directory
- `class_alignments/` - Contains multiple sequence alignments for each GPCR class

### Class F Special Analysis
- `CR_detection_classF/` - Special analysis folder for Class F GPCRs
  - `FZD_CR_detection_aligortihm.py` - Algorithm for detecting CRs in FZD receptors
  - `classF_conservation_FZD7.tsv` - Conservation data for class F
  - `classF_conservation_FZD7_CRs.tsv` - Identified CRs for class F

## Data Format

### Label Files Structure
Each label file contains tab-separated data with columns:
- `HRH2_res_no` - Residue number (using HRH2 as reference)
- `conservation` - Family-wide conservation percentage
- `entropy` - Shannon entropy value
- `annotation` - Functional annotation (Ligand, Transducer, Known Motifs, etc.)

### Functional Annotations
- **Ligand** - Residues involved in ligand binding
- **Transducer** - Residues involved in G protein/transducer binding
- **Known Motifs** - Conserved motifs (CWxP, NPxxY, PIF, DRY, Na+ Pocket)
- **Disulfide Bridge** - Cysteine residues forming disulfide bonds
- **Cholesterol** - Cholesterol binding site residues
- **ICL2** - Intracellular loop 2 residues (Class T specific)
- **Tethered Agonist** - Residues involved in tethered agonist binding (Class B2)
- **VFT Ligand** - Venus Flytrap domain ligand binding residues (Class C)
- **Allosteric Modulator** - Residues involved in allosteric modulation
- **WNT Binding** - Residues involved in WNT protein binding (Class F)
- **Other** - Residues without specific functional annotation

## Class F Special Analysis

The `CR_detection_classF/` folder contains specialized analysis for Class F GPCRs (Frizzled receptors), which required specific criteria due to their unique evolutionary patterns and high sequence similarity within the family.

## Usage

### Generating Scatter Plots
1. Run `scatter_plot_generate.py` with appropriate class-specific thresholds
2. Adjust conservation and entropy thresholds for each GPCR class
3. Output SVG files show CR/SR distribution with color-coded annotations

### Analyzing Conservation Data
1. Use `family_wide_conservation_entropy_calculate.py` to process alignment data
2. Apply class-specific thresholds to identify CRs and SRs
3. Generate label files with functional annotations

### Class F Analysis
1. Use `FZD_CR_detection_aligortihm.py` for specialized Class F analysis
2. Apply stricter criteria to identify most critical CRs
3. Consider both identical and similar amino acid conservation 