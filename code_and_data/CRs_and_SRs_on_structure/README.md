# CRs and SRs on Structure - GPCR Conservation Analysis - Decoding Functional Specialization in GPCRs through Evolution-Guided Residue Profiling

This directory contains structural visualizations and analysis files for Common Residues (CRs) and Selective Residues (SRs) across different GPCR classes, supporting the manuscript "Decoding Functional Specialization in GPCRs through Evolution-Guided Residue Profiling".

## Overview

This analysis distinguishes between two types of evolutionarily conserved residues in GPCRs:

- **Common Residues (CRs)**: Low entropy residues that support family-wide structural and functional foundations
- **Selective Residues (SRs)**: High entropy residues that confer specific functional characteristics within families

## File Organization

### Structural Models and Visualizations

#### Class A GPCRs
- `HRH2_classA_alphafold.pdb` - AlphaFold2 model of human HRH2 (histamine receptor)
- `HRH2_classA_alphafold.pse` - PyMOL session file with CR/SR mapping on HRH2
- `classA_CRs_and_SRs.png` - Visual representation of CR/SR mapping on HRH2 structure

#### Class B GPCRs
- `PTH1R_classB1_alphafold.pse` - PyMOL session for Class B1 (PTH1R model)
- `ClassB1_CRs_and_SRs.png` - Visual representation of CR/SR mapping on PTH1R structure
- `ADGRL3_classB2_alphafold.pse` - PyMOL session for Class B2 (ADGRL3 model)
- `classB2_CRs_and_SRs.png` - Visual representation of CR/SR mapping on ADGRL3 structure

#### Class C GPCRs
- `CASR_alphafold.pdb` - AlphaFold2 model of human CASR (calcium-sensing receptor)
- `CASR_ClassC_alphafold.pse` - PyMOL session with CR/SR mapping
- `ClassC_TMD_CR_SR.png` - Visual representation of transmembrane domain CR/SR mapping
- `ClassC_TMD_cholesterol.png` - Visual representation of cholesterol binding site
- `ClassC_TMD_motif.png` - Visual representation of conserved motifs
- `ClassC_VFT_CRs_SRs.png` - Visual representation of Venus Flytrap domain CR/SR mapping

#### Class F GPCRs
- `FZD7_alphafold.pdb` - AlphaFold2 model of human FZD7 (Frizzled-7)
- `FZD7_classF_alphafold.pse` - PyMOL session with CR/SR mapping
- `ClassF_cholesterol.png` - Visual representation of cholesterol binding site
- `ClassF_TMD_CR_SR.png` - Visual representation of transmembrane domain CR/SR mapping
- `ClassF_TMD_CR_SR_intracellular.png` - Visual representation of intracellular region CR/SR mapping
- `ClassF_WNT_CR_SR.png` - Visual representation of WNT binding interface CR/SR mapping

#### Class T (Taste) GPCRs
- `T2R39_alphafold.pdb` - AlphaFold2 model of human T2R39 (bitter taste receptor)
- `T2R39_classT_alphafold.pse` - PyMOL session with CR/SR mapping
- `ClassT_CRs_and_SRs.png` - Visual representation of CR/SR mapping on T2R39 structure

#### Olfactory Receptors
- `O52I2_alphafold.pdb` - AlphaFold2 model of human O52I2 (olfactory receptor)
- `O52I2_classOlf_alphafold.pse` - PyMOL session with CR/SR mapping
- `classOlf_CRs_and_SRs.png` - Visual representation of CR/SR mapping on O52I2 structure

### Experimental Structures
- `4f0a.cif` - WNT-FZD complex structure (PDB: 4F0A)
- `9epo.cif` - FZD7 structure with cholesterol (PDB: 9EPO)

## File Types

### PyMOL Session Files (.pse)
- Interactive 3D structural visualizations
- CRs represented as blue spheres
- SRs represented as red spheres
- Can be opened in PyMOL for detailed examination

### PNG Image Files (.png)
- Static visual representations of the PyMOL sessions
- Used as figures in the manuscript
- Show the same CR/SR mapping as the corresponding .pse files

### PDB Structure Files (.pdb)
- AlphaFold2 predicted structures
- Used as templates for CR/SR mapping
- Can be opened in any molecular visualization software

### CIF Structure Files (.cif)
- Experimental structures from PDB
- Used for validation and overlay analysis
- Provide experimental validation of predicted structures

## Usage

### Viewing PyMOL Sessions
1. Open .pse files in PyMOL to view interactive 3D structures
2. Blue spheres represent CRs, red spheres represent SRs
3. Can overlay experimental structures (.cif files) for validation

### Viewing Static Images
1. PNG files show the same information as PyMOL sessions
2. Used as figures in the manuscript
3. Provide quick reference for structural analysis
