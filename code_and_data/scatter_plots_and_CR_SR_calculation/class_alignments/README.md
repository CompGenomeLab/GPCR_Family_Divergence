# Class Alignments - GPCR Multiple Sequence Alignments - Decoding Functional Specialization in GPCRs through Evolution-Guided Residue Profiling

This directory contains multiple sequence alignments (MSAs) for each GPCR class, used for family-wide conservation analysis and receptor mapping, supporting the manuscript "Decoding Functional Specialization in GPCRs through Evolution-Guided Residue Profiling". These alignments are the foundation for identifying Common Residues (CRs) and Selective Residues (SRs) across GPCR families.

## Alignment Generation Process

### Representative Sequence Selection
The alignments are generated using a multi-step process:

1. **CD-HIT Clustering**: Ortholog sequences are clustered using CD-HIT with 65% identity threshold
2. **Top 5 Representatives**: The 5 largest clusters are selected as representative sequences ("reps")
3. **Human Sequence Addition**: Human sequences are added to the representative set
4. **MAFFT Alignment**: Multiple sequence alignment is performed using MAFFT
5. **Human Extraction**: Human sequences are extracted from the final alignment for family-wide analysis

### Script: `get_human_and_reps.py`

This script orchestrates the alignment generation process:

**Key Functions:**
- `get_reps()`: Extracts representative sequences from CD-HIT clusters
- `human_seq_get()`: Retrieves human sequences for each receptor
- `dict_to_fasta()`: Converts sequence dictionaries to FASTA format
- `get_cluster_sizes()`: Sorts clusters by size to select top representatives

**Input Requirements:**
- CD-HIT cluster files (`.clstr` format)
- Receptor lists for each GPCR class
- Human sequence files from ortholog pipeline

**Output:**
- FASTA files containing human sequences + representative sequences
- Files are organized by GPCR class and alignment type

## File Naming Conventions

### Alignment Types

#### 1. **Full Alignments** (with representatives)
- `classA_top5_reps_whuman.fasta` - Class A with representatives + humans
- `classB_top5_reps_whuman.fasta` - Class B with representatives + humans
- `classC_top5_reps_whuman.fasta` - Class C with representatives + humans
- `classF_top5_reps_whuman.fasta` - Class F with representatives + humans
- `classT_top5_reps_whuman.fasta` - Class T with representatives + humans
- `classOlf_top5_reps_whuman.fasta` - Olfactory with representatives + humans

#### 2. **Refined Alignments** (MAFFT L-INS-i)
- `classA_top5_reps_whuman_linsi.fasta` - Class A refined with L-INS-i
- `classB_top5_reps_whuman_linsi.fasta` - Class B refined with L-INS-i
- `classC_top5_reps_whuman_linsi.fasta` - Class C refined with L-INS-i
- `classF_top5_reps_whuman_linsi.fasta` - Class F refined with L-INS-i
- `classT_top5_reps_whuman_linsi.fasta` - Class T refined with L-INS-i
- `classOlf_top5_reps_whuman_linsi.fasta` - Olfactory refined with L-INS-i

#### 3. **Human-Only Alignments** (extracted for analysis)
- `classA_top5_reps_whuman_linsi_onlyhuman.fasta` - Class A human sequences only
- `classB_top5_reps_whuman_linsi_onlyhuman.fasta` - Class B human sequences only
- `classC_top5_reps_whuman_linsi_onlyhuman.fasta` - Class C human sequences only
- `classF_top5_reps_whuman_linsi_onlyhuman.fasta` - Class F human sequences only
- `classT_top5_reps_whuman_linsi_onlyhuman.fasta` - Class T human sequences only
- `classOlf_top5_reps_whuman_linsi_onlyhuman.fasta` - Olfactory human sequences only

#### 4. **Domain-Specific Alignments**
- `classC_top5_reps_whuman_linsi_onlyhumanVFT.fasta` - Class C VFT domain only
- `classC_top5_reps_whuman_linsi_onlyhumanVFTandCRD.fasta` - Class C VFT + CRD domains

#### 5. **Class-Specific Subsets**
- `classB1_humans_MSA.fasta` - Class B1 receptors only
- `classB2_humans_MSA.fasta` - Class B2 receptors only
- `classB_top5_reps_whuman_linsi_onlyhuman_onlyclassB1.fasta` - Class B1 subset
- `classB_top5_reps_whuman_linsi_onlyhuman_onlyclassBadh.fasta` - Class B2 subset

## File Sizes and Content

### Large Files (Full Alignments)
- `classA_top5_reps_whuman_linsi.fasta` (44MB) - Complete Class A alignment
- `classOlf_top5_reps_whuman_linsi.fasta` (19MB) - Complete Olfactory alignment
- `classB_top5_reps_whuman_linsi.fasta` (5.0MB) - Complete Class B alignment
- `classC_top5_reps_whuman_linsi.fasta` (912KB) - Complete Class C alignment

### Medium Files (Human-Only Alignments)
- `classA_top5_reps_whuman_linsi_onlyhuman.fasta` (1.4MB) - Class A humans only
- `classOlf_top5_reps_whuman_linsi_onlyhuman.fasta` (577KB) - Olfactory humans only
- `classB_top5_reps_whuman_linsi_onlyhuman.fasta` (417KB) - Class B humans only

### Small Files (Domain-Specific)
- `classF_top5_reps_whuman_linsi_onlyhuman.fasta` (13KB) - Class F humans only
- `classT_top5_reps_whuman_linsi_onlyhuman.fasta` (10KB) - Class T humans only
- `classC_top5_reps_whuman_linsi_onlyhumanVFT.fasta` (50KB) - Class C VFT domain
- `classC_top5_reps_whuman_linsi_onlyhumanVFTandCRD.fasta` (43KB) - Class C VFT+CRD

## Usage in Analysis

### 1. **Family-Wide Conservation Analysis**
- Human-only alignments are used for calculating family-wide conservation percentages
- Each position in the alignment represents a homologous residue across family members
- Conservation is calculated as the percentage of receptors with conserved residues at each position

### 2. **Entropy Calculations**
- Human-only alignments provide the receptor set for entropy calculations
- Shannon entropy is calculated for each position across the family
- High entropy positions indicate selective residues (SRs)
- Low entropy positions indicate common residues (CRs)

### 3. **Receptor Mapping**
- Alignment positions determine residue numbering and conservation calculations
- Different alignments produce different outcomes due to varying receptor sets
- Domain-specific alignments focus on particular functional regions

### 4. **Structural Analysis**
- Alignments are used to map conservation patterns onto 3D structures
- CR/SR classifications are based on alignment-derived conservation metrics
- Structural visualizations use alignment-derived residue positions

## Impact of Alignment Choice

### Receptor Set Variation
- Different representative sequences change receptor composition
- Class-specific subsets (B1 vs B2) provide different evolutionary perspectives
- Domain-specific alignments focus on particular functional regions

### Conservation Calculations
- Alignment positions determine residue numbering
- Conservation percentages depend on the receptor set included
- Entropy values vary based on sequence diversity in the alignment

### Analysis Scope
- Full alignments provide complete evolutionary context
- Human-only alignments enable family-wide analysis
- Domain-specific alignments reveal region-specific conservation patterns

## Technical Details

### Sequence Headers
- Human sequences: `sp|UniProtID|GeneName_HUMAN|9606`
- Representative sequences: Include cluster information and source species
- Taxonomic ID 9606 indicates human sequences

### Alignment Quality
- MAFFT L-INS-i provides high-quality alignments
- Gap removal and trimming applied for analysis
- Quality control ensures reliable conservation calculations

### File Formats
- All files in FASTA format
- Aligned sequences with gaps represented as "-"
- Headers contain protein identification and taxonomic information

## Dependencies

- **CD-HIT**: For sequence clustering and representative selection
- **MAFFT**: For multiple sequence alignment (FFT-NS and L-INS-i algorithms)
- **Python**: For script execution and data processing
- **BioPython**: For sequence manipulation and analysis

This alignment strategy ensures robust family-wide analysis while maintaining computational efficiency through representative sequence selection. 