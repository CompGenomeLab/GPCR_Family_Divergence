# GPCR Ortholog Pipeline - Decoding Functional Specialization in GPCRs through Evolution-Guided Residue Profiling

This directory contains the complete pipeline for identifying and analyzing GPCR orthologs across different classes (A, B, C, F, T), supporting the manuscript "Decoding Functional Specialization in GPCRs through Evolution-Guided Residue Profiling". The pipeline consists of Python scripts in the `Scripts/` folder and shell scripts for SLURM cluster execution.

## Directory Structure

```text
ortholog_pipeline/
├── Scripts/                    # Python scripts for data processing
├── human_sequences/           # Human protein sequences by GPCR class
│   ├── classA_fasta/         # Class A GPCR human sequences
│   ├── classB_fasta/         # Class B GPCR human sequences
│   ├── classC_fasta/         # Class C GPCR human sequences
│   ├── classF_fasta/         # Class F GPCR human sequences
│   └── classT_fasta/         # Class T GPCR human sequences
├── Pipeline_Part1.sh         # Initial BLAST and gene clade identification
├── Pipline_Part2.sh          # Refined phylogenetic analysis
├── Pipeline_Part3_orthologMSA_get.sh  # Final alignment and quality control
└── README.md                 # This documentation file
```

## Human Sequences Folder

The `human_sequences/` folder contains human protein sequences organized by GPCR class:

- **classA_fasta/**: Contains 730 human Class A/Olfactory GPCR sequences (e.g., 5HT1A_HUMAN.fasta, DRD2_HUMAN.fasta)
- **classB_fasta/**: Contains 48 human Class B GPCR sequences (e.g., AGRA1_HUMAN.fasta)
- **classC_fasta/**: Contains 22 human Class C GPCR sequences (e.g., CASR_HUMAN.fasta, GABR1_HUMAN.fasta)
- **classF_fasta/**: Contains 11 human Class F GPCR sequences (e.g., FZD1_HUMAN.fasta, SMO_HUMAN.fasta)
- **classT_fasta/**: Contains 25 human Class T GPCR sequences (e.g., T2R10_HUMAN.fasta)

These sequences serve as query proteins for BLAST searches against the UniProt reference proteomes database to identify orthologs across species.

## Shell Scripts Overview

### Pipeline_Part1.sh
**Purpose**: Performs initial BLAST search and gene clade identification for a single protein.

**Key Features**:
- Uses SLURM array jobs to process multiple proteins in parallel
- Contains protein arrays for different GPCR classes (Class A, B, C, F, T)
- Executes BLAST search against UniProt reference proteomes
- Generates initial phylogenetic trees using FastTree
- Identifies gene-specific clades and representative sequences

**Input**: Human protein sequence from specific GPCR class
**Output**: Gene clade tree and representative sequences

### Pipline_Part2.sh
**Purpose**: Performs refined phylogenetic analysis and paralog removal.

**Key Features**:
- Uses MAFFT L-INS-i for high-quality alignment refinement
- Performs model selection using IQ-Tree2 Modelfinder
- Constructs maximum likelihood trees using RAxML-NG
- Removes paralogous sequences using statistical criteria
- Reorders trees to prioritize human sequences

**Input**: Gene clade sequences from Part 1
**Output**: Cleaned ortholog tree and sequences

### Pipeline_Part3_orthologMSA_get.sh
**Purpose**: Extracts final aligned sequences and performs quality control.

**Key Features**:
- Extracts aligned sequences from MSA based on final tree
- Removes columns with all gaps using TrimAl
- Generates final high-quality ortholog alignment

**Input**: Ortholog tree and MSA from Part 2
**Output**: Final ortholog alignment ready for analysis

## Script Descriptions

### 1. FastaObtainerV2.py
**Purpose**: Extracts FASTA sequences from BLAST results and proteome databases.

**Inputs**:
- `--blastout`: BLAST output file in tabular format
- `--num`: Number of target sequences to extract
- `--taxid`: Taxonomic ID to filter sequences
- `--out`: Output FASTA file path

**Functionality**:
- Parses BLAST results to identify proteins from specific taxonomic groups
- Extracts corresponding sequences from proteome databases
- Removes duplicate sequences within the same taxonomic group
- Generates a clean FASTA file with filtered sequences

### 2. gene_clade_find.py
**Purpose**: Identifies gene clades in phylogenetic trees and determines appropriate outgroups.

**Inputs**:
- `--tree`: Phylogenetic tree file in Newick format
- `--blastout`: BLAST result file for detecting target gene and outgroup
- `--out`: Output tree file path

**Functionality**:
- Analyzes phylogenetic trees to identify gene-specific clades
- Determines appropriate outgroup sequences for rooting
- Handles duplications within target species (e.g., human)
- Outputs a pruned tree with identified gene clade

### 3. ete3order.py
**Purpose**: Reorders phylogenetic trees to place human sequences at the top.

**Inputs**:
- `--tree`: Input tree file in Newick format
- `--out`: Output tree file path

**Functionality**:
- Traverses phylogenetic tree nodes
- Identifies human sequences (taxID 9606)
- Swaps child nodes to place human sequences at the top of the tree
- Outputs reordered tree with human sequences prioritized

### 4. ete3trimmer.py
**Purpose**: Removes paralogous sequences from phylogenetic trees using multiple criteria.

**Inputs**:
- `--tree`: Phylogenetic tree file in Newick format
- `--fasta`: FASTA file containing raw sequences
- `--protein`: Protein name of interest
- `--out`: Output tree file path

**Functionality**:
- Compares sequence similarity using BLOSUM62 scoring
- Analyzes taxonomic relationships and evolutionary distances
- Performs statistical tests to identify significant differences between clades
- Removes paralogous sequences based on similarity, distance, and taxonomic criteria
- Prioritizes clades containing human sequences
- Outputs a cleaned ortholog tree

### 5. subtree_alignment.py
**Purpose**: Extracts aligned sequences from a multiple sequence alignment (MSA) based on tree leaf names.

**Inputs**:
- `--tree`: Tree file in Newick format
- `--msa`: Source multiple sequence alignment file
- `--out`: Output FASTA file path

**Functionality**:
- Extracts leaf names from phylogenetic tree
- Searches for corresponding sequences in the MSA
- Creates a new FASTA file with only the sequences present in the tree
- Reports statistics on found vs. missing sequences

### 6. subtree_fasta.py
**Purpose**: Extracts FASTA sequences from proteome databases based on tree leaf names.

**Inputs**:
- `--tree`: Tree file in Newick format
- `--out`: Output FASTA file path

**Functionality**:
- Extracts leaf names from phylogenetic tree
- Searches for corresponding sequences in proteome databases
- Parses taxonomic IDs from sequence headers
- Creates a new FASTA file with sequences matching tree leaves
- Reports statistics on found vs. missing sequences

### 7. iqtree_model_get.py
**Purpose**: Extracts the best substitution model from IQ-Tree2 Modelfinder log files.

**Inputs**:
- `--logfile`: IQ-Tree2 modelfinder log file
- `--modelout`: Output text file for the best model

**Functionality**:
- Parses IQ-Tree2 log files for Bayesian Information Criterion results
- Extracts the best-fitting substitution model
- Outputs the model name to a text file for use in downstream analyses

## Dependencies

The scripts require the following Python packages:
- `ete3`: For phylogenetic tree manipulation
- `Bio`: For sequence alignment and scoring
- `scipy`: For statistical tests
- `numpy`: For numerical operations
- `pickle`: For data serialization

## Usage Notes

1. **File Paths**: Some scripts contain hardcoded database paths that may need to be updated for your system
2. **Taxonomic IDs**: Scripts use NCBI taxonomic IDs (e.g., 9606 for human)
3. **Tree Formats**: Scripts expect Newick format trees with specific node naming conventions
4. **Sequence Headers**: FASTA headers should include taxonomic information in specific formats

## Pipeline Workflow

The complete ortholog identification pipeline consists of three main parts, each executed as separate SLURM jobs:

### Part 1: Initial BLAST and Gene Clade Identification (`Pipeline_Part1.sh`)

1. **BLAST Search**: Performs BLASTp search against UniProt reference proteomes database
   - Query: Human protein sequence from specific GPCR class
   - Output: Tabular BLAST results with taxonomic information

2. **Sequence Extraction**: Uses `FastaObtainerV2.py` to extract sequences from BLAST results
   - Filters by taxonomic ID (9606 for human)
   - Removes duplicates within same taxonomic group
   - Generates initial FASTA file

3. **Multiple Sequence Alignment**: Uses MAFFT for initial alignment
   - Applies FFT-NS algorithm for speed
   - Trims alignment using ClipKit with kpic-gappy method

4. **Fast Tree Construction**: Uses FastTree for initial phylogenetic tree
   - Generates approximate maximum likelihood tree

5. **Gene Clade Identification**: Uses `gene_clade_find.py` to identify gene-specific clades
   - Determines appropriate outgroup sequences
   - Handles duplications within target species
   - Outputs pruned tree with identified gene clade

6. **Representative Sequence Selection**: Uses CD-HIT to reduce redundancy
   - Clusters sequences at 70% and 65% identity thresholds
   - Generates representative sequence sets

### Part 2: Refined Phylogenetic Analysis (`Pipeline_Part2.sh`)

1. **Refined Alignment**: Uses MAFFT L-INS-i for high-quality alignment
   - Applies iterative refinement with up to 1000 iterations
   - Generates improved multiple sequence alignment

2. **Model Selection**: Uses IQ-Tree2 Modelfinder to determine best substitution model
   - Tests various amino acid substitution models
   - Uses `iqtree_model_get.py` to extract best model from log file

3. **Maximum Likelihood Tree**: Uses RAxML-NG for robust tree construction
   - Applies best-fitting substitution model
   - Uses parsimony and random starting trees
   - Generates final phylogenetic tree

4. **Paralog Removal**: Uses `ete3trimmer.py` to remove paralogous sequences
   - Compares sequence similarity using BLOSUM62 scoring
   - Analyzes taxonomic relationships and evolutionary distances
   - Performs statistical tests to identify significant differences
   - Prioritizes clades containing human sequences

5. **Tree Reordering**: Uses `ete3order.py` to place human sequences at top
   - Swaps child nodes to prioritize human sequences (taxID 9606)

6. **Final Sequence Extraction**: Uses `subtree_fasta.py` to extract final ortholog sequences
   - Creates FASTA file with cleaned ortholog set

### Part 3: Final Alignment and Quality Control (`Pipeline_Part3_orthologMSA_get.sh`)

1. **Ortholog Alignment Extraction**: Uses `subtree_alignment.py` to extract aligned sequences
   - Extracts sequences from MSA based on final tree leaves
   - Creates aligned FASTA file with ortholog set

2. **Alignment Quality Control**: Uses TrimAl to remove problematic columns
   - Removes columns with all gaps
   - Generates final high-quality ortholog alignment

### Key Features of the Pipeline:

- **Multi-stage filtering**: BLAST → Gene clade identification → Paralog removal
- **Quality control**: Multiple alignment refinement steps and gap removal
- **Statistical rigor**: Model selection and statistical tests for paralog identification
- **Human-centric**: Prioritizes human sequences and maintains human orthologs
- **Scalable**: Designed for SLURM cluster execution with array jobs
- **Reproducible**: Each step generates intermediate files for validation

### Input Requirements:

- **Human Sequences**: The `human_sequences/` folder provides the starting point for each analysis
  - Each GPCR class has its own subfolder (classA_fasta/, classB_fasta/, etc.)
  - Individual FASTA files contain human protein sequences (e.g., 5HT1A_HUMAN.fasta)
  - These sequences serve as BLAST queries to identify orthologs across species
- **UniProt Reference Proteomes Database**: Comprehensive protein database for BLAST searches
- **Taxonomic Information**: Sequence headers must include taxonomic IDs for proper filtering
- **SLURM Cluster Environment**: Required modules include BLAST+, MAFFT, IQ-Tree2, RAxML-NG, etc.

### Human Sequences Utilization in Pipeline:

1. **Part 1**: Human sequences from `human_sequences/` are used as BLAST queries
   - Script reads from appropriate class folder (e.g., classA_fasta/ for Class A GPCRs)
   - Each protein gets its own analysis directory (e.g., 5HT1A_26-10-2022_classA/)
   - BLAST search identifies orthologs across species using human sequence as query

2. **Part 2**: Human sequences are prioritized during paralog removal
   - `ete3trimmer.py` specifically looks for human sequences (taxID 9606)
   - Clades containing human sequences are preserved during paralog filtering
   - `ete3order.py` places human sequences at the top of phylogenetic trees

3. **Part 3**: Final alignments maintain human sequences as reference
   - Human sequences serve as anchors for alignment quality assessment
   - Final ortholog sets always include the human reference sequence

### Output Files:

- Ortholog phylogenetic trees (Newick format)
- Multiple sequence alignments (FASTA format)
- Representative sequence sets
- Log files with processing statistics
- Best substitution models for each protein family