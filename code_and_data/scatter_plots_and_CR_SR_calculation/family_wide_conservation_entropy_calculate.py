# -*- coding: utf-8 -*-
"""
Created on Thu May 15 14:24:15 2025

@author: selcuk.1
"""

import os
import json
import pandas as pd
import math
from collections import Counter

# --- helper functions ---

def load_metadata(json_path: str) -> dict:
    """Load JSON metadata and return a mapping geneName -> metadata dict."""
    with open(json_path) as jf:
        meta_list = json.load(jf)
    return {entry['geneName']: entry for entry in meta_list}



def filter_and_rename_genes(raw_genes: list, meta_by_gene: dict) -> (list, list):
    """Filter genes by ortholog count >= 26"""
    included = [g for g in raw_genes if meta_by_gene.get(g, {}).get('numOrthologs', 0) >= 26]
    excluded = [g for g in raw_genes if g not in included]
    return included, excluded


def parse_aa_conservation_file(file_path: str) -> dict:
    """Parse an AAconservation.txt file into pos(int) -> (aa_str, label) map."""
    pos_map = {}
    if not os.path.exists(file_path):
        return pos_map
    with open(file_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 4 or not parts[0].isdigit():
                continue
            pos = int(parts[0])
            aa_str = parts[2]
            label = parts[3]
            pos_map[pos] = (aa_str, label)
    return pos_map


def calculate_conservation_entropy(cons_list: list, aa_list: list) -> tuple:
    """
    Given a list of 1/0 conservation flags and corresponding amino acids, return:
      - percent conserved (float)
      - entropy (float)
    """
    n = len(cons_list)
    conserved_count = sum(cons_list)
    pct_conserved = (conserved_count / n * 100) if n > 0 else 0.0
    freqs = Counter(aa_list)
    entropy = 0.0
    for aa, count in freqs.items():
        p = count / n
        if p > 0:
            entropy -= p * math.log2(p)
    return round(pct_conserved, 2), round(entropy, 2)

from Bio import SeqIO

def build_alignment_map(alignment_fasta: str, reference: str) -> dict:
    """
    Parse a FASTA MSA and build mapping from each seq’s residue positions
    to the reference sequence positions.

    Parameters
    ----------
    alignment_fasta : str
        FASTA file containing reference + other receptors.
    reference : str
        ID of the reference in FASTA (must match a header-derived gene name).

    Returns
    -------
    dict
        { gene_name: { seq_pos: ref_pos, … }, … }
    """
    # read all sequences into dict { gene_name: sequence }
    seqs = {}
    for rec in SeqIO.parse(alignment_fasta, "fasta"):

        gene_name = rec.id.split("|")[2].split("_")[0]  
        seqs[gene_name] = str(rec.seq)

    if reference not in seqs:
        raise KeyError(f"Reference '{reference}' not found in alignment headers")

    ref_seq = seqs[reference]
    alignment_map = {}
    # build mapping for each sequence
    for gene, seq in seqs.items():
        ref_pos = 0
        seq_pos = 0
        mapping = {}
        for r_aa, s_aa in zip(ref_seq, seq):
            if r_aa != "-":
                ref_pos += 1
            if s_aa != "-":
                seq_pos += 1
            if r_aa != "-" and s_aa != "-":
                mapping[ref_pos] = seq_pos
        alignment_map[gene] = mapping

    return alignment_map

# --- main analysis function ---

def analyze_conservation_labels(alignment_map,
                                 json_path: str,
                                 base_dir: str,
                                 date: str,
                                 class_name: str,
                                 reference_gene: str,
                                 output_tsv_path: str):
    """
    Orchestrates conservation/entropy analysis across receptors.
    """
    # Load and filter genes
    meta_by_gene = load_metadata(json_path)

    raw_genes = list(alignment_map.keys())
    included, excluded = filter_and_rename_genes(raw_genes, meta_by_gene)
    # Extract reference residue order
    ref_map = alignment_map.get(reference_gene)
    if ref_map is None:
        raise ValueError(f"Reference gene '{reference_gene}' not found in alignment_map")
    reference_positions = list(ref_map.keys())
    n_positions = len(reference_positions)

    # Prepare master matrices keyed by alignment index
    cons_mat = {i: [] for i in range(n_positions)}
    aa_mat = {i: [] for i in range(n_positions)}
    gap_C_cases = []
    aa_list_all = []

    # Iterate each receptor, aligned by reference_positions
    for gene in included:
        # Build per-reference mapping list using direct gene key
        gene_map = alignment_map[gene]
        mapping_list = [gene_map.get(ref_res, '-') for ref_res in reference_positions]

        if gene=="CML2":
            # Load conservation labels for this gene
            file_path = os.path.join(base_dir,
                                     f"GPR1_{date}_{class_name}",
                                     f"GPR1_{date}_AAconservation.txt")
        else:
            # Load conservation labels for this gene
            file_path = os.path.join(base_dir,
                                     f"{gene}_{date}_{class_name}",
                                     f"{gene}_{date}_AAconservation.txt")
        pos_map = parse_aa_conservation_file(file_path)
        # Aggregate data at each reference index
        for idx, resno in enumerate(mapping_list):
            if resno!="-":
                aa_str, label = pos_map.get(resno)
            else:
                aa_str, label = ('-', 'NC')

            first_aa = aa_str.split('/')[0]
            cons_mat[idx].append(1 if label == 'C' else 0)
            aa_mat[idx].append(first_aa)
            if label == 'C' and first_aa == '-':
                gap_C_cases.append((gene, idx))
            if first_aa not in aa_list_all:
                aa_list_all.append(first_aa)

    # Summarize per reference position
    results = []
    for idx, ref_res in enumerate(reference_positions):
        cons_flags = cons_mat[idx]
        aas = aa_mat[idx]
        pct, ent = calculate_conservation_entropy(cons_flags, aas)
        results.append({
            'Residue_Map': ref_res,
            'Percent_Conserved': pct,
            'Entropy': ent
        })

    out_df = pd.DataFrame(results)
    out_df.to_csv(output_tsv_path, sep='\t', index=False)

    return out_df, {'included_genes': included, 'excluded_genes': excluded, "gap_conservation_error":gap_C_cases}
#Dir for the pipeline files
base_dir=r"PATH_For_Pipeline_Files" #You can use the zenodo link to download
#Path for the receptors.json file. Used to filter based on number of orthologs.
jsonf="receptors.json"
#Date of the pipeline run, necessary for folder name match
date="26-10-2022"
#%%
reference_gene="HRH2"
class_name="classA"
class_fasta=r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\GPCR_classes_manuscript\Supplementary Information\code_and_data\generation_of_scatter_plots\class_alignments\classA_top5_reps_whuman_linsi_onlyhuman.fasta"
output_tsv_path=r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\classA_HRH2_conservation_and_entropy.tsv"
#%%
reference_gene="O52I2"
class_name="classA"
class_fasta=r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\GPCR_classes_manuscript\Supplementary Information\code_and_data\generation_of_scatter_plots\class_alignments\classOlf_top5_reps_whuman_linsi_onlyhuman.fasta"
output_tsv_path=r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\classOlf_O52I2_conservation_and_entropy.tsv"
#%%
reference_gene="T2R39"
class_name="classT"
class_fasta=r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\GPCR_classes_manuscript\Supplementary Information\code_and_data\generation_of_scatter_plots\class_alignments\classT_top5_reps_whuman_linsi_onlyhuman.fasta"
output_tsv_path=r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\classT_T2R39_conservation_and_entropy.tsv"
#%%
reference_gene="PTH1R"
class_name="classB"
class_fasta=r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\GPCR_classes_manuscript\Supplementary Information\code_and_data\generation_of_scatter_plots\class_alignments\classB_top5_reps_whuman_linsi_onlyhuman_onlyclassB1.fasta"
output_tsv_path=r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\classB1_PTH1R_conservation_and_entropy.tsv"
#%%
reference_gene="AGRL3"
class_name="classB"
class_fasta=r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\GPCR_classes_manuscript\Supplementary Information\code_and_data\generation_of_scatter_plots\class_alignments\classB_top5_reps_whuman_linsi_onlyhuman_onlyclassBadh.fasta"
output_tsv_path=r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\classBadh_AGRL3_conservation_and_entropy.tsv"
#%%
reference_gene="FZD7"
class_name="classF"
class_fasta=r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\GPCR_classes_manuscript\Supplementary Information\code_and_data\generation_of_scatter_plots\class_alignments\classF_top5_reps_whuman_linsi_onlyhuman.fasta"
output_tsv_path=r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\classF_FZD7_conservation_and_entropy.tsv"
#%%
reference_gene="CASR"
class_name="classC"
class_fasta=r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\GPCR_classes_manuscript\Supplementary Information\code_and_data\generation_of_scatter_plots\class_alignments\classC_top5_reps_whuman_linsi_onlyhumanVFT.fasta"
output_tsv_path=r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\classC_CASR_conservation_and_entropy.tsv"
#%% For the CRD analysis
reference_gene="CASR"
class_name="classC"
class_fasta=r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\GPCR_classes_manuscript\Supplementary Information\code_and_data\scatter_plots_and_CR_SR_calculation\class_alignments\classC_top5_reps_whuman_linsi_onlyhumanVFTandCRD.fasta"
output_tsv_path=r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\classC_CASR_conservation_and_entropy_VFTandCRD.tsv"
#%%
alignment_map=build_alignment_map(class_fasta, reference_gene)
df, summary = analyze_conservation_labels(
    alignment_map,
    jsonf,
    base_dir,
    date,
    class_name,
    reference_gene,
    output_tsv_path
)

