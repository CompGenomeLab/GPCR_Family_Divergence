# -*- coding: utf-8 -*-
"""
Created on Mon Apr 21 14:03:03 2025

@author: selcuk.1
"""

import pandas as pd
import re

def filter_rows_by_conservation(
    tsv_path: str,
    conservation_threshold,
    min_fraction
) -> pd.DataFrame:
    """
    Read a TSV with alternating 'Residue'/'Conservation' columns,
    parse out the numeric conservation percentages, and return only
    those rows where at least `min_fraction` of the conservation
    columns are >= `conservation_threshold`.

    Parameters
    ----------
    tsv_path : str
        Path to the input TSV file.
    conservation_threshold : float, optional
        The conservation percentage threshold (default 90.0).
    min_fraction : float, optional
        Minimum fraction (0–1) of conservation columns that must
        meet/exceed `conservation_threshold` (default 0.5).

    Returns
    -------
    pd.DataFrame
        Subset of the original DataFrame meeting the criterion.
    """
    # 1) Load all columns as strings
    df = pd.read_csv(tsv_path, sep='\t', dtype=str)

    # 2) Identify the conservation columns
    cons_cols = [c for c in df.columns if c.startswith('Conservation')]

    # 3) Parse out the numeric percentages (others → NaN)
    def parse_pct(s):
        if isinstance(s, str):
            m = re.search(r'(\d+(?:\.\d+)?)\s*%', s)
            if m:
                return float(m.group(1))
        return float('nan')

    pct_df = df[cons_cols].applymap(parse_pct)

    # 4) Compute for each row the fraction of columns ≥ threshold
    meets = pct_df.ge(conservation_threshold)
    fraction = meets.sum(axis=1) / len(cons_cols)

    # 5) Tag & filter
    df['_high_cons_frac'] = fraction
    print(df)
    return df[df['_high_cons_frac'] >= min_fraction].drop(columns=['_high_cons_frac'])

from collections import Counter
import blosum as bl

def filter_identical_and_similar_conservation(
    df: "pd.DataFrame",
    blosum_matrix=80,
    score_threshold: int = 1
) -> "pd.DataFrame":
    """
    Keep only those rows where the conserved AAs across receptors are either
      1) all identical, or
      2) all similar to the most frequent one (BLOSUM score > score_threshold).

    Parameters
    ----------
    df : pd.DataFrame
        Output of filter_rows_by_conservation (with your Conservation cols).
    blosum_matrix : int or str, optional
        If int, loads BLOSUM{blosum_matrix} (45, 50, 62, 80, 90).
        If str, treats as a file path to load a custom matrix.
        Default is 80.
    score_threshold : int, optional
        Minimum BLOSUM score to consider “similar” (default 1).
    """
    # load the matrix
    matrix = bl.BLOSUM(blosum_matrix)
    print(matrix)

    # identify conservation columns
    cons_cols = [c for c in df.columns if c.startswith("Conservation")]

    # helper: get first AA letter (None for “--” or missing)
    def first_aa(s):
        aa_part = s.split()[0]        # e.g. "L/M"
        return aa_part.split("/")[0]  # take the first

    aa_df = df[cons_cols].applymap(first_aa)

    keep_idxs = []
    for idx, row in aa_df.iterrows():
        letters = [a for a in row if a is not None]
        if not letters:
            continue

        # most common AA
        modal_aa, _ = Counter(letters).most_common(1)[0]

        # all identical?
        if all(a == modal_aa for a in letters):
            keep_idxs.append(idx)
            continue

        # else check similarity
        def sim(a, b):
            # matrix[a][b] returns the score
            if {a,b}=={"D","E"} or {a,b}=={"W","F"}:  #Exceptions to blosum80 similarity
                return 2
            return matrix[a][b]

        if all(sim(modal_aa, a) > score_threshold for a in letters if a != modal_aa):
            keep_idxs.append(idx)

    return df.loc[keep_idxs].reset_index(drop=True)

#%%
cons_data=r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\GPCR_classes_manuscript\Supplementary Information\code_and_data\CR_detection_classF\classF_conservation_FZD7.tsv"
filtered = filter_rows_by_conservation(
    cons_data,
    conservation_threshold=90.0,
    min_fraction=0.72
)

refined = filter_identical_and_similar_conservation(filtered)
#%%
print(refined["FZD7 Residue"])
refined.to_csv(r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\GPCR_classes_manuscript\Supplementary Information\code_and_data\CR_detection_classF\classF_conservation_FZD7_CRs.tsv",sep="\t")
