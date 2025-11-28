# -*- coding: utf-8 -*-
"""
Created on Wed Nov 12 10:34:36 2025

@author: selcuk.1
"""

import pandas as pd
import re
from pathlib import Path

# ========= YOU EDIT THESE =========
GENE_TSV = r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\GPCR_classes_manuscript\Supplementary Information\code_and_data\scatter_plots_and_CR_SR_calculation\uniprot_hgcn_mapping.txt"                # columns: uniprot, hgnc
CLINVAR_SUMMARY = r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\GPCR_classes_manuscript\Supplementary Information\code_and_data\scatter_plots_and_CR_SR_calculation\variant_summary.txt" # .txt or .txt.gz
OUT_FILE = r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\GPCR_classes_manuscript\Supplementary Information\code_and_data\scatter_plots_and_CR_SR_calculation\variant_summary_classA_GPCR.txt"    # .txt or .txt.gz
CHUNK_SIZE = 200_000                                 # adjust if needed
# ==================================

def _normalize_hgnc_value(x: str):
    """Return ('ID','HGNC:####') or ('SYMBOL','SYMBOL') or None."""
    if x is None:
        return None
    s = str(x).strip()
    if not s:
        return None
    s_up = s.upper()
    if re.fullmatch(r"HGNC:\d+", s_up):
        return ("ID", s_up)
    if re.fullmatch(r"\d+", s_up):
        return ("ID", f"HGNC:{s_up}")
    return ("SYMBOL", s_up)

def _load_gene_list(tsv_path: str):
    df = pd.read_csv(tsv_path, sep=None, engine="python", dtype=str)  # auto-detect sep
    if "hgnc" not in {c.lower() for c in df.columns}:
        raise ValueError("Your TSV must contain a column named 'hgnc' (case-insensitive).")
    hcol = [c for c in df.columns if c.lower() == "hgnc"][0]
    ids, symbols = set(), set()
    for v in df[hcol].dropna().astype(str):
        norm = _normalize_hgnc_value(v)
        if not norm:
            continue
        kind, val = norm
        if kind == "ID":
            ids.add(val)
        else:
            symbols.add(val)
    return ids, symbols

def filter_variant_summary_by_hgnc(
    gene_tsv: str,
    variant_summary_path: str,
    out_path: str,
    chunksize: int = 200_000
):
    hgnc_ids, hgnc_symbols = _load_gene_list(gene_tsv)
    if not hgnc_ids and not hgnc_symbols:
        raise ValueError("No usable HGNC values parsed from your TSV.")

    # Prepare output file (ensure dir exists)
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)

    # We read header first to keep same column order
    head_df = pd.read_csv(variant_summary_path, sep="\t", nrows=1, dtype=str, compression="infer")
    required_cols = ["HGNC_ID", "GeneSymbol"]
    for col in required_cols:
        if col not in head_df.columns:
            raise ValueError(f"Column '{col}' not found in variant_summary file.")
    # write header to OUT (overwrite if exists)
    head_df.iloc[0:0].to_csv(out_path, sep="\t", index=False, compression="infer")

    kept_rows = 0
    matched_ids, matched_syms = set(), set()

    # Stream in chunks
    for chunk in pd.read_csv(
        variant_summary_path,
        sep="\t",
        dtype=str,
        chunksize=chunksize,
        compression="infer",
        keep_default_na=False
    ):
        # Normalize columns used for matching
        chunk["HGNC_ID_norm"] = chunk["HGNC_ID"].str.upper()
        chunk["GeneSymbol_norm"] = chunk["GeneSymbol"].str.upper()

        mask = False
        if hgnc_ids:
            mask = chunk["HGNC_ID_norm"].isin(hgnc_ids)
        if hgnc_symbols:
            mask = mask | chunk["GeneSymbol_norm"].isin(hgnc_symbols)

        sub = chunk.loc[mask, head_df.columns]  # keep original column order

        if not sub.empty:
            # Append without header
            sub.to_csv(out_path, sep="\t", index=False, mode="a", header=False, compression="infer")
            kept_rows += len(sub)

            # Track which requested genes we actually saw
            matched_ids.update(set(sub["HGNC_ID"].str.upper()) & hgnc_ids)
            matched_syms.update(set(sub["GeneSymbol"].str.upper()) & hgnc_symbols)

    present_gene_count = len(matched_ids | matched_syms)
    requested_gene_count = len(hgnc_ids) + len(hgnc_symbols)

    print(f"[clinvar-filter] Output rows written: {kept_rows}")
    print(f"[clinvar-filter] Genes requested: {requested_gene_count}")
    print(f"[clinvar-filter] Genes present in subset: {present_gene_count}")

    # Optional: show missing counts (comment out if noisy)
    missing_ids = hgnc_ids - matched_ids
    missing_syms = hgnc_symbols - matched_syms
    print(f"[clinvar-filter] Missing HGNC IDs: {len(missing_ids)}")
    print(f"[clinvar-filter] Missing symbols: {len(missing_syms)}")

# ===== run it =====
filter_variant_summary_by_hgnc(GENE_TSV, CLINVAR_SUMMARY, OUT_FILE, CHUNK_SIZE)

#%%
import pandas as pd
import re
from pathlib import Path

# =========================
# -------- FUNCTIONS -------
# =========================

def _label_from_clinsig(clinsig: str) -> str | None:
    """
    Map ClinicalSignificance text to a binary label:
      - 'pathogenic' if Pathogenic or Likely_pathogenic (incl combos)
      - 'benign'      if Benign or Likely_benign (incl combos)
      - None          otherwise (uncertain/conflicting/risk/association/drug_response/protective/affects)
    """
    if not isinstance(clinsig, str):
        return None
    s = clinsig.lower().replace(" ", "_")
    if any(tok in s for tok in ["uncertain", "conflict", "risk", "association", "drug_response", "protective", "affects"]):
        return None
    has_path = "pathogenic" in s
    has_ben  = "benign" in s
    if has_path and not has_ben:
        return "pathogenic"
    if has_ben and not has_path:
        return "benign"
    return None  # mixed benign/pathogenic → ambiguous → drop

# --- missense detection helpers ---

# true missense pattern like p.Arg131Gly / p.R131G (optionally wrapped in parentheses)
_MISSENSE_RE = re.compile(r"^p\.?\(?([A-Za-z]{1,3})(\d+)([A-Za-z]{1,3})\)?$", flags=re.IGNORECASE)
_STOP_TOKENS = {"*", "TER", "X"}
AA3 = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G',
    'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
    'THR':'T','TRP':'W','TYR':'Y','VAL':'V','SEC':'U','PYL':'O'
}

def _norm_aa(tok: str):
    t = tok.upper()
    if len(t) == 1: return t
    return AA3.get(t, None)

def _is_synonymous(ps: str) -> bool:
    """Return True for p.Ser23= / p.(Ser23=) / Ser23= styles."""
    if not isinstance(ps, str): 
        return False
    s = ps.strip()
    if s.startswith("p."): s = s[2:]
    s = s.strip("()").strip()
    return s.endswith("=")

def _is_true_missense(ps: str) -> bool:
    """True AA→AA substitution (no fs/*/= or del/ins/dup/ext)."""
    if not isinstance(ps, str) or not ps.strip():
        return False
    ps2 = ps.strip()
    if ps2.startswith("p."): ps2 = ps2[2:]
    ps2 = ps2.strip("()").strip()
    low = ps2.lower()
    if any(x in low for x in ["fs", "del", "ins", "dup", "ext"]):  # indels/frameshift/extension
        return False
    if ps2.endswith("="):  # synonymous
        return False
    m = _MISSENSE_RE.match("p." + ps2 if not ps2.startswith("p.") else ps2)
    if not m:
        return False
    ref, _, alt = m.groups()
    ref1 = _norm_aa(ref); alt1 = _norm_aa(alt)
    if ref1 is None or alt1 is None: return False
    if alt1 in _STOP_TOKENS: return False  # nonsense like * / Ter / X
    return ref1 != alt1

# residue extractor (works for missense/nonsense/fs/etc.)
_AA_POS_RE = re.compile(r"p\.[A-Za-z]{1,3}(\d+)[A-Za-z\*]{0,20}", flags=re.IGNORECASE)

def _extract_protein_sub(row) -> str | None:
    """
    Prefer 'Protein_change'; else parse from 'Name' like:
    'NM_...(...):c.... (p.Arg131Gly)'
    """
    val = row.get("Protein_change", "")
    if isinstance(val, str) and val.strip().startswith("p."):
        return val.strip()
    name = row.get("Name", "")
    if isinstance(name, str):
        m = re.search(r"\(p\.[^)]+\)", name)
        if m:
            return m.group(0).strip("()")
    return None

def _extract_residue_number(protein_sub: str) -> int | None:
    if not isinstance(protein_sub, str):
        return None
    m = _AA_POS_RE.search(protein_sub)
    return int(m.group(1)) if m else None

def process_clinvar_subset_keep_missense(
    in_path: str,
    out_path: str,
    chunksize: int = 200_000
):
    """
    Start from your existing ClinVar **subfile** (already gene-filtered).
    Keep ONLY missense rows with clear benign/likely_benign or
    pathogenic/likely_pathogenic, and write:

      GeneSymbol, HGNC_ID, ProteinSub, res_num, label, ClinicalSignificance
    """
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)

    head = pd.read_csv(in_path, sep="\t", nrows=1, dtype=str, compression="infer")
    required = ["GeneSymbol", "HGNC_ID", "ClinicalSignificance"]
    for c in required:
        if c not in head.columns:
            raise ValueError(f"Required column '{c}' not found in input.")
    has_protein_change = "Protein_change" in head.columns
    has_name = "Name" in head.columns
    if not (has_protein_change or has_name):
        raise ValueError("Neither 'Protein_change' nor 'Name' is present; cannot derive ProteinSub.")
    has_mc = "MolecularConsequence" in head.columns

    out_cols = ["GeneSymbol", "HGNC_ID", "ProteinSub", "res_num", "label", "ClinicalSignificance"]
    pd.DataFrame(columns=out_cols).to_csv(out_path, sep="\t", index=False)

    kept_rows = 0
    seen_genes = set()
    label_counts = {"benign": 0, "pathogenic": 0}

    for chunk in pd.read_csv(
        in_path, sep="\t", dtype=str, chunksize=chunksize, compression="infer", keep_default_na=False
    ):
        # pull columns we use
        need = ["GeneSymbol", "HGNC_ID", "ClinicalSignificance"]
        if has_protein_change: need.append("Protein_change")
        if has_name: need.append("Name")
        if has_mc: need.append("MolecularConsequence")
        sub = chunk[need].copy()

        # derive ProteinSub
        sub["ProteinSub"] = sub.apply(_extract_protein_sub, axis=1)

        # --- robust missense filter ---
        # If ProteinSub present → must be true missense (and not synonymous)
        ps_present = sub["ProteinSub"].notna() & (sub["ProteinSub"].astype(str).str.len() > 0)
        ps_true_miss = sub["ProteinSub"].apply(_is_true_missense)

        # If ProteinSub missing → allow MC to admit missense rows
        if has_mc:
            mc_has_missense = chunk["MolecularConsequence"].astype(str).str.lower().str.contains("missense_variant")
        else:
            mc_has_missense = pd.Series(False, index=sub.index)

        missense_mask = (ps_present & ps_true_miss) | (~ps_present & mc_has_missense)
        sub = sub.loc[missense_mask].copy()
        if sub.empty:
            continue

        # hard guard: drop any stray '=' that slipped through for safety
        bad_syn = sub["ProteinSub"].astype(str).str.contains(r"=", regex=True, na=False)
        if bad_syn.any():
            sub = sub.loc[~bad_syn].copy()
            if sub.empty():
                continue

        # label filter
        sub["label"] = sub["ClinicalSignificance"].apply(_label_from_clinsig)
        sub = sub.loc[sub["label"].notna()].copy()
        if sub.empty:
            continue

        # residue number
        sub["res_num"] = sub["ProteinSub"].apply(_extract_residue_number)
        sub["res_num"] = pd.to_numeric(sub["res_num"], errors="coerce").astype("Int64")
        sub = sub.loc[sub["res_num"].notna()].copy()
        if sub.empty:
            continue

        # write skinny
        sub = sub[["GeneSymbol", "HGNC_ID", "ProteinSub", "res_num", "label", "ClinicalSignificance"]]
        sub.to_csv(out_path, sep="\t", index=False, mode="a", header=False)

        kept_rows += len(sub)
        seen_genes.update(sub["GeneSymbol"].unique())
        for k, v in sub["label"].value_counts().to_dict().items():
            label_counts[k] = label_counts.get(k, 0) + v

    print(f"[missense-only] Wrote {kept_rows} rows to: {out_path}")
    print(f"[missense-only] Unique genes present: {len(seen_genes)}")
    print(f"[missense-only] Label counts → benign: {label_counts.get('benign',0)}, pathogenic: {label_counts.get('pathogenic',0)}")

# =========================
# --------- CALLS ---------
# =========================

# ====== YOU EDIT THESE ======
IN_SUBFILE  = r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\GPCR_classes_manuscript\Supplementary Information\code_and_data\scatter_plots_and_CR_SR_calculation\variant_summary_classA_GPCR.txt"
OUT_TSV     = r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\GPCR_classes_manuscript\Supplementary Information\code_and_data\scatter_plots_and_CR_SR_calculation\variant_summary_classA_GPCR_missense_simplified.txt"
CHUNK_SIZE  = 200_000
# ==========================

process_clinvar_subset_keep_missense(IN_SUBFILE, OUT_TSV, CHUNK_SIZE)
#%%

import pandas as pd
from pathlib import Path
import re

# =========================
# -------- FUNCTIONS -------
# =========================

def _norm_hgnc_key(x: str):
    """
    Normalize a mapping 'hgnc' cell:
      - 'HGNC:1234' -> ('ID','HGNC:1234')
      - '1234'      -> ('ID','HGNC:1234')
      - 'ADRB2'     -> ('SYMBOL','ADRB2')
    """
    if x is None:
        return (None, None)
    s = str(x).strip()
    if not s:
        return (None, None)
    u = s.upper()
    if re.fullmatch(r"HGNC:\d+", u):
        return ("ID", u)
    if re.fullmatch(r"\d+", u):
        return ("ID", f"HGNC:{u}")
    return ("SYMBOL", u)

def _prep_mapping_df(mapping_path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Load simple mapping TSV (columns: uniprot, hgnc) and return two lookup tables:
      - map_id (HGNC_ID_norm -> UniProt)
      - map_sym (GeneSymbol_norm -> UniProt)
    If multiple UniProt per key, they are semicolon-joined.
    """
    m = pd.read_csv(mapping_path, sep=None, engine="python", dtype=str)
    cols = {c.lower(): c for c in m.columns}
    if "uniprot" not in cols or "hgnc" not in cols:
        raise ValueError("Mapping file must have columns: 'uniprot' and 'hgnc' (case-insensitive).")
    ucol, hcol = cols["uniprot"], cols["hgnc"]

    # normalize keys
    kinds, keys = zip(*m[hcol].map(_norm_hgnc_key))
    m["__kind"] = kinds
    m["__key"]  = keys
    m["UniProt"] = m[ucol].astype(str).str.strip()

    # by HGNC ID
    map_id = (m[m["__kind"]=="ID"]
                .dropna(subset=["__key"])
                .groupby("__key")["UniProt"]
                .apply(lambda s: ";".join(sorted(set([x for x in s if x]))))
                .reset_index()
                .rename(columns={"__key":"HGNC_ID_norm"}))

    # by Gene symbol
    map_sym = (m[m["__kind"]=="SYMBOL"]
                .dropna(subset=["__key"])
                .groupby("__key")["UniProt"]
                .apply(lambda s: ";".join(sorted(set([x for x in s if x]))))
                .reset_index()
                .rename(columns={"__key":"GeneSymbol_norm"}))

    return map_id, map_sym

def add_uniprot_column(
    in_tsv: str,
    mapping_tsv: str,
    out_tsv: str
):
    """
    Read your simplified missense TSV (must have GeneSymbol, HGNC_ID),
    add a 'UniProt' column using mapping TSV (columns: uniprot, hgnc),
    and write to out_tsv.
    """
    Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(in_tsv, sep="\t", dtype=str)
    need = {"GeneSymbol","HGNC_ID"}
    if not need.issubset(df.columns):
        raise ValueError(f"Input file must contain columns: {need}")

    # prep mapping
    map_id, map_sym = _prep_mapping_df(mapping_tsv)

    # normalize keys in df for joining
    df["HGNC_ID_norm"] = df["HGNC_ID"].astype(str).str.upper().str.replace(r"^\s*(\d+)\s*$", r"HGNC:\1", regex=True)
    df["GeneSymbol_norm"] = df["GeneSymbol"].astype(str).str.upper()

    # join by HGNC_ID first
    out = df.merge(map_id, on="HGNC_ID_norm", how="left")

    # fill missing UniProt via GeneSymbol
    if "UniProt" in out.columns:
        missing = out["UniProt"].isna() | (out["UniProt"].astype(str).str.len()==0)
    else:
        missing = pd.Series(True, index=out.index)

    if missing.any():
        out2 = out.loc[missing].merge(map_sym, on="GeneSymbol_norm", how="left", suffixes=("","_sym"))
        out.loc[missing, "UniProt"] = out2["UniProt"]

    # tidy columns: put UniProt first
    preferred = ["UniProt","GeneSymbol","HGNC_ID"]
    keep_existing = [c for c in out.columns if c not in preferred and c not in ["HGNC_ID_norm","GeneSymbol_norm"]]
    out = out[preferred + keep_existing]

    out.to_csv(out_tsv, sep="\t", index=False)
    matched = out["UniProt"].notna().sum()
    print(f"[uniprot] Added UniProt for {matched} rows. Wrote: {out_tsv}")

# =========================
# --------- CALLS ---------
# =========================

# >>> EDIT THESE PATHS <<<
IN_SIMPLIFIED = r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\GPCR_classes_manuscript\Supplementary Information\code_and_data\scatter_plots_and_CR_SR_calculation\variant_summary_classA_GPCR_missense_simplified.txt"
MAP_TSV       = r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\GPCR_classes_manuscript\Supplementary Information\code_and_data\scatter_plots_and_CR_SR_calculation\uniprot_hgcn_mapping.txt"                             # columns: uniprot, hgnc
OUT_WITH_UP   = r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\GPCR_classes_manuscript\Supplementary Information\code_and_data\scatter_plots_and_CR_SR_calculation\variant_summary_classA_GPCR_missense_simplified_with_uniprot.txt"

add_uniprot_column(IN_SIMPLIFIED, MAP_TSV, OUT_WITH_UP)

#%%
import pandas as pd
from pathlib import Path

# =========================
# -------- FUNCTIONS -------
# =========================

def load_variants(variant_tsv: str) -> dict:
    """
    Read simplified missense TSV with columns: UniProt, res_num, label (benign/deleterious).
    Returns dict UniProt -> list of (res_num:int, sign:+1/-1).
    """
    df = pd.read_csv(variant_tsv, sep="\t", dtype=str)
    need = {"UniProt", "res_num", "label"}
    if not need.issubset(df.columns):
        raise ValueError(f"Variant file must contain columns {need}")
    # normalize
    df["UniProt"] = df["UniProt"].astype(str).str.strip()
    df["res_num"] = pd.to_numeric(df["res_num"], errors="coerce").astype("Int64")
    df = df[df["res_num"].notna()].copy()
    df["label"] = df["label"].str.lower().str.strip()
    sign_map = {"pathogenic": 1, "benign": -1}
    df = df[df["label"].isin(sign_map.keys())].copy()
    df["sign"] = df["label"].map(sign_map)

    by_up = {}
    for up, sub in df.groupby("UniProt"):
        by_up[up] = list(zip(sub["res_num"].astype(int).tolist(), sub["sign"].tolist()))
    return by_up

def parse_fasta_with_uniprot_third_field(msa_fasta: str):
    """
    Read FASTA MSA. Header is split by '|', and header_parts[2] is UniProt accession.
    Returns:
      seqs: list of dicts: {"header":str, "uniprot":str, "seq":str}
      up_index: dict UniProt -> index in seqs
    """
    seqs = []
    with open(msa_fasta, "r", encoding="utf-8") as fh:
        header = None
        chunks = []
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    seqs.append({"header": header, "uniprot": _uniprot_from_header(header), "seq": "".join(chunks)})
                header = line.strip()[1:]
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            seqs.append({"header": header, "uniprot": _uniprot_from_header(header), "seq": "".join(chunks)})

    up_index = {}
    for i, rec in enumerate(seqs):
        up = rec["uniprot"]
        if up:  # first wins if duplicates
            up_index.setdefault(up, i)
    return seqs, up_index

def _uniprot_from_header(header: str) -> str | None:
    """
    Split header by '|' and return 3rd field as UniProt accession if present.
    """
    parts = header.split("|")
    if len(parts) >= 3:
        return parts[2].strip()
    return None

def residue_to_column_map(seq: str) -> dict:
    """
    For a single aligned sequence, return mapping: residue_index(1-based) -> alignment_column(0-based).
    Gaps ('-') do not advance the residue index.
    """
    r = 0
    mapping = {}
    for col, aa in enumerate(seq):
        if aa != "-":
            r += 1
            mapping[r] = col
    return mapping

def find_ref_col_to_resnum(seqs: list, ref_name_substring: str = "HRH2_HUMAN") -> dict:
    """
    Find the sequence whose header contains ref_name_substring (case-insensitive),
    and build a mapping: alignment_column(0-based) -> HRH2_HUMAN residue number (1-based),
    with None where HRH2_HUMAN has a gap.
    """
    ref_idx = None
    for i, rec in enumerate(seqs):
        if ref_name_substring.lower() in rec["header"].lower():
            ref_idx = i
            break
    if ref_idx is None:
        raise ValueError(f"Could not find reference sequence containing '{ref_name_substring}' in headers.")

    seq = seqs[ref_idx]["seq"]
    col_to_res = {}
    res = 0
    for col, aa in enumerate(seq):
        if aa == "-":
            col_to_res[col] = None
        else:
            res += 1
            col_to_res[col] = res
    return col_to_res

def aggregate_alignment_scores_no_evidence_filtered(msa_fasta: str, variant_tsv: str, ref_name_substring: str, out_tsv: str):
    """
    - Reads MSA (3rd header field = UniProt accession) and variant table (UniProt, res_num, label).
    - Maps variants to alignment columns and sums +1/-1 per column.
    - Converts alignment columns to HRH2_HUMAN residue numbers.
    - Writes TSV with only residues that have >=1 mapped variant (no-evidence columns are dropped).
      Keeps residues even if total_score == 0 (as long as there was evidence).
    Output columns: HRH2_HUMAN_resnum, total_score, n_evidence
    """
    # Load
    variants = load_variants(variant_tsv)  # {UniProt: [(res_num, +1/-1), ...]}
    seqs, up_index = parse_fasta_with_uniprot_third_field(msa_fasta)

    # Precompute residue->column map for needed sequences
    res2col_by_up = {}
    for up, idx in up_index.items():
        if up in variants:
            res2col_by_up[up] = residue_to_column_map(seqs[idx]["seq"])

    # Per-column totals and evidence counts
    n_cols = len(seqs[0]["seq"]) if seqs else 0
    col_sum = [0] * n_cols
    col_evidence = [0] * n_cols  # number of variant rows mapped here

    missing_map = 0
    used_pairs = 0
    for up, items in variants.items():
        rmap = res2col_by_up.get(up)
        if rmap is None:
            continue
        for resnum, sign in items:
            col = rmap.get(resnum)
            if col is None:
                missing_map += 1
                continue
            col_sum[col] += int(sign)
            col_evidence[col] += 1
            used_pairs += 1

    # Map columns → HRH2_HUMAN residue numbers
    col_to_hrh2 = find_ref_col_to_resnum(seqs, ref_name_substring=ref_name_substring)

    # Build rows: include only columns with HRH2 residue AND evidence>0
    rows = []
    for col, score in enumerate(col_sum):
        res = col_to_hrh2.get(col)
        if res is not None and col_evidence[col] > 0:
            rows.append((res, score, col_evidence[col]))

    # Write
    from pathlib import Path
    import pandas as pd
    Path(out_tsv).parent.mkdir(parents=True, exist_ok=True)
    out_df = pd.DataFrame(rows, columns=["HRH2_HUMAN_resnum", "total_score", "n_evidence"]).sort_values("HRH2_HUMAN_resnum")
    out_df.to_csv(out_tsv, sep="\t", index=False)

    print(f"[msa] Variants mapped to columns: {used_pairs}")
    print(f"[msa] Variants with no column (gaps/isoform mismatch): {missing_map}")
    print(f"[msa] Residues written (evidence > 0): {len(out_df)} → {out_tsv}")

# =========================
# --------- CALLS ---------
# =========================

# >>> EDIT THESE PATHS <<<
MSA_FASTA   = r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\GPCRevolution\public\alignments\classA_humans_MSA.fasta" # FASTA; headers like '>sp|Q9...|P12345|...'
VAR_TSV     = r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\GPCR_classes_manuscript\Supplementary Information\code_and_data\scatter_plots_and_CR_SR_calculation\benchmark\variant_summary_classA_GPCR_missense_simplified_with_uniprot.txt"  # from previous step (has UniProt, res_num, label)
REF_NAME    = "HRH2_HUMAN"   # substring to identify the reference sequence in FASTA headers
OUT_TSV     = r"C:\Users\selcuk.1\OneDrive - The Ohio State University\Desktop\GPCR_classes_manuscript\Supplementary Information\code_and_data\scatter_plots_and_CR_SR_calculation\benchmark\ClinVar_aggregated_pathogenicity.txt"

aggregate_alignment_scores_no_evidence_filtered(MSA_FASTA, VAR_TSV, REF_NAME, OUT_TSV)
