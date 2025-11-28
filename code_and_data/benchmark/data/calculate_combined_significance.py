#!/usr/bin/env python3
"""
Calculate combined mutational intolerance across multiple deep mutational scanning experiments.
This script:
1. Reads data from B2AR, MC4R, GPR68, and V2R experiments
2. Maps to HRH2 positions
3. Normalizes/ranks values (reversing sensitivity data)
4. Combines scores to get cumulative significance
"""

import csv
import statistics

def read_tsv(filepath):
    """Read TSV file and return list of dicts."""
    with open(filepath, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        return list(reader)

def read_excel_xls(filepath):
    """Read Excel file using xlrd for old .xls format."""
    import xlrd
    wb = xlrd.open_workbook(filepath)
    ws = wb.sheet_by_index(0)
    headers = [ws.cell_value(0, col) for col in range(ws.ncols)]
    data = []
    for row in range(1, ws.nrows):
        row_data = {headers[col]: ws.cell_value(row, col) for col in range(ws.ncols)}
        data.append(row_data)
    return data

def parse_num(value):
    """Parse numeric value, return None if invalid."""
    if value is None or value == '' or value == '-':
        return None
    try:
        return float(value)
    except (ValueError, TypeError):
        return None

def rank_normalize(values_dict):
    """
    Rank normalize a dictionary of {key: value} to 0-1 scale.
    Higher rank = higher significance (more intolerant to mutation).
    Returns dict with same keys.
    """
    # Filter valid values
    valid_items = [(k, v) for k, v in values_dict.items() if v is not None]
    if len(valid_items) == 0:
        return {k: None for k in values_dict.keys()}
    
    # Sort by value and assign ranks
    sorted_items = sorted(valid_items, key=lambda x: x[1])
    ranks = {}
    for i, (k, v) in enumerate(sorted_items):
        ranks[k] = i
    
    # Normalize to 0-1
    max_rank = len(sorted_items) - 1
    if max_rank == 0:
        normalized = {k: 0.5 for k in ranks.keys()}
    else:
        normalized = {k: r / max_rank for k, r in ranks.items()}
    
    # Fill in None for missing values
    result = {k: normalized.get(k, None) for k in values_dict.keys()}
    return result

# Read consolidated mapping file
print("Reading consolidated GPCR mapping file...")
consolidated_map = read_tsv('../Consolidated_GPCR_mapping.tsv')

# Build mapping dictionaries for each receptor -> HRH2
adrb2_map_dict = {}
mc4r_map_dict = {}
gpr68_map_dict = {}
v2r_map_dict = {}

for row in consolidated_map:
    hrh2_pos = parse_num(row.get('HRH2_resNum'))
    if hrh2_pos:
        hrh2_pos = int(hrh2_pos)
        
        # ADRB2 mapping
        adrb2_pos = parse_num(row.get('ADRB2_resNum'))
        if adrb2_pos:
            adrb2_map_dict[int(adrb2_pos)] = hrh2_pos
        
        # MC4R mapping
        mc4r_pos = parse_num(row.get('MC4R_resNum'))
        if mc4r_pos:
            mc4r_map_dict[int(mc4r_pos)] = hrh2_pos
        
        # GPR68 mapping
        gpr68_pos = parse_num(row.get('GPR68_resNum'))
        if gpr68_pos:
            gpr68_map_dict[int(gpr68_pos)] = hrh2_pos
        
        # V2R mapping
        v2r_pos = parse_num(row.get('V2R_resNum'))
        if v2r_pos:
            v2r_map_dict[int(v2r_pos)] = hrh2_pos

print(f"  ADRB2->HRH2 mappings: {len(adrb2_map_dict)}")
print(f"  MC4R->HRH2 mappings: {len(mc4r_map_dict)}")
print(f"  GPR68->HRH2 mappings: {len(gpr68_map_dict)}")
print(f"  V2R->HRH2 mappings: {len(v2r_map_dict)}")

# Read HRH2 labels file for reference positions
print("\nReading HRH2 reference data...")
hrh2_data = read_tsv('HRH2_entropy_conservation.txt')
hrh2_positions = sorted(set(int(row['HRH2_res_no']) for row in hrh2_data if row.get('HRH2_res_no')))

# Initialize results dictionary
results = {pos: {
    'HRH2_resNum': pos,
    'B2AR_score': None,
    'MC4R_score': None,
    'GPR68_score': None,
    'V2R_score': None,
    'Emax_score': None,
    'BArr_score': None,
    'Tolerance': None,
    'sens_score': None,
    'pH5.5': None,
    'expression_score': None,
    'Emax': None,
    'BArr_norm': None
} for pos in hrh2_positions}

# ==================== B2AR Tolerance ====================
print("\nProcessing B2AR tolerance data...")
try:
    # Read B2AR tolerance data from Excel
    adrb2_tol = read_excel_xls('elife-54895-supp3-v2 (1).xls')
    print(f"  Total rows in B2AR file: {len(adrb2_tol)}")
    
    # Filter for condition 0.625
    adrb2_tol_filtered = []
    for row in adrb2_tol:
        cond = parse_num(row.get('Condition'))
        if cond is not None and abs(cond - 0.625) < 1e-6:
            adrb2_tol_filtered.append(row)
    print(f"  Rows with Condition=0.625: {len(adrb2_tol_filtered)}")
    
    # Use consolidated mapping (already loaded above)
    # adrb2_map_dict is already available
    
    # Build tolerance dict by HRH2 position
    tolerance_by_hrh2 = {}
    for row in adrb2_tol_filtered:
        adrb2_pos = parse_num(row.get('Pos'))
        tolerance = parse_num(row.get('Tolerance'))
        if adrb2_pos and tolerance is not None:
            hrh2_pos = adrb2_map_dict.get(int(adrb2_pos))
            if hrh2_pos and hrh2_pos in results:
                tolerance_by_hrh2[hrh2_pos] = tolerance
                results[hrh2_pos]['Tolerance'] = tolerance
    
    # For tolerance: LOWER values = MORE sensitive to mutation = HIGHER significance
    # So we REVERSE the ranking (1 - rank_normalize)
    normalized = rank_normalize(tolerance_by_hrh2)
    for hrh2_pos, score in normalized.items():
        if score is not None:
            results[hrh2_pos]['B2AR_score'] = 1 - score
    
    print(f"  Found {len(tolerance_by_hrh2)} HRH2 positions with B2AR data")
except Exception as e:
    print(f"  Warning: Could not process B2AR data: {e}")
    import traceback
    traceback.print_exc()

# ==================== MC4R Sensitivity ====================
print("\nProcessing MC4R sensitivity data...")
try:
    mc4r_sens = read_tsv('mc4r_residue_sensitivity_table.tsv')
    
    # Use consolidated mapping (already loaded above)
    # mc4r_map_dict is already available
    
    # Build sensitivity dict by HRH2 position
    sens_by_hrh2 = {}
    for row in mc4r_sens:
        mc4r_pos = parse_num(row.get('pos'))
        sens = parse_num(row.get('sens_score'))
        if mc4r_pos and sens is not None:
            hrh2_pos = mc4r_map_dict.get(int(mc4r_pos))
            if hrh2_pos and hrh2_pos in results:
                sens_by_hrh2[hrh2_pos] = sens
                results[hrh2_pos]['sens_score'] = sens
    
    # For sensitivity: HIGHER values = MORE sensitive to mutation = HIGHER significance
    normalized = rank_normalize(sens_by_hrh2)
    for hrh2_pos, score in normalized.items():
        if score is not None:
            results[hrh2_pos]['MC4R_score'] = score
    
    print(f"  Found {len(sens_by_hrh2)} HRH2 positions with MC4R data")
except Exception as e:
    print(f"  Warning: Could not process MC4R data: {e}")

# ==================== GPR68 pH 5.5 ====================
print("\nProcessing GPR68 pH 5.5 data...")
try:
    gpr68_data = read_tsv('GPR68_normalized_activation_scores.tsv')
    
    # Use consolidated mapping (already loaded above)
    # gpr68_map_dict is already available
    
    # Build pH 5.5 dict by HRH2 position
    ph55_by_hrh2 = {}
    for row in gpr68_data:
        gpr68_pos = parse_num(row.get('pos'))
        ph55 = parse_num(row.get('pH5.5'))
        if gpr68_pos and ph55 is not None:
            hrh2_pos = gpr68_map_dict.get(int(gpr68_pos))
            if hrh2_pos and hrh2_pos in results:
                ph55_by_hrh2[hrh2_pos] = ph55
                results[hrh2_pos]['pH5.5'] = ph55
    
    # For pH 5.5: LOWER values = loss of function = HIGHER significance
    # So we REVERSE the ranking
    normalized = rank_normalize(ph55_by_hrh2)
    for hrh2_pos, score in normalized.items():
        if score is not None:
            results[hrh2_pos]['GPR68_score'] = 1 - score
    
    print(f"  Found {len(ph55_by_hrh2)} HRH2 positions with GPR68 data")
except Exception as e:
    print(f"  Warning: Could not process GPR68 data: {e}")

# ==================== V2R Expression ====================
print("\nProcessing V2R expression data...")
try:
    v2r_expr = read_tsv('average_V2R_expression_per_residue.tsv')
    
    # Use consolidated mapping (already loaded above)
    # v2r_map_dict is already available
    
    # Build expression dict by HRH2 position
    expr_by_hrh2 = {}
    for row in v2r_expr:
        v2r_pos = None
        expr_score = None
        for key, val in row.items():
            if 'residue' in key.lower() and 'number' in key.lower():
                v2r_pos = parse_num(val)
            if key.lower() == 'score':
                expr_score = parse_num(val)
        
        if v2r_pos and expr_score is not None:
            hrh2_pos = v2r_map_dict.get(int(v2r_pos))
            if hrh2_pos and hrh2_pos in results:
                expr_by_hrh2[hrh2_pos] = expr_score
                results[hrh2_pos]['expression_score'] = expr_score
    
    # For expression: LOWER values = loss of expression = HIGHER significance
    # So we REVERSE the ranking
    normalized = rank_normalize(expr_by_hrh2)
    for hrh2_pos, score in normalized.items():
        if score is not None:
            results[hrh2_pos]['V2R_score'] = 1 - score
    
    print(f"  Found {len(expr_by_hrh2)} HRH2 positions with V2R data")
except Exception as e:
    print(f"  Warning: Could not process V2R data: {e}")

# ==================== B2AR Emax (G protein) ====================
print("\nProcessing B2AR Emax (G protein) data...")
try:
    emax_data = read_tsv('Emax_mean_by_position_QCpassed.tsv')
    
    # Read mapping (ADRB2 to HRH2) - reuse from B2AR tolerance
    # adrb2_map_dict already exists from above
    
    # Build Emax dict by HRH2 position
    emax_by_hrh2 = {}
    for row in emax_data:
        adrb2_pos = parse_num(row.get('pos'))
        emax = parse_num(row.get('Emax'))
        if adrb2_pos and emax is not None:
            hrh2_pos = adrb2_map_dict.get(int(adrb2_pos))
            if hrh2_pos and hrh2_pos in results:
                emax_by_hrh2[hrh2_pos] = emax
                results[hrh2_pos]['Emax'] = emax
    
    # For Emax: LOWER values = reduced efficacy = HIGHER significance
    # So we REVERSE the ranking
    normalized = rank_normalize(emax_by_hrh2)
    for hrh2_pos, score in normalized.items():
        if score is not None:
            results[hrh2_pos]['Emax_score'] = 1 - score
    
    print(f"  Found {len(emax_by_hrh2)} HRH2 positions with Emax data")
except Exception as e:
    print(f"  Warning: Could not process Emax data: {e}")

# ==================== B2AR Beta-Arrestin (Normalized) ====================
print("\nProcessing B2AR Beta-Arrestin (normalized) data...")
try:
    barr_data = read_tsv('BArr_EC100_mean_by_position.tsv')
    
    # Read mapping (ADRB2 to HRH2) - reuse from B2AR tolerance
    # adrb2_map_dict already exists from above
    
    # Build Beta-Arrestin dict by HRH2 position
    barr_by_hrh2 = {}
    for row in barr_data:
        adrb2_pos = parse_num(row.get('pos'))
        barr_norm = parse_num(row.get('BArr_EC100_norm_mean'))
        if adrb2_pos and barr_norm is not None:
            hrh2_pos = adrb2_map_dict.get(int(adrb2_pos))
            if hrh2_pos and hrh2_pos in results:
                barr_by_hrh2[hrh2_pos] = barr_norm
                results[hrh2_pos]['BArr_norm'] = barr_norm
    
    # For Beta-Arrestin: LOWER values = reduced recruitment = HIGHER significance
    # So we REVERSE the ranking
    normalized = rank_normalize(barr_by_hrh2)
    for hrh2_pos, score in normalized.items():
        if score is not None:
            results[hrh2_pos]['BArr_score'] = 1 - score
    
    print(f"  Found {len(barr_by_hrh2)} HRH2 positions with Beta-Arrestin data")
except Exception as e:
    print(f"  Warning: Could not process Beta-Arrestin data: {e}")

# ==================== Calculate Combined Score ====================
print("\nCalculating combined significance scores...")

score_cols = ['B2AR_score', 'MC4R_score', 'GPR68_score', 'V2R_score', 'Emax_score', 'BArr_score']

for pos, data in results.items():
    # Count how many datasets have this position
    data['n_datasets'] = sum(1 for col in score_cols if data[col] is not None)
    
    # Fill missing with 0 (position doesn't exist/align in that receptor)
    # Always divide by 6 to reward positions present across multiple datasets
    scores = [data[col] if data[col] is not None else 0.0 for col in score_cols]
    data['combined_score'] = sum(scores) / 6.0
    data['combined_significance'] = data['combined_score']

# Save results
output_file = 'combined_mutational_intolerance.tsv'
with open(output_file, 'w', newline='', encoding='utf-8') as f:
    fieldnames = ['HRH2_resNum', 'n_datasets', 'combined_significance', 'combined_score',
                  'B2AR_score', 'MC4R_score', 'GPR68_score', 'V2R_score', 'Emax_score', 'BArr_score',
                  'Tolerance', 'sens_score', 'pH5.5', 'expression_score', 'Emax', 'BArr_norm']
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    for pos in sorted(results.keys()):
        writer.writerow(results[pos])

print(f"Results saved to: {output_file}")

# Filter results (≥2 datasets) for plotting
results_filtered = {pos: data for pos, data in results.items() if data['n_datasets'] >= 2}
output_filtered = 'combined_mutational_intolerance_filtered.tsv'
with open(output_filtered, 'w', newline='', encoding='utf-8') as f:
    fieldnames = ['HRH2_resNum', 'n_datasets', 'combined_significance', 'combined_score',
                  'B2AR_score', 'MC4R_score', 'GPR68_score', 'V2R_score', 'Emax_score', 'BArr_score']
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()
    for pos in sorted(results_filtered.keys()):
        row = {k: results_filtered[pos][k] for k in fieldnames}
        writer.writerow(row)

print(f"Filtered results saved to: {output_filtered}")

# Summary statistics
total_positions = len(results)
positions_2plus = len(results_filtered)
positions_all4 = sum(1 for data in results.values() if data['n_datasets'] == 4)

print(f"\nSummary:")
print(f"  Total HRH2 positions analyzed: {total_positions}")
print(f"  Positions with ≥2 datasets: {positions_2plus}")
print(f"  Positions with all 4 datasets: {positions_all4}")

# Print top 20 most significant positions
print("\nTop 20 most mutational intolerant positions (highest combined significance):")
sorted_positions = sorted(
    [(pos, data) for pos, data in results_filtered.items()],
    key=lambda x: x[1]['combined_significance'] if x[1]['combined_significance'] is not None else -1,
    reverse=True
)[:20]

print(f"{'HRH2':<8} {'Score':<8} {'N':<4} {'B2AR':<8} {'MC4R':<8} {'GPR68':<8} {'V2R':<8}")
print("-" * 60)
for pos, data in sorted_positions:
    print(f"{pos:<8} {data['combined_significance']:.4f}   {data['n_datasets']:<4} "
          f"{data['B2AR_score'] if data['B2AR_score'] is not None else 'N/A':<8} "
          f"{data['MC4R_score'] if data['MC4R_score'] is not None else 'N/A':<8} "
          f"{data['GPR68_score'] if data['GPR68_score'] is not None else 'N/A':<8} "
          f"{data['V2R_score'] if data['V2R_score'] is not None else 'N/A':<8}")

# Statistics
all_scores = [data['combined_significance'] for data in results_filtered.values() 
              if data['combined_significance'] is not None]
if all_scores:
    print(f"\nCombined significance statistics:")
    print(f"  Mean: {statistics.mean(all_scores):.4f}")
    print(f"  Median: {statistics.median(all_scores):.4f}")
    print(f"  Min: {min(all_scores):.4f}")
    print(f"  Max: {max(all_scores):.4f}")
    print(f"  StdDev: {statistics.stdev(all_scores):.4f}")
