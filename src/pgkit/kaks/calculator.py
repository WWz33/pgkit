"""Ka/Ks calculator wrappers and fallback estimators."""

import os
import shutil
import subprocess

from pgkit.core.utils import log
from pgkit.kaks.alignment import align_proteins, back_translate, generate_axt
from pgkit.kaks.constants import CODON_TABLES

# Ka/Ks estimation (Python fallback - Nei-Gojobori)
# ============================================================
def estimate_kaks_ng(cds1, cds2, genetic_code=1):
    """Estimate Ka/Ks using Nei-Gojobori method (simplified)"""
    if len(cds1) != len(cds2) or len(cds1) % 3 != 0:
        return None, None, None

    codon_table = CODON_TABLES.get(genetic_code, CODON_TABLES[1])
    n_codons = len(cds1) // 3

    S_sites, N_sites = 0.0, 0.0
    S_diff, N_diff = 0.0, 0.0

    for i in range(n_codons):
        pos = i * 3
        c1 = cds1[pos:pos+3].upper()
        c2 = cds2[pos:pos+3].upper()

        if '-' in c1 or '-' in c2:
            continue
        if c1 not in codon_table or c2 not in codon_table:
            continue

        aa1, aa2 = codon_table[c1], codon_table[c2]
        if aa1 == '*' or aa2 == '*':
            continue

        # Site counting (simplified)
        S_sites += 1.0
        N_sites += 2.0

        diff_positions = [j for j in range(3) if c1[j] != c2[j]]
        n_diff = len(diff_positions)

        if n_diff == 0:
            continue
        elif n_diff == 1:
            if aa1 == aa2:
                S_diff += 1.0
            else:
                N_diff += 1.0
        else:
            for j in diff_positions:
                test_c1 = c1[:j] + c2[j] + c1[j+1:]
                if test_c1 in codon_table:
                    if codon_table[test_c1] == aa1:
                        S_diff += 0.5
                    else:
                        N_diff += 0.5
                else:
                    N_diff += 0.5

    if S_sites == 0 or N_sites == 0:
        return None, None, None

    pS = S_diff / S_sites
    pN = N_diff / N_sites

    # Jukes-Cantor correction
    Ka = -0.75 * (1.0 - (4.0/3.0) * pN) if 0 < pN < 0.75 else (0.0 if pN == 0 else None)
    Ks = -0.75 * (1.0 - (4.0/3.0) * pS) if 0 < pS < 0.75 else (0.0 if pS == 0 else None)

    if Ka is not None and Ka < 0: Ka = 0.0
    if Ks is not None and Ks < 0: Ks = 0.0

    if Ka is not None and Ks is not None and Ks > 0:
        kaks = Ka / Ks
    elif Ka is not None and Ka > 0 and (Ks is None or Ks == 0):
        kaks = None
    elif Ka == 0 and Ks == 0:
        kaks = 0.0
    else:
        kaks = None

    return Ka, Ks, kaks


# ============================================================
# KaKs_Calculator wrapper
# ============================================================
def run_kakscalculator(axt_file, output_file, method='MA', genetic_code=1,
                       calculator_path=None):
    """
    Run KaKs_Calculator (v2.0 or v3.0)

    Parameters:
        axt_file: Input AXT file
        output_file: Output file
        method: Calculation method (NG, YN, MYN, MA, MS, etc.)
        genetic_code: Genetic code table (1-33)
        calculator_path: Path to KaKs_Calculator executable
    """
    # Expand ~ in path
    if calculator_path:
        calculator_path = os.path.expanduser(calculator_path)

    if calculator_path is None:
        # Try common paths
        for path in ['KaKs_Calculator', 'KaKs', 'KaKs.exe']:
            if shutil.which(path):
                calculator_path = path
                break

        if calculator_path is None:
            # Try local paths
            for path in [
                './KaKs_Calculator3.0/bin/KaKs.exe',
                './KaKs_Calculator2.0/bin/KaKs_Calculator',
                'KaKs_Calculator3.0/bin/KaKs.exe',
            ]:
                if os.path.exists(path):
                    calculator_path = path
                    break

    if calculator_path is None:
        return False, "KaKs_Calculator not found"

    cmd = [
        calculator_path,
        '-i', axt_file,
        '-o', output_file,
        '-m', method,
        '-c', str(genetic_code),
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        # Debug: log command and output if failed
        if result.returncode != 0:
            log(f"KaKs_Calculator failed: {' '.join(cmd)}")
            log(f"  stderr: {result.stderr[:200]}")
        return result.returncode == 0, result.stderr
    except (FileNotFoundError, subprocess.TimeoutExpired) as e:
        return False, str(e)


# ============================================================
# Parse KaKs_Calculator output (supports both v2.0 and v3.0)
# ============================================================
def parse_kakscalculator_output(output_file):
    """
    Parse KaKs_Calculator output file

    v2.0 format (9 columns):
        Sequence  Method  Ka  Ks  Ka/Ks  Length  S-Sites  N-Sites  Fold-Sites(0:2:4)

    v3.0 format (21 columns):
        Sequence  Method  Ka  Ks  Ka/Ks  P-Value(Fisher)  Length  S-Sites  N-Sites
        Fold-Sites(0:2:4)  Substitutions  Syn-Subs  Nonsyn-Subs
        Fold-Syn-Subs(0:2:4)  Fold-Nonsyn-Subs(0:2:4)  Divergence-Distance
        Substitution-Rate-Ratio  GC(1:2:3)  ML-Score  AICc  Akaike-Weight  Model
    """
    results = []

    try:
        with open(output_file, 'r') as f:
            header = f.readline().strip().split('\t')

            # Detect version by number of columns
            n_cols = len(header)
            is_v3 = n_cols >= 15

            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 5:
                    continue

                entry = {
                    'Sequence': parts[0],
                    'Method': parts[1],
                    'Ka': _safe_float(parts[2]),
                    'Ks': _safe_float(parts[3]),
                    'Ka_Ks': _safe_float(parts[4]),
                    'Length': _safe_int(parts[6]) if len(parts) > 6 else None,
                    'S_Sites': _safe_float(parts[7]) if len(parts) > 7 else None,
                    'N_Sites': _safe_float(parts[8]) if len(parts) > 8 else None,
                }

                # v3.0 additional fields
                if is_v3:
                    entry['P_Value'] = _safe_float(parts[5]) if len(parts) > 5 else None
                    entry['Substitutions'] = _safe_float(parts[10]) if len(parts) > 10 else None
                    entry['Syn_Subs'] = _safe_float(parts[11]) if len(parts) > 11 else None
                    entry['Nonsyn_Subs'] = _safe_float(parts[12]) if len(parts) > 12 else None
                    entry['Divergence'] = _safe_float(parts[15]) if len(parts) > 15 else None
                    entry['GC_Content'] = parts[17] if len(parts) > 17 else None
                    entry['ML_Score'] = _safe_float(parts[18]) if len(parts) > 18 else None
                    entry['AICc'] = _safe_float(parts[19]) if len(parts) > 19 else None
                    entry['Akaike_Weight'] = _safe_float(parts[20]) if len(parts) > 20 else None
                    entry['Model'] = parts[21] if len(parts) > 21 else None

                results.append(entry)
    except Exception as e:
        pass

    return results


def _safe_float(s):
    """Safely convert string to float"""
    try:
        return float(s) if s != '-' and s != 'NA' else None
    except (ValueError, IndexError):
        return None


def _safe_int(s):
    """Safely convert string to int"""
    try:
        return int(s) if s != '-' and s != 'NA' else None
    except (ValueError, IndexError):
        return None


# ============================================================
# Single pair Ka/Ks calculation (for multiprocessing)
# ============================================================
def calculate_single_pair(args_tuple):
    """
    Calculate Ka/Ks for a single gene pair

    Used by multiprocessing Pool
    """
    (og_id, g1, g2, protein_dict, cds_dict, tmp_dir,
     use_calculator, method, genetic_code, calculator_path) = args_tuple

    if g1 not in protein_dict or g2 not in protein_dict:
        return None
    if g1 not in cds_dict or g2 not in cds_dict:
        return None

    # Create unique tmp dir for this pair
    pair_tmp = os.path.join(tmp_dir, f"{g1}_{g2}")
    os.makedirs(pair_tmp, exist_ok=True)

    try:
        # Align proteins
        aligned_prot = align_proteins(g1, g2, protein_dict[g1], protein_dict[g2], pair_tmp)

        if aligned_prot is None or len(aligned_prot) < 2:
            return None

        # Back-translate
        cds_pair = {g1: cds_dict[g1], g2: cds_dict[g2]}
        aligned_cds = back_translate(aligned_prot, cds_pair)

        if g1 not in aligned_cds or g2 not in aligned_cds:
            return None

        cds1, cds2 = aligned_cds[g1], aligned_cds[g2]

        if use_calculator:
            # Use KaKs_Calculator
            axt_file = os.path.join(pair_tmp, f"{g1}_{g2}.axt")
            kaks_file = os.path.join(pair_tmp, f"{g1}_{g2}.kaks")

            generate_axt(g1, g2, cds1, cds2, axt_file)

            success, msg = run_kakscalculator(axt_file, kaks_file, method,
                                               genetic_code, calculator_path)

            if success and os.path.exists(kaks_file):
                parsed = parse_kakscalculator_output(kaks_file)
                if parsed:
                    r = parsed[0]
                    return {
                        'Orthogroup': og_id,
                        'Gene1': g1,
                        'Gene2': g2,
                        'Ka': r['Ka'],
                        'Ks': r['Ks'],
                        'Ka_Ks': r['Ka_Ks'],
                        'Method': r.get('Method', method),
                        'P_Value': r.get('P_Value'),
                        'AICc': r.get('AICc'),
                        'Akaike_Weight': r.get('Akaike_Weight'),
                        'Model': r.get('Model'),
                    }

            # Calculator failed - return None (no fallback)
            return None
        else:
            # Fallback: Python Nei-Gojobori
            Ka, Ks, kaks = estimate_kaks_ng(cds1, cds2, genetic_code)

            # Debug: log if calculation failed
            if Ka is None and Ks is None and kaks is None:
                pass  # Don't log every failure

            return {
                'Orthogroup': og_id,
                'Gene1': g1,
                'Gene2': g2,
                'Ka': Ka,
                'Ks': Ks,
                'Ka_Ks': kaks,
                'Method': 'NG-Python',
                'P_Value': None,
                'AICc': None,
                'Akaike_Weight': None,
                'Model': None,
            }

    finally:
        # Cleanup
        if os.path.exists(pair_tmp):
            shutil.rmtree(pair_tmp, ignore_errors=True)


def _fmt(v):
    """Format value for output"""
    if v is None:
        return 'NA'
    if isinstance(v, float):
        return f"{v:.6f}"
    return str(v)
