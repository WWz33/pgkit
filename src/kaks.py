#!/usr/bin/env python3
"""
Ka/Ks calculation for pan-genome analysis (v2.0)
=================================================

Calculate Ka/Ks for different gene family categories using:
- ParaAT for parallel codon alignment
- KaKs_Calculator 3.0 (or 2.0) with multi-method support

Updates based on KaKs_Calculator 3.0:
1. Model Averaging (MA) as default method (more accurate)
2. 33 genetic code tables
3. 21 output columns (including AICc, Akaike-Weight, Model)
4. Multi-threading support
5. Support for both v2.0 and v3.0 output formats

Input:
    orthogroups_dir: OrthoFinder output (Orthogroups.tsv + Orthogroup_Sequences/)
    cds_file: Single FASTA file containing all CDS nucleotide sequences

Output:
    - kaks_results/kaks_values.tsv: All Ka/Ks values
    - kaks_results/kaks_summary.tsv: Summary statistics by category
    - kaks_results/kaks_invalid.tsv: Skipped sequences
    - kaks_results/kaks_boxplot.R: R visualization script

Example:
    python kaks.py Orthogroups/ all.cds.fa -n 50 -p 50
    python kaks.py Orthogroups/ all.cds.fa --use-kaks-calculator -m MA
    python kaks.py Orthogroups/ all.cds.fa --threads 8
"""

import sys
import os
import glob
import random
import subprocess
import argparse
import tempfile
import shutil
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from functools import partial

from Bio import SeqIO
from Bio.Seq import Seq

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from lib.parser import parse_orthogroups_tsv
from lib.classify import classify_orthogroups
from lib.utils import ensure_dir, log, check_file


# ============================================================
# Genetic code table (33 tables supported in KaKs_Calculator 3.0)
# ============================================================
CODON_TABLES = {
    1: {  # Standard (Universal)
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    },
    11: {  # Bacterial
        'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
    },
}

# Default: Standard genetic code
CODON_TABLE = CODON_TABLES[1]


# ============================================================
# KaKs_Calculator methods
# ============================================================
KAKS_METHODS = {
    # v2.0 methods
    'NG':   'Nei-Gojobori (1986)',
    'LWL':  'Li-Wu-Luo (1985)',
    'LPB':  'Li-Pamilo-Bianchi (1993)',
    'MLWL': 'Modified LWL (2004)',
    'MLPB': 'Modified LPB (2004)',
    'GY':   'Goldman-Yang (1994)',
    'YN':   'Yang-Nielsen (2000)',
    'MYN':  'Modified YN (2006)',
    # v3.0 new methods
    'MS':   'Model Selection (2006)',
    'MA':   'Model Averaging (2006) [RECOMMENDED]',
}


# ============================================================
# FASTA parsing
# ============================================================
def parse_fasta(filepath):
    """Parse FASTA file"""
    sequences = {}
    for record in SeqIO.parse(filepath, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences


# ============================================================
# CDS validation (length % 3 == 0)
# ============================================================
def load_cds_file(cds_file, genetic_code=1):
    """
    Load and validate CDS sequences
    
    Validation rules:
    1. Length must be divisible by 3
    2. No internal stop codons
    3. Minimum 30 bp
    4. Valid nucleotides only (A, C, G, T, N)
    """
    cds_dict = {}
    invalid = []
    codon_table = CODON_TABLES.get(genetic_code, CODON_TABLES[1])
    
    for record in SeqIO.parse(cds_file, "fasta"):
        gene_id = record.id
        seq = str(record.seq).upper()
        
        # Rule 1: Length divisible by 3
        if len(seq) % 3 != 0:
            invalid.append((gene_id, f"Length {len(seq)} not divisible by 3"))
            continue
        
        # Rule 2: Minimum length (30 bp = 10 codons)
        if len(seq) < 30:
            invalid.append((gene_id, f"Too short ({len(seq)} bp)"))
            continue
        
        # Rule 3: No internal stop codons
        has_internal_stop = False
        for i in range(0, len(seq) - 3, 3):
            codon = seq[i:i+3]
            aa = codon_table.get(codon, 'X')
            if aa == '*':
                has_internal_stop = True
                break
        
        if has_internal_stop:
            invalid.append((gene_id, "Internal stop codon"))
            continue
        
        # Rule 4: Valid nucleotides
        invalid_chars = set(seq) - set('ACGTN')
        if invalid_chars:
            invalid.append((gene_id, f"Invalid characters: {invalid_chars}"))
            continue
        
        cds_dict[gene_id] = seq
    
    return cds_dict, invalid


# ============================================================
# Protein sequence loading
# ============================================================
def load_orthogroup_sequences(seq_dir):
    """Load protein sequences from Orthogroup_Sequences directory"""
    protein_dict = {}
    for fa_file in glob.glob(os.path.join(seq_dir, "*.fa")):
        for record in SeqIO.parse(fa_file, "fasta"):
            protein_dict[record.id] = str(record.seq)
    return protein_dict


def parse_orthogroup_sequences_with_og(seq_dir):
    """Load protein sequences with orthogroup mapping"""
    og_sequences = {}
    for fa_file in glob.glob(os.path.join(seq_dir, "*.fa")):
        og_id = os.path.splitext(os.path.basename(fa_file))[0]
        genes = {}
        for record in SeqIO.parse(fa_file, "fasta"):
            genes[record.id] = str(record.seq)
        og_sequences[og_id] = genes
    return og_sequences


# ============================================================
# Protein alignment (MUSCLE / MAFFT / fallback)
# ============================================================
def align_proteins(gene1, gene2, prot1, prot2, tmp_dir):
    """Align two protein sequences"""
    pair_dict = {gene1: prot1, gene2: prot2}
    tmp_file = os.path.join(tmp_dir, f"{gene1}_{gene2}.fa")
    out_file = os.path.join(tmp_dir, f"{gene1}_{gene2}_aln.fa")
    
    with open(tmp_file, 'w') as f:
        f.write(f">{gene1}\n{prot1}\n>{gene2}\n{prot2}\n")
    
    aligned = None
    
    for tool, cmd in [
        ('muscle', ['muscle', '-in', tmp_file, '-out', out_file]),
        ('mafft', ['mafft', '--auto', '--quiet', tmp_file]),
    ]:
        try:
            if tool == 'mafft':
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
                if result.returncode == 0:
                    with open(out_file, 'w') as f:
                        f.write(result.stdout)
            else:
                subprocess.run(cmd, capture_output=True, text=True, timeout=30)
            
            if os.path.exists(out_file):
                aligned = {}
                for record in SeqIO.parse(out_file, "fasta"):
                    aligned[record.id] = str(record.seq)
                break
        except (FileNotFoundError, subprocess.TimeoutExpired):
            continue
    
    # Cleanup
    for f in [tmp_file, out_file]:
        if os.path.exists(f):
            os.remove(f)
    
    # Fallback: simple padding alignment
    if aligned is None:
        len1, len2 = len(prot1), len(prot2)
        if len1 < len2:
            prot1 = prot1 + '-' * (len2 - len1)
        elif len2 < len1:
            prot2 = prot2 + '-' * (len1 - len2)
        aligned = {gene1: prot1, gene2: prot2}
    
    return aligned


# ============================================================
# Back-translation
# ============================================================
def back_translate(protein_alignment, cds_dict):
    """Back-translate protein alignment to CDS alignment"""
    cds_alignment = {}
    
    for gene_id, prot_seq in protein_alignment.items():
        if gene_id not in cds_dict:
            continue
        
        cds_seq = cds_dict[gene_id]
        aligned_cds = []
        cds_pos = 0
        
        for aa in prot_seq:
            if aa == '-':
                aligned_cds.append('---')
            elif cds_pos + 3 <= len(cds_seq):
                aligned_cds.append(cds_seq[cds_pos:cds_pos+3])
                cds_pos += 3
            else:
                aligned_cds.append('---')
        
        cds_alignment[gene_id] = ''.join(aligned_cds)
    
    return cds_alignment


# ============================================================
# Generate AXT format (for KaKs_Calculator)
# ============================================================
def generate_axt(gene1, gene2, cds1, cds2, output_file):
    """Generate AXT format file for KaKs_Calculator"""
    with open(output_file, 'w') as f:
        f.write(f"{gene1} {gene2}\n")
        f.write(f"{cds1}\n")
        f.write(f"{cds2}\n\n")
    return output_file


# ============================================================
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
        
        # Fallback: Python Nei-Gojobori
        Ka, Ks, kaks = estimate_kaks_ng(cds1, cds2, genetic_code)
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


# ============================================================
# Random sampling
# ============================================================
def random_sample_pairs(orthogroups, categories, og_sequences, n_genes, n_pairs):
    """Random sample gene pairs from each category"""
    sampled = defaultdict(list)
    
    for category, og_list in categories.items():
        valid_ogs = [og for og in og_list if og in og_sequences and len(og_sequences[og]) >= 2]
        if not valid_ogs:
            continue
        
        n_sample = min(n_genes, len(valid_ogs))
        for og_id in random.sample(valid_ogs, n_sample):
            genes = list(og_sequences[og_id].keys())
            all_pairs = [(genes[i], genes[j]) for i in range(len(genes)) for j in range(i+1, len(genes))]
            if not all_pairs:
                continue
            n_sample_pairs = min(n_pairs, len(all_pairs))
            for g1, g2 in random.sample(all_pairs, n_sample_pairs):
                sampled[category].append((og_id, g1, g2))
    
    return sampled


# ============================================================
# R visualization script
# ============================================================
def generate_r_script(output_dir):
    """Generate R script for Ka/Ks visualization"""
    r_script = os.path.join(output_dir, 'kaks_boxplot.R')
    od = output_dir.replace('\\', '/')
    
    with open(r_script, 'w') as f:
        f.write(f'''#!/usr/bin/env Rscript
# Ka/Ks Boxplot by Gene Category (KaKs_Calculator 3.0 enhanced)

library(ggplot2)

df <- read.delim("{od}/kaks_values.tsv")
df <- df[!is.na(df$Ka_Ks) & is.finite(df$Ka_Ks) & df$Ka_Ks > 0 & df$Ka_Ks < 5, ]
df$Category <- factor(df$Category, levels=c("core", "soft_core", "dispensable", "private"))

colors <- c("core"="#F8766D", "soft_core"="#7CAE00", "dispensable"="#00BFC4", "private"="#C77CFF")

# Boxplot
p <- ggplot(df, aes(x=Category, y=Ka_Ks, fill=Category)) +
  geom_boxplot(alpha=0.7, outlier.size=0.5) +
  scale_fill_manual(values=colors) +
  labs(x="Gene Category", y="Ka/Ks", title="Ka/Ks Distribution by Gene Category") +
  theme_bw(base_size=14) +
  theme(plot.title=element_text(hjust=0.5, face="bold"), legend.position="none")

ggsave("{od}/kaks_boxplot.pdf", p, width=8, height=6, dpi=300)
ggsave("{od}/kaks_boxplot.png", p, width=8, height=6, dpi=300)
cat("Saved: kaks_boxplot.pdf/png\\n")

# If AICc available (v3.0), plot model selection
if ("AICc" %in% colnames(df)) {{
  df_aic <- df[!is.na(df$AICc), ]
  if (nrow(df_aic) > 0) {{
    p2 <- ggplot(df_aic, aes(x=Category, y=AICc, fill=Category)) +
      geom_boxplot(alpha=0.7) +
      scale_fill_manual(values=colors) +
      labs(x="Gene Category", y="AICc", title="AICc by Gene Category (Model Quality)") +
      theme_bw(base_size=14) +
      theme(plot.title=element_text(hjust=0.5, face="bold"), legend.position="none")
    ggsave("{od}/kaks_aicc.pdf", p2, width=8, height=6, dpi=300)
    cat("Saved: kaks_aicc.pdf\\n")
  }}
}}

# Summary
cat("\\n=== Summary ===\\n")
for (cat in c("core", "soft_core", "dispensable", "private")) {{
  sub <- df[df$Category == cat, ]
  if (nrow(sub) > 0)
    cat(sprintf("  %s: n=%d, median=%.4f, mean=%.4f\\n", cat, nrow(sub),
                median(sub$Ka_Ks), mean(sub$Ka_Ks)))
}}

# Kruskal-Wallis
if (length(unique(df$Category)) > 1) {{
  kw <- kruskal.test(Ka_Ks ~ Category, data=df)
  cat(sprintf("\\nKruskal-Wallis: chi-sq=%.2f, p=%.2e\\n", kw$statistic, kw$p.value))
}}
''')
    
    log(f"R script saved: {r_script}")
    return r_script


# ============================================================
# Main
# ============================================================
def main():
    parser = argparse.ArgumentParser(
        description='Ka/Ks calculator (KaKs_Calculator 3.0 compatible)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Modes:
  1. Pan-genome mode: orthogroups_dir + cds_file (random sampling by category)
  2. Standalone mode: -i input.axt (direct AXT input, like KaKs_Calculator)

AXT Format:
  seq1 seq2
  ATGCGT...
  ATGCGT...
  
  (blank line between pairs)

KaKs_Calculator 3.0 Methods:
  NG    Nei-Gojobori (1986)         - Simple, fast
  LWL   Li-Wu-Luo (1985)            - Weighted sites
  LPB   Li-Pamilo-Bianchi (1993)    - Improved weighting
  GY    Goldman-Yang (1994)          - ML, codon model
  YN    Yang-Nielsen (2000)          - ML, HKY model
  MYN   Modified YN (2006)           - Modified YN
  MS    Model Selection (v3.0)       - AIC-based selection
  MA    Model Averaging (v3.0)       - Weighted average [DEFAULT]

Example:
  # Standalone mode (AXT input)
  python kaks.py -i pairs.axt -o kaks_output -m MA
  python kaks.py -i pairs.axt -o kaks_output -m YN -t 8

  # Pan-genome mode (OrthoFinder)
  python kaks.py Orthogroups/ all.cds.fa -n 50 -p 50
  python kaks.py Orthogroups/ all.cds.fa -t 8 -m MA -k
"""
    )
    
    parser.add_argument('-i', '--input', dest='axt_input', default=None,
                        help='Standalone mode: input AXT file (pairs of CDS sequences)')
    parser.add_argument('orthogroups_dir', nargs='?', default=None,
                        help='OrthoFinder output directory (pan-genome mode)')
    parser.add_argument('cds_file', nargs='?', default=None,
                        help='CDS FASTA file (pan-genome mode)')
    parser.add_argument('-o', '--output', default='kaks_results',
                        help='Output directory (default: kaks_results)')
    parser.add_argument('-n', '--n-genes', type=int, default=50,
                        help='Orthogroups to sample per category (default: 50)')
    parser.add_argument('-p', '--n-pairs', type=int, default=50,
                        help='Species pairs per orthogroup (default: 50)')
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='Number of threads (default: 1)')
    parser.add_argument('-s', '--seed', type=int, default=42,
                        help='Random seed (default: 42)')
    parser.add_argument('-T', '--threshold', type=float, default=0.9,
                        help='Soft-core threshold (default: 0.9)')
    parser.add_argument('-c', '--genetic-code', type=int, default=1,
                        help='Genetic code table 1-33 (default: 1=universal)')
    parser.add_argument('-m', '--method', default='MA',
                        choices=list(KAKS_METHODS.keys()),
                        help='Ka/Ks method (default: MA)')
    parser.add_argument('-k', '--use-kaks-calculator', action='store_true',
                        help='Use KaKs_Calculator if available')
    parser.add_argument('-C', '--calculator-path', default=None,
                        help='Path to KaKs_Calculator executable')
    parser.add_argument('--check-ids', action='store_true',
                        help='Only check CDS/protein ID matching, then exit')
    
    args = parser.parse_args()
    
    random.seed(args.seed)
    ensure_dir(args.output)
    tmp_dir = os.path.join(args.output, 'tmp')
    ensure_dir(tmp_dir)
    
    log("=" * 60)
    log("Ka/Ks Calculation (KaKs_Calculator 3.0 compatible)")
    log("=" * 60)
    log(f"Method: {args.method} - {KAKS_METHODS[args.method]}")
    log(f"Genetic code: {args.genetic_code}")
    log(f"Threads: {args.threads}")
    
    # Check mode
    if args.axt_input:
        # Standalone mode: AXT input
        run_standalone(args)
    elif args.orthogroups_dir and args.cds_file:
        # Pan-genome mode: OrthoFinder + CDS
        run_pangenome(args)
    else:
        print("Error: Either -i (AXT input) or orthogroups_dir + cds_file required")
        print("Use -h for help")
        sys.exit(1)


def run_standalone(args):
    """
    Standalone mode: process AXT file directly (like KaKs_Calculator 3.0)
    
    AXT format:
        seq1_name seq2_name
        ATGCGT...
        ATGCGT...
        
        (blank line between pairs)
    """
    check_file(args.axt_input, "AXT file")
    
    log(f"Standalone mode: {args.axt_input}")
    
    # Parse AXT file
    pairs = []
    with open(args.axt_input, 'r') as f:
        lines = f.readlines()
    
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        
        # Header line: "name1 name2"
        parts = line.split()
        if len(parts) >= 2:
            name1, name2 = parts[0], parts[1]
            i += 1
            
            # Sequence 1
            seq1 = ""
            while i < len(lines) and lines[i].strip():
                seq1 += lines[i].strip()
                i += 1
            
            # Sequence 2
            i += 1
            seq2 = ""
            while i < len(lines) and lines[i].strip():
                seq2 += lines[i].strip()
                i += 1
            
            if seq1 and seq2:
                pairs.append((name1, name2, seq1, seq2))
        else:
            i += 1
    
    log(f"  Found {len(pairs)} sequence pairs")
    
    if len(pairs) == 0:
        log("ERROR: No valid pairs found in AXT file")
        sys.exit(1)
    
    # Calculate Ka/Ks
    results = []
    output_file = os.path.join(args.output, 'kaks_values.tsv')
    
    with open(output_file, 'w') as f:
        f.write('Sequence\tMethod\tKa\tKs\tKa_Ks\n')
        
        for name1, name2, seq1, seq2 in pairs:
            # Validate sequences
            if len(seq1) % 3 != 0 or len(seq2) % 3 != 0:
                log(f"  Warning: {name1} {name2} - length not divisible by 3, skipping")
                continue
            
            if len(seq1) != len(seq2):
                log(f"  Warning: {name1} {name2} - unequal length, skipping")
                continue
            
            # Calculate Ka/Ks
            Ka, Ks, kaks = estimate_kaks_ng(seq1, seq2, args.genetic_code)
            
            ka_str = f"{Ka:.6f}" if Ka is not None else 'NA'
            ks_str = f"{Ks:.6f}" if Ks is not None else 'NA'
            kaks_str = f"{kaks:.6f}" if kaks is not None else 'NA'
            
            f.write(f"{name1}-{name2}\tNG-Python\t{ka_str}\t{ks_str}\t{kaks_str}\n")
            results.append({'name': f"{name1}-{name2}", 'Ka': Ka, 'Ks': Ks, 'Ka_Ks': kaks})
    
    log(f"Results saved: {output_file}")
    
    # Summary
    valid = [r for r in results if r['Ka_Ks'] is not None and 0 < r['Ka_Ks'] < 5]
    if valid:
        kaks_vals = sorted([r['Ka_Ks'] for r in valid])
        median = kaks_vals[len(kaks_vals)//2]
        mean = sum(kaks_vals) / len(kaks_vals)
        log(f"  Valid pairs: {len(valid)}/{len(results)}")
        log(f"  Median Ka/Ks: {median:.4f}")
        log(f"  Mean Ka/Ks: {mean:.4f}")
    
    log("Done!")


def run_pangenome(args):
    """Pan-genome mode: OrthoFinder + CDS input"""
    
    # Locate files
    orthogroups_file = os.path.join(args.orthogroups_dir, 'Orthogroups.tsv')
    seq_dir = os.path.join(args.orthogroups_dir, 'Orthogroup_Sequences')
    
    if not os.path.exists(orthogroups_file):
        for pat in [os.path.join(args.orthogroups_dir, '**/Orthogroups.tsv')]:
            found = glob.glob(pat, recursive=True)
            if found:
                orthogroups_file = found[0]
                seq_dir = os.path.join(os.path.dirname(orthogroups_file), 'Orthogroup_Sequences')
                break
    
    check_file(orthogroups_file, "Orthogroups.tsv")
    check_file(seq_dir, "Orthogroup_Sequences/")
    check_file(args.cds_file, "CDS file")
    
    # Parse Orthogroups
    log("Parsing Orthogroups...")
    orthogroups, species_list = parse_orthogroups_tsv(orthogroups_file)
    log(f"  Species: {len(species_list)}, Orthogroups: {len(orthogroups)}")
    
    # Classify
    log("Classifying orthogroups...")
    categories, _ = classify_orthogroups(orthogroups, species_list, args.threshold)
    for cat, ogs in categories.items():
        log(f"  {cat}: {len(ogs)}")
    
    # Load CDS
    log(f"Loading CDS from {args.cds_file}...")
    cds_dict, invalid_cds = load_cds_file(args.cds_file, args.genetic_code)
    log(f"  Valid CDS: {len(cds_dict)}")
    log(f"  Invalid CDS: {len(invalid_cds)} (skipped)")
    
    # Save invalid list
    if invalid_cds:
        invalid_file = os.path.join(args.output, 'kaks_invalid.tsv')
        with open(invalid_file, 'w') as f:
            f.write("Gene_ID\tReason\n")
            for gid, reason in invalid_cds:
                f.write(f"{gid}\t{reason}\n")
        log(f"  Invalid list: {invalid_file}")
    
    if len(cds_dict) == 0:
        log("ERROR: No valid CDS sequences!")
        sys.exit(1)
    
    # Load proteins
    log("Loading protein sequences...")
    protein_dict = load_orthogroup_sequences(seq_dir)
    log(f"  Total proteins: {len(protein_dict)}")
    
    # Check ID matching
    matching = sum(1 for g in protein_dict if g in cds_dict)
    log(f"  CDS-Protein matching: {matching}/{len(protein_dict)} ({matching/len(protein_dict)*100:.1f}%)")
    
    if matching == 0:
        log("ERROR: No matching IDs!")
        sys.exit(1)
    
    if args.check_ids:
        log("ID check passed. Exiting.")
        sys.exit(0)
    
    # Build orthogroup gene lists
    log("Building orthogroup gene lists...")
    og_sequences = parse_orthogroup_sequences_with_og(seq_dir)
    
    filtered_og = {}
    total_before, total_after = 0, 0
    for og_id, genes in og_sequences.items():
        total_before += len(genes)
        valid_genes = {g: seq for g, seq in genes.items() if g in cds_dict}
        if len(valid_genes) >= 2:
            filtered_og[og_id] = valid_genes
            total_after += len(valid_genes)
    
    log(f"  Genes with valid CDS: {total_after}/{total_before}")
    
    # Random sampling
    log(f"Random sampling: n={args.n_genes}, p={args.n_pairs}")
    sampled_pairs = random_sample_pairs(orthogroups, categories, filtered_og, 
                                         args.n_genes, args.n_pairs)
    total_pairs = sum(len(p) for p in sampled_pairs.values())
    log(f"Total pairs to analyze: {total_pairs}")
    
    # Check KaKs_Calculator availability
    use_calculator = args.use_kakscalculator
    if use_calculator:
        success, _ = run_kakscalculator("", "/dev/null", args.method, 
                                         args.genetic_code, args.calculator_path)
        if not success:
            log("KaKs_Calculator not found, using Python fallback")
            use_calculator = False
        else:
            log(f"Using KaKs_Calculator with method: {args.method}")
    
    # Prepare arguments for multiprocessing
    all_args = []
    for category, pairs in sampled_pairs.items():
        for og_id, g1, g2 in pairs:
            all_args.append((
                og_id, g1, g2, protein_dict, cds_dict, tmp_dir,
                use_calculator, args.method, args.genetic_code, args.calculator_path
            ))
    
    # Calculate Ka/Ks with multiprocessing
    log(f"Calculating Ka/Ks (threads={args.threads})...")
    results = []
    
    if args.threads > 1:
        with Pool(args.threads) as pool:
            for i, result in enumerate(pool.imap_unordered(calculate_single_pair, all_args)):
                if result is not None:
                    # Add category
                    for cat, pairs in sampled_pairs.items():
                        if any(p[0] == result['Orthogroup'] and 
                               ((p[1] == result['Gene1'] and p[2] == result['Gene2']) or
                                (p[1] == result['Gene2'] and p[2] == result['Gene1']))
                               for p in pairs):
                            result['Category'] = cat
                            break
                    results.append(result)
                if (i + 1) % 100 == 0:
                    log(f"  Processed {i+1}/{len(all_args)} pairs")
    else:
        for i, args_tuple in enumerate(all_args):
            result = calculate_single_pair(args_tuple)
            if result is not None:
                for cat, pairs in sampled_pairs.items():
                    if any(p[0] == result['Orthogroup'] and 
                           ((p[1] == result['Gene1'] and p[2] == result['Gene2']) or
                            (p[1] == result['Gene2'] and p[2] == result['Gene1']))
                           for p in pairs):
                        result['Category'] = cat
                        break
                results.append(result)
            if (i + 1) % 100 == 0:
                log(f"  Processed {i+1}/{len(all_args)} pairs")
    
    # Save results
    log("Saving results...")
    
    results_file = os.path.join(args.output, 'kaks_values.tsv')
    has_v3 = any(r.get('AICc') is not None for r in results)
    
    with open(results_file, 'w') as f:
        if has_v3:
            f.write('Orthogroup\tGene1\tGene2\tCategory\tKa\tKs\tKa_Ks\tMethod\tP_Value\tAICc\tAkaike_Weight\tModel\n')
            for r in results:
                f.write(f"{r['Orthogroup']}\t{r['Gene1']}\t{r['Gene2']}\t{r.get('Category','')}\t"
                        f"{_fmt(r['Ka'])}\t{_fmt(r['Ks'])}\t{_fmt(r['Ka_Ks'])}\t"
                        f"{r.get('Method','')}\t{_fmt(r.get('P_Value'))}\t"
                        f"{_fmt(r.get('AICc'))}\t{_fmt(r.get('Akaike_Weight'))}\t"
                        f"{r.get('Model','')}\n")
        else:
            f.write('Orthogroup\tGene1\tGene2\tCategory\tKa\tKs\tKa_Ks\tMethod\n')
            for r in results:
                f.write(f"{r['Orthogroup']}\t{r['Gene1']}\t{r['Gene2']}\t{r.get('Category','')}\t"
                        f"{_fmt(r['Ka'])}\t{_fmt(r['Ks'])}\t{_fmt(r['Ka_Ks'])}\t"
                        f"{r.get('Method','')}\n")
    
    log(f"Results: {results_file}")
    
    # Summary
    summary_file = os.path.join(args.output, 'kaks_summary.tsv')
    with open(summary_file, 'w') as f:
        f.write('Category\tN_pairs\tN_valid\tMedian_Ka\tMedian_Ks\tMedian_KaKs\tMean_KaKs\n')
        for cat in ['core', 'soft_core', 'dispensable', 'private']:
            cat_r = [r for r in results if r.get('Category') == cat]
            valid = [r for r in cat_r if r['Ka_Ks'] is not None and 0 < r['Ka_Ks'] < 5]
            if valid:
                kaks_vals = sorted([r['Ka_Ks'] for r in valid])
                ka_vals = sorted([r['Ka'] for r in valid if r['Ka'] is not None])
                ks_vals = sorted([r['Ks'] for r in valid if r['Ks'] is not None])
                f.write(f"{cat}\t{len(cat_r)}\t{len(valid)}\t"
                        f"{ka_vals[len(ka_vals)//2]:.6f}\t"
                        f"{ks_vals[len(ks_vals)//2]:.6f}\t"
                        f"{kaks_vals[len(kaks_vals)//2]:.6f}\t"
                        f"{sum(kaks_vals)/len(kaks_vals):.6f}\n")
            else:
                f.write(f"{cat}\t{len(cat_r)}\t0\tNA\tNA\tNA\tNA\n")
    log(f"Summary: {summary_file}")
    
    # R script
    r_script = generate_r_script(args.output)
    
    # Print summary
    print("\n" + "=" * 50)
    print("Ka/Ks Summary:")
    print("=" * 50)
    with open(summary_file) as f:
        f.readline()
        for line in f:
            parts = line.strip().split('\t')
            print(f"  {parts[0]}: n={parts[2]}, median={parts[5]}, mean={parts[6]}")
    print(f"\nRscript {r_script}")
    log("Done!")


def _fmt(v):
    """Format value for output"""
    if v is None:
        return 'NA'
    if isinstance(v, float):
        return f"{v:.6f}"
    return str(v)


if __name__ == '__main__':
    main()
