#!/usr/bin/env python3
"""Ka/Ks command-line workflow."""

import argparse
import glob
import os
import random
import shutil
import subprocess
import sys
from multiprocessing import Pool

from pgkit.core.classify import classify_orthogroups
from pgkit.core.parser import parse_orthogroups_tsv
from pgkit.core.utils import check_file, ensure_dir, log
from pgkit.kaks.calculator import _fmt, calculate_single_pair, estimate_kaks_ng
from pgkit.kaks.constants import KAKS_METHODS
from pgkit.kaks.io import load_cds_file, load_orthogroup_sequences, parse_orthogroup_sequences_with_og
from pgkit.kaks.report import generate_r_script
from pgkit.kaks.sampling import random_sample_pairs

# ============================================================
# Subcommand registration (for pgkit.py)
# ============================================================
def register(subparsers):
    """Register kaks subcommand"""
    p = subparsers.add_parser('kaks', help='Ka/Ks calculation (KaKs_Calculator 3.0 compatible)')
    p.add_argument('-i', '--input', dest='axt_input', default=None,
                   help='Standalone mode: input AXT file')
    p.add_argument('orthogroups_dir', nargs='?', default=None,
                   help='OrthoFinder output directory (pan-genome mode)')
    p.add_argument('cds_file', nargs='?', default=None,
                   help='CDS FASTA file (pan-genome mode)')
    p.add_argument('-o', '--output', default='kaks_results',
                   help='Output directory (default: kaks_results)')
    p.add_argument('-n', '--n-genes', type=int, default=50,
                   help='Orthogroups to sample per category (default: 50)')
    p.add_argument('-p', '--n-pairs', type=int, default=50,
                   help='Species pairs per orthogroup (default: 50)')
    p.add_argument('-t', '--threads', type=int, default=1,
                   help='Number of threads (default: 1)')
    p.add_argument('-s', '--seed', type=int, default=42,
                   help='Random seed (default: 42)')
    p.add_argument('-T', '--threshold', type=float, default=0.9,
                   help='Soft-core threshold (default: 0.9)')
    p.add_argument('-c', '--genetic-code', type=int, default=1,
                   help='Genetic code table 1-33 (default: 1=universal)')
    p.add_argument('-m', '--method', default='MA',
                   choices=list(KAKS_METHODS.keys()),
                   help='Ka/Ks method (default: MA)')
    p.add_argument('-k', '--use-kaks-calculator', action='store_true',
                   help='Use KaKs_Calculator if available')
    p.add_argument('-C', '--calculator-path', default=None,
                   help='Path to KaKs_Calculator executable')
    p.add_argument('--check-ids', action='store_true',
                   help='Only check CDS/protein ID matching, then exit')
    p.set_defaults(func=run_cli)


def run_cli(args):
    """Entry point for subcommand"""
    # Convert argparse namespace to work with main logic
    if args.axt_input:
        run_standalone(args)
    elif args.orthogroups_dir and args.cds_file:
        run_pangenome(args)
    else:
        print("Error: Either -i (AXT input) or orthogroups_dir + cds_file required")
        sys.exit(1)


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

    # Setup directories
    ensure_dir(args.output)
    tmp_dir = os.path.join(args.output, 'tmp')
    ensure_dir(tmp_dir)

    # Locate Orthogroups.tsv
    orthogroups_file = os.path.join(args.orthogroups_dir, 'Orthogroups.tsv')

    if not os.path.exists(orthogroups_file):
        for pat in [os.path.join(args.orthogroups_dir, '**/Orthogroups.tsv')]:
            found = glob.glob(pat, recursive=True)
            if found:
                orthogroups_file = found[0]
                break

    check_file(orthogroups_file, "Orthogroups.tsv")

    # Locate Orthogroup_Sequences/ (same level as Orthogroups/)
    # OrthoFinder structure:
    #   Orthogroups/
    #   Orthogroup_Sequences/
    orthogroups_dir = os.path.dirname(orthogroups_file)
    parent_dir = os.path.dirname(orthogroups_dir)

    seq_dir = None
    for candidate in [
        os.path.join(parent_dir, 'Orthogroup_Sequences'),  # same level as Orthogroups/
        os.path.join(args.orthogroups_dir, 'Orthogroup_Sequences'),  # input dir itself
    ]:
        if os.path.isdir(candidate):
            seq_dir = candidate
            break

    if seq_dir is None:
        log("ERROR: Orthogroup_Sequences/ not found!")
        log(f"  Searched: {parent_dir}/Orthogroup_Sequences")
        log(f"  Searched: {args.orthogroups_dir}/Orthogroup_Sequences")
        sys.exit(1)

    log(f"  Orthogroup_Sequences: {seq_dir}")
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
    use_calculator = args.use_kaks_calculator
    calc_path = None

    if use_calculator:
        # Find calculator path
        calc_path = args.calculator_path
        if calc_path:
            calc_path = os.path.expanduser(calc_path)

        if calc_path is None:
            for path in ['KaKs_Calculator', 'KaKs', 'KaKs.exe']:
                resolved = shutil.which(path)
                if resolved:
                    calc_path = resolved
                    break
        else:
            resolved = shutil.which(calc_path)
            if resolved:
                calc_path = resolved

        if calc_path:
            log(f"Using KaKs_Calculator: {calc_path} (method: {args.method})")
        else:
            log("KaKs_Calculator not found! Use -C to specify path")
            log("  Example: -C /path/to/KaKs_Calculator-3.0/bin/KaKs")
            log("  Falling back to Python NG method")
            use_calculator = False

    # Prepare arguments for multiprocessing
    all_args = []
    for category, pairs in sampled_pairs.items():
        for og_id, g1, g2 in pairs:
            all_args.append((
                og_id, g1, g2, protein_dict, cds_dict, tmp_dir,
                use_calculator, args.method, args.genetic_code, calc_path
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

    # Count valid results
    n_valid = sum(1 for r in results if r.get('Ka_Ks') is not None and r['Ka_Ks'] > 0)
    log(f"  Total results: {len(results)}, Valid Ka/Ks: {n_valid}")

    # Debug: show category distribution
    for cat in ['Core', 'Softcore', 'Dispensable', 'Private']:
        cat_count = sum(1 for r in results if r.get('Category') == cat)
        log(f"  {cat}: {cat_count} results")

    # Debug: show method distribution
    methods = {}
    for r in results:
        m = r.get('Method', 'unknown')
        methods[m] = methods.get(m, 0) + 1
    log("  Methods used:")
    for m, count in sorted(methods.items()):
        log(f"    {m}: {count}")

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
        for cat in ['Core', 'Softcore', 'Dispensable', 'Private']:
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

    # Auto run R script
    log("Generating Ka/Ks boxplot...")
    out_prefix = os.path.join(args.output, 'kaks')
    try:
        result = subprocess.run(
            ['Rscript', r_script, results_file, out_prefix],
            capture_output=True, text=True, timeout=120
        )
        if result.returncode == 0:
            print(result.stdout)
        else:
            log("R script failed, run manually:")
            log(f"  Rscript {r_script} {results_file} {out_prefix}")
    except FileNotFoundError:
        log("Rscript not found, run manually:")
        log(f"  Rscript {r_script} {results_file} {out_prefix}")

    # Print summary
    print("\n" + "=" * 50)
    print("Ka/Ks Summary:")
    print("=" * 50)
    with open(summary_file) as f:
        f.readline()
        for line in f:
            parts = line.strip().split('\t')
            print(f"  {parts[0]}: n={parts[2]}, median={parts[5]}, mean={parts[6]}")

    log("Done!")


if __name__ == '__main__':
    main()
