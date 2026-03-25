"""
PAV matrix construction + gene family classification + auto visualization

Input: OrthoFinder output directory or Orthogroups.tsv file
Auto-detect and merge Orthogroups_UnassignedGenes.tsv

Merged: PAV matrix + frequency table + gene count matrix + classification + visualization
"""
import sys
import os
import glob
import subprocess
import shutil

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.parser import parse_orthogroups_tsv, parse_unassigned_genes
from lib.classify import classify_orthogroups, build_og_to_category, count_genes_per_species_per_category
from lib.utils import ensure_dir, log, check_file, write_tsv


def register(subparsers):
    """Register subcommand"""
    p = subparsers.add_parser('pav', help='Build PAV matrix + classification + auto visualization')
    p.add_argument('input', help='OrthoFinder output directory or Orthogroups.tsv file')
    p.add_argument('-o', '--output', default='pgkit_output', help='Output directory')
    p.add_argument('-t', '--threshold', type=float, default=0.9,
                   help='Soft-core threshold (default: 0.9)')
    p.add_argument('-f', '--format', choices=['png', 'pdf', 'svg'], default='png',
                   help='Output image format (default: png)')
    p.add_argument('-r', '--save-r', action='store_true',
                   help='Save R scripts to output directory for later modification')
    p.add_argument('-n', '--no-plot', action='store_true',
                   help='Skip visualization step')
    p.add_argument('-s', '--stacked', action='store_true',
                   help='Use stacked bar chart')
    p.set_defaults(func=run)


def detect_files(input_path):
    """Auto-detect OrthoFinder output files"""
    orthogroups_file = None
    unassigned_file = None

    if os.path.isdir(input_path):
        log(f"Scanning directory: {input_path}")
        patterns = [
            os.path.join(input_path, 'Orthogroups.tsv'),
            os.path.join(input_path, '**/Orthogroups.tsv'),
        ]
        for pat in patterns:
            found = glob.glob(pat, recursive=True)
            if found:
                orthogroups_file = found[0]
                break

        if orthogroups_file:
            dir_path = os.path.dirname(orthogroups_file)
            unassigned = os.path.join(dir_path, 'Orthogroups_UnassignedGenes.tsv')
            if os.path.exists(unassigned):
                unassigned_file = unassigned

    elif os.path.isfile(input_path):
        orthogroups_file = input_path
        dir_path = os.path.dirname(input_path)
        unassigned = os.path.join(dir_path, 'Orthogroups_UnassignedGenes.tsv')
        if os.path.exists(unassigned):
            unassigned_file = unassigned
    else:
        print(f"Error: Input path does not exist: {input_path}")
        sys.exit(1)

    return orthogroups_file, unassigned_file


def parse_unassigned_to_orthogroups_format(unassigned_file):
    """Convert UnassignedGenes to Orthogroups format"""
    orthogroups = {}
    species_list = []

    with open(unassigned_file, 'r', encoding='utf-8') as f:
        header = f.readline().strip().split('\t')
        species_list = header[1:]

        idx = 0
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue

            for i, sp in enumerate(species_list):
                genes_str = parts[i + 1].strip() if i + 1 < len(parts) else ''
                if genes_str and genes_str != '-':
                    genes = [g.strip() for g in genes_str.split(',')]
                    for gene in genes:
                        if gene:
                            idx += 1
                            og_id = f'OG_Unassigned_{idx:07d}'
                            og_data = {s: '' for s in species_list}
                            og_data[sp] = gene
                            orthogroups[og_id] = og_data

    return orthogroups, species_list


def run_r_script(script_name, args_list):
    """Run an R script"""
    r_script = os.path.join(os.path.dirname(__file__), '..', 'scripts', script_name)
    cmd = ['Rscript', os.path.abspath(r_script)] + args_list
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if result.returncode == 0:
            print(result.stdout)
            return True
        else:
            print(result.stderr)
            log(f"R script {script_name} failed")
            return False
    except FileNotFoundError:
        log("Rscript not found, please install R")
        log(f"Manual run: {' '.join(cmd)}")
        return False


def copy_r_scripts(output_dir):
    """Copy R scripts to output directory"""
    scripts_dir = os.path.join(os.path.dirname(__file__), '..', 'scripts')
    dst_dir = os.path.join(output_dir, 'r_scripts')
    ensure_dir(dst_dir)

    for fname in ['plot_pie.R', 'plot_bar.R', 'plot_heatmap.R', 'plot_curve.R']:
        src = os.path.join(scripts_dir, fname)
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(dst_dir, fname))
    log(f"R scripts saved to: {dst_dir}")


def run(args):
    """Execute PAV construction and classification"""
    ensure_dir(args.output)

    # Auto-detect files
    orthogroups_file, unassigned_file = detect_files(args.input)

    if not orthogroups_file:
        print("Error: Orthogroups.tsv not found")
        sys.exit(1)

    log(f"Orthogroups file: {orthogroups_file}")
    if unassigned_file:
        log(f"UnassignedGenes file: {unassigned_file}")

    # Parse Orthogroups
    log("Parsing Orthogroups file...")
    orthogroups, species_list = parse_orthogroups_tsv(orthogroups_file)
    log(f"  Clustered gene families: {len(orthogroups)}")

    # Parse UnassignedGenes
    if unassigned_file:
        log("Parsing UnassignedGenes file...")
        unassigned_ogs, _ = parse_unassigned_to_orthogroups_format(unassigned_file)
        log(f"  Unassigned genes: {len(unassigned_ogs)}")
        orthogroups.update(unassigned_ogs)

    total_ogs = len(orthogroups)
    total_species = len(species_list)
    log(f"Species: {total_species}, Total orthogroups: {total_ogs}")

    # Classify
    log("Classifying gene families...")
    categories, species_counts = classify_orthogroups(
        orthogroups, species_list, args.threshold
    )

    # Save classification lists
    for cat, ogs in categories.items():
        filepath = os.path.join(args.output, f'{cat}_orthogroups.txt')
        with open(filepath, 'w') as f:
            for og in sorted(ogs):
                f.write(f"{og}\n")
        log(f"  {cat}: {len(ogs)} ({len(ogs)/total_ogs*100:.1f}%)")

    # Build PAV matrix
    log("Building PAV matrix...")
    pav_rows = []
    for og_id in sorted(orthogroups.keys()):
        og_data = orthogroups[og_id]
        row = [og_id]
        for sp in species_list:
            val = og_data.get(sp, '')
            row.append('1' if val and val.strip() else '0')
        pav_rows.append(row)
    pav_file = os.path.join(args.output, 'pav_matrix.tsv')
    write_tsv(pav_file, ['Orthogroup'] + species_list, pav_rows)

    # Build frequency table
    log("Building frequency table...")
    og_to_cat = build_og_to_category(categories)
    freq_rows = []
    for og_id in sorted(orthogroups.keys()):
        count = species_counts[og_id]
        freq = count / total_species
        cat = og_to_cat.get(og_id, 'unknown')
        freq_rows.append([og_id, f"{freq:.4f}", cat, count])
    freq_file = os.path.join(args.output, 'frequency_table.tsv')
    write_tsv(freq_file, ['Orthogroup', 'Frequency', 'Category', 'Species_Count'], freq_rows)

    # Build gene-category-species table
    log("Building gene-category-species table...")
    gene_cat_rows = []
    for og_id in sorted(orthogroups.keys()):
        og_data = orthogroups[og_id]
        cat = og_to_cat.get(og_id, 'unknown')
        for sp in species_list:
            genes_str = og_data.get(sp, '')
            if genes_str and genes_str.strip():
                genes = [g.strip() for g in genes_str.split(',')]
                for gene in genes:
                    if gene:
                        gene_cat_rows.append([og_id, gene, cat, sp])
    gene_cat_file = os.path.join(args.output, 'gene_category.tsv')
    write_tsv(gene_cat_file, ['Orthogroup', 'Gene', 'Category', 'Species'], gene_cat_rows)
    log(f"  Total genes: {len(gene_cat_rows)}")

    # Build gene count matrix
    log("Building gene count matrix...")
    gene_counts = count_genes_per_species_per_category(orthogroups, categories, species_list)
    count_rows = []
    for sp in species_list:
        core = gene_counts[sp]['core']
        soft = gene_counts[sp]['soft_core']
        disp = gene_counts[sp]['dispensable']
        priv = gene_counts[sp]['private']
        total = core + soft + disp + priv
        count_rows.append([sp, core, soft, disp, priv, total])
    count_file = os.path.join(args.output, 'gene_count_matrix.tsv')
    write_tsv(count_file, ['Species', 'Core', 'Soft-core', 'Dispensable', 'Private', 'Total'], count_rows)

    # Save R scripts
    if args.save_r:
        copy_r_scripts(args.output)

    # Auto visualization
    if not args.no_plot:
        fmt = args.format
        out_prefix = os.path.join(args.output, 'pgkit')
        stacked = "TRUE" if args.stacked else "FALSE"

        log(f"Generating visualizations (format: {fmt})...")

        log("  Generating pie chart...")
        run_r_script('plot_pie.R', [freq_file, out_prefix, fmt])

        log("  Generating bar chart...")
        run_r_script('plot_bar.R', [count_file, out_prefix, stacked, fmt])

        log("  Generating heatmap...")
        run_r_script('plot_heatmap.R', [pav_file, out_prefix, freq_file, fmt])

    # Summary
    log("Output files:")
    log(f"  PAV matrix: {pav_file}")
    log(f"  Frequency table: {freq_file}")
    log(f"  Gene-category-species: {gene_cat_file}")
    log(f"  Gene count matrix: {count_file}")
    log(f"  Classification lists: {args.output}/*_orthogroups.txt")
    if not args.no_plot:
        log(f"  Visualizations: {args.output}/pgkit.*.{fmt}")
    if args.save_r:
        log(f"  R scripts: {args.output}/r_scripts/")
    log("Done!")
