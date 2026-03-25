"""
Statistics report generation
"""
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.parser import parse_frequency_table
from lib.utils import ensure_dir, log, check_file


def register(subparsers):
    p = subparsers.add_parser('stats', help='Generate statistics report')
    p.add_argument('frequency', help='Frequency table file (frequency_table.tsv)')
    p.add_argument('-g', '--gene-count', help='Gene count matrix (gene_count_matrix.tsv)')
    p.add_argument('-o', '--output', default='pgkit_output', help='Output directory')
    p.set_defaults(func=run)


def run(args):
    ensure_dir(args.output)

    log("Generating statistics report...")
    check_file(args.frequency, "Frequency table")

    data = parse_frequency_table(args.frequency)
    total_ogs = len(data)

    # Category statistics
    cat_counts = {}
    for row in data:
        cat = row['Category']
        cat_counts[cat] = cat_counts.get(cat, 0) + 1

    # Species count
    n_species = 0
    for row in data:
        sc = row['Species_Count']
        if isinstance(sc, str):
            sc = int(sc)
        n_species = max(n_species, sc)

    report_file = os.path.join(args.output, 'statistics_report.txt')
    with open(report_file, 'w', encoding='utf-8') as f:
        f.write("=" * 60 + "\n")
        f.write("          Pan-genome Analysis Statistics Report\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Total orthogroups: {total_ogs}\n")
        f.write(f"Max species presence: {n_species}\n\n")

        f.write("-" * 40 + "\n")
        f.write("Category Statistics:\n")
        f.write("-" * 40 + "\n\n")

        cat_order = ['core', 'soft_core', 'dispensable', 'private']
        cat_names = {
            'core': 'Core (present in all samples)',
            'soft_core': 'Soft-core (present in >=90% samples)',
            'dispensable': 'Dispensable (present in some samples)',
            'private': 'Private (present in single sample)'
        }

        for cat in cat_order:
            count = cat_counts.get(cat, 0)
            pct = count / total_ogs * 100 if total_ogs > 0 else 0
            f.write(f"  {cat_names.get(cat, cat)}:\n")
            f.write(f"    Count: {count}\n")
            f.write(f"    Percentage: {pct:.2f}%\n\n")

        # Gene count statistics
        if args.gene_count and os.path.exists(args.gene_count):
            f.write("-" * 40 + "\n")
            f.write("Per-species Gene Count Statistics:\n")
            f.write("-" * 40 + "\n\n")

            with open(args.gene_count, 'r') as cf:
                header = cf.readline().strip().split('\t')
                f.write(f"  {'Species':<20} {'Core':<8} {'Soft':<8} {'Disp':<8} {'Priv':<8} {'Total':<8}\n")
                f.write("  " + "-" * 60 + "\n")
                for line in cf:
                    parts = line.strip().split('\t')
                    f.write(f"  {parts[0]:<20} {parts[1]:<8} {parts[2]:<8} {parts[3]:<8} {parts[4]:<8} {parts[5]:<8}\n")

        f.write("\n" + "=" * 60 + "\n")

    log(f"Report saved: {report_file}")

    # Print summary
    print("\n" + "=" * 40)
    print("Summary:")
    print("=" * 40)
    for cat in cat_order:
        count = cat_counts.get(cat, 0)
        pct = count / total_ogs * 100 if total_ogs > 0 else 0
        print(f"  {cat}: {count} ({pct:.1f}%)")
    print()
    log("Done!")
