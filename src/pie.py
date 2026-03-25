"""
Pie chart - gene family category proportions
"""
import sys
import os
import subprocess

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.parser import parse_frequency_table
from lib.utils import ensure_dir, log, check_file


def register(subparsers):
    p = subparsers.add_parser('pie', help='Pie chart: gene family category proportions')
    p.add_argument('frequency', help='Frequency table file (frequency_table.tsv)')
    p.add_argument('-o', '--output', default='pgkit_output', help='Output directory')
    p.set_defaults(func=run)


def run(args):
    ensure_dir(args.output)

    log("Reading frequency table...")
    check_file(args.frequency, "Frequency table")
    data = parse_frequency_table(args.frequency)

    # Count categories
    cat_counts = {}
    for row in data:
        cat = row['Category']
        cat_counts[cat] = cat_counts.get(cat, 0) + 1

    total = sum(cat_counts.values())
    log(f"Total: {total} orthogroups")
    for cat, count in cat_counts.items():
        log(f"  {cat}: {count} ({count/total*100:.1f}%)")

    # Call R script
    r_script = os.path.join(os.path.dirname(__file__), 'scripts', 'plot_pie.R')
    out_prefix = os.path.join(args.output, 'pie')

    log("Generating pie chart...")
    try:
        result = subprocess.run(
            ['Rscript', os.path.abspath(r_script), args.frequency, out_prefix],
            capture_output=True, text=True, timeout=120
        )
        if result.returncode == 0:
            print(result.stdout)
            log("Done!")
        else:
            print(result.stderr)
            log("R script failed, check R environment")
    except FileNotFoundError:
        log("Rscript not found, please install R")
        log(f"Manual run: Rscript {os.path.abspath(r_script)} {args.frequency} {out_prefix}")
