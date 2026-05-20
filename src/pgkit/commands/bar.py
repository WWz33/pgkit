"""
Bar chart - gene family composition by sample
"""
import sys
import os
import subprocess

from pgkit.core.utils import check_file, ensure_dir, log, script_path


def register(subparsers):
    p = subparsers.add_parser('bar', help='Bar chart: gene family composition by sample')
    p.add_argument('gene_count', help='Gene count matrix (gene_count_matrix.tsv)')
    p.add_argument('-o', '--output', default='pgkit_output', help='Output directory')
    p.add_argument('-s', '--stacked', action='store_true', help='Use stacked bar chart')
    p.set_defaults(func=run)


def run(args):
    ensure_dir(args.output)

    log("Reading gene count matrix...")
    check_file(args.gene_count, "Gene count matrix")

    with open(args.gene_count, 'r') as f:
        f.readline()
        n_samples = sum(1 for _ in f)
    log(f"Samples: {n_samples}")

    # Call R script
    r_script = script_path('plot_bar.R')
    out_prefix = os.path.join(args.output, 'bar')
    stacked = "TRUE" if args.stacked else "FALSE"

    log("Generating bar chart...")
    try:
        result = subprocess.run(
            ['Rscript', os.path.abspath(r_script), args.gene_count, out_prefix, stacked],
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
        log(f"Manual run: Rscript {os.path.abspath(r_script)} {args.gene_count} {out_prefix} {stacked}")
