"""
Heatmap - PAV matrix visualization
"""
import sys
import os
import subprocess

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from lib.utils import ensure_dir, log, check_file


def register(subparsers):
    p = subparsers.add_parser('heatmap', help='Heatmap: PAV matrix visualization')
    p.add_argument('pav', help='PAV matrix file (pav_matrix.tsv)')
    p.add_argument('-f', '--frequency', help='Frequency table for row annotation')
    p.add_argument('-P', '--pop', help='Population file (2-col TSV: species, group)')
    p.add_argument('-o', '--output', default='pgkit_output', help='Output directory')
    p.set_defaults(func=run)


def run(args):
    ensure_dir(args.output)

    log("Reading PAV matrix...")
    check_file(args.pav, "PAV matrix")

    with open(args.pav, 'r') as f:
        header = f.readline().strip().split('\t')
        n_genes = sum(1 for _ in f)
        n_samples = len(header) - 1
    log(f"Matrix dimensions: {n_genes} genes x {n_samples} samples")

    # Build command
    r_script = os.path.join(os.path.dirname(__file__), 'scripts', 'plot_heatmap.R')
    out_prefix = os.path.join(args.output, 'heatmap')

    cmd = ['Rscript', os.path.abspath(r_script), args.pav, out_prefix]

    # Add frequency file
    if args.frequency and os.path.exists(args.frequency):
        cmd.append(args.frequency)
    else:
        cmd.append("NULL")

    # Add pop file
    if args.pop and os.path.exists(args.pop):
        cmd.append(args.pop)
        log(f"Population file: {args.pop}")

    log("Generating heatmap...")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode == 0:
            print(result.stdout)
            log("Done!")
        else:
            print(result.stderr)
            log("R script failed, check R environment")
    except FileNotFoundError:
        log("Rscript not found, please install R")
        log(f"Manual run: Rscript {' '.join(cmd[1:])}")
