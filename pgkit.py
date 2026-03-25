#!/usr/bin/env python3
"""
pgkit - Pan-genome analysis toolkit

Pan-genome analysis based on OrthoFinder output:
- PAV matrix construction and gene family classification
- Core/Pan gene family saturation curve
- Visualization (pie, bar, heatmap)
- Statistics report
- Ka/Ks calculation

Usage:
    pgkit <command> [options]

Commands:
    pav     Build PAV matrix + classification + auto visualization
    curve   Core/Pan gene family saturation curve
    pie     Pie chart: gene family category proportions
    bar     Bar chart: gene family composition by sample
    heatmap Heatmap: PAV matrix visualization
    stats   Generate statistics report
    kaks    Ka/Ks calculation by gene category

Example:
    pgkit pav Orthogroups/ -o results
    pgkit curve results/pav_matrix.tsv -o results -s 100
    pgkit kaks Orthogroups/ all.cds.fa -t 8 -m MA
"""

import argparse
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from src import pav, curve, pie, bar, heatmap, stats, kaks


def main():
    parser = argparse.ArgumentParser(
        prog='pgkit',
        description='Pan-gene Family Analysis Toolkit',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Workflow:
  1. pav     Build PAV matrix + classification + visualization
  2. curve   Saturation curve analysis (optional)
  3. pie     Pie chart visualization
  4. bar     Bar chart visualization
  5. heatmap Heatmap visualization
  6. stats   Statistics report
  7. kaks    Ka/Ks calculation

Examples:
  pgkit pav Orthogroups/ -o results
  pgkit pav Orthogroups/ -o results -f pdf -r
  pgkit curve results/pav_matrix.tsv -o results -s 100
  pgkit stats results/frequency_table.tsv -g results/gene_count_matrix.tsv -o results
  pgkit kaks Orthogroups/ all.cds.fa -t 8 -m MA -k
"""
    )

    subparsers = parser.add_subparsers(dest='command', help='Subcommand')

    pav.register(subparsers)
    curve.register(subparsers)
    pie.register(subparsers)
    bar.register(subparsers)
    heatmap.register(subparsers)
    stats.register(subparsers)
    kaks.register(subparsers)

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    args.func(args)


if __name__ == '__main__':
    main()
