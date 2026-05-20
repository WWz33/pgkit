"""Command-line interface for pgkit."""

from __future__ import annotations

import argparse
import sys
from typing import List, Optional

from pgkit.commands import bar, curve, heatmap, kaks, pav, pie, run, stats


def build_parser() -> argparse.ArgumentParser:
    """Build the top-level pgkit argument parser."""
    parser = argparse.ArgumentParser(
        prog="pgkit",
        description="Pan-gene Family Analysis Toolkit",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Commands:
  run     Full pipeline: PAV + curve + visualization
  pav     Build PAV matrix + classification + visualization
  curve   Saturation curve analysis
  pie     Pie chart
  bar     Bar chart
  heatmap Heatmap visualization (supports --pop for population annotation)
  stats   Statistics report
  kaks    Ka/Ks calculation

Examples:
  pgkit run Orthogroups/ -o results
  pgkit pav Orthogroups/ -o results -f pdf
  pgkit heatmap results/pav_matrix.tsv -f results/frequency_table.tsv -P pop.tsv
  pgkit kaks Orthogroups/ all.cds.fa -t 8 -m MA -k
""",
    )

    subparsers = parser.add_subparsers(dest="command", help="Subcommand")
    run.register(subparsers)
    pav.register(subparsers)
    curve.register(subparsers)
    pie.register(subparsers)
    bar.register(subparsers)
    heatmap.register(subparsers)
    stats.register(subparsers)
    kaks.register(subparsers)
    return parser


def main(argv: Optional[List[str]] = None) -> None:
    """Run the pgkit command-line interface."""
    parser = build_parser()
    args = parser.parse_args(argv)

    if not args.command:
        parser.print_help()
        sys.exit(1)

    args.func(args)
