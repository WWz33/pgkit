"""
Run full pipeline: PAV + saturation curve
"""
import sys
import os

from pgkit.commands import curve, group, pav
from pgkit.core.utils import log


def register(subparsers):
    p = subparsers.add_parser('run', help='Full pipeline: PAV + classification + curve + visualization')
    p.add_argument('input', help='OrthoFinder output directory or Orthogroups.tsv')
    p.add_argument('-o', '--output', default='pgkit_output', help='Output directory')
    p.add_argument('-t', '--threshold', type=float, default=0.9, help='Soft-core threshold')
    p.add_argument('-f', '--format', choices=['png', 'pdf', 'svg'], default='png', help='Image format')
    p.add_argument('-s', '--simulations', type=int, default=100, help='Simulations for curve')
    p.add_argument('-r', '--save-r', action='store_true', help='Save R scripts')
    p.add_argument('-n', '--no-plot', action='store_true', help='Skip visualization')
    p.add_argument('-m', '--metadata', help='Optional sample metadata TSV for group summaries')
    p.add_argument('--sample-col', default=None, help='Sample column name in metadata')
    p.add_argument('--group-col', default=None, help='Group column name in metadata')
    p.set_defaults(func=run)


def run(args):
    log("=" * 60)
    log("Running full pipeline: PAV + Curve")
    log("=" * 60)

    # Create a namespace for pav
    class PavArgs:
        pass

    pav_args = PavArgs()
    pav_args.input = args.input
    pav_args.output = args.output
    pav_args.threshold = args.threshold
    pav_args.format = args.format
    pav_args.save_r = args.save_r
    pav_args.no_plot = args.no_plot
    pav_args.stacked = False
    pav_args.func = pav.run

    # Run PAV
    log("\n[Step 1/2] PAV matrix construction...")
    pav.run(pav_args)

    # Run Curve
    if not args.no_plot:
        log("\n[Step 2/2] Saturation curve...")

        class CurveArgs:
            pass

        curve_args = CurveArgs()
        curve_args.pav = os.path.join(args.output, 'pav_matrix.tsv')
        curve_args.output = args.output
        curve_args.simulations = args.simulations
        curve_args.func = curve.run

        if os.path.exists(curve_args.pav):
            curve.run(curve_args)
        else:
            log("Warning: PAV matrix not found, skipping curve")

    if args.metadata:
        log("\n[Optional] Group-level PAV summaries...")

        class GroupArgs:
            pass

        group_args = GroupArgs()
        group_args.pav = os.path.join(args.output, 'pav_matrix.tsv')
        group_args.metadata = args.metadata
        group_args.output = args.output
        group_args.frequency = os.path.join(args.output, 'frequency_table.tsv')
        group_args.sample_col = args.sample_col
        group_args.group_col = args.group_col
        group_args.specific_min = 1.0
        group_args.outside_max = 0.0
        group_args.func = group.run

        if os.path.exists(group_args.pav):
            group.run(group_args)
        else:
            log("Warning: PAV matrix not found, skipping group summaries")

    log("\n" + "=" * 60)
    log("Pipeline complete!")
    log("=" * 60)
