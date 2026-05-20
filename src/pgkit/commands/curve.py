"""
Saturation curve analysis
Core/Pan gene families vs. number of accessions
"""
import sys
import os
import itertools
import math
import random
import subprocess

from pgkit.core.parser import parse_pav_matrix
from pgkit.core.utils import check_file, ensure_dir, log, script_path


def register(subparsers):
    p = subparsers.add_parser('curve', help='Core/Pan gene family saturation curve')
    p.add_argument('pav', help='PAV matrix file (pav_matrix.tsv)')
    p.add_argument('-o', '--output', default='pgkit_output', help='Output directory')
    p.add_argument('-s', '--simulations', type=int, default=100,
                   help='Number of simulations per sample count (default: 100)')
    p.set_defaults(func=run)


def calc_core_pan(pav_dict, sample_subset):
    """Calculate core and pan counts for a subset of samples"""
    n = len(sample_subset)
    core = 0
    pan = 0
    for og_id, og_data in pav_dict.items():
        presence = sum(og_data.get(sp, 0) for sp in sample_subset)
        if presence == n:
            core += 1
        if presence >= 1:
            pan += 1
    return core, pan


def iter_sample_sets(species_list, k, simulations, max_exact=10000):
    """Yield sample subsets using the original exact-or-random policy."""
    if simulations <= 1:
        n_combinations = math.comb(len(species_list), k)
        if n_combinations > max_exact:
            log(f"  n={k}: {n_combinations} combinations too many, using random sampling")
            for _ in range(200):
                yield tuple(random.sample(species_list, k))
        else:
            yield from itertools.combinations(species_list, k)
    else:
        for _ in range(simulations):
            yield tuple(random.sample(species_list, k))


def run(args):
    ensure_dir(args.output)

    log("Reading PAV matrix...")
    check_file(args.pav, "PAV matrix")
    pav_dict, species_list = parse_pav_matrix(args.pav)
    n = len(species_list)

    log(f"Species: {n}, Orthogroups: {len(pav_dict)}")

    results = []

    for k in range(1, n + 1):
        core_vals = []
        pan_vals = []
        for combo in iter_sample_sets(species_list, k, args.simulations):
            c, p = calc_core_pan(pav_dict, combo)
            core_vals.append(c)
            pan_vals.append(p)

        core_avg = sum(core_vals) / len(core_vals)
        pan_avg = sum(pan_vals) / len(pan_vals)

        core_sd = (sum((v - core_avg) ** 2 for v in core_vals) / len(core_vals)) ** 0.5
        pan_sd = (sum((v - pan_avg) ** 2 for v in pan_vals) / len(pan_vals)) ** 0.5

        results.append((k, core_avg, pan_avg, core_sd, pan_sd))
        log(f"  n={k}: core={core_avg:.1f}+/-{core_sd:.1f}, pan={pan_avg:.1f}+/-{pan_sd:.1f}")

    # Save data with SD
    data_file = os.path.join(args.output, 'saturation_curve.tsv')
    with open(data_file, 'w') as f:
        f.write('n_accessions\tcore\tpan\tcore_sd\tpan_sd\n')
        for k, core, pan, csd, psd in results:
            f.write(f'{k}\t{core:.1f}\t{pan:.1f}\t{csd:.1f}\t{psd:.1f}\n')
    log(f"Data saved: {data_file}")

    # Call R script
    r_script = script_path('plot_curve.R')
    out_prefix = os.path.join(args.output, 'curve')

    log("Generating saturation curve plot...")
    try:
        result = subprocess.run(
            ['Rscript', os.path.abspath(r_script), data_file, out_prefix],
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
        log(f"Manual run: Rscript {os.path.abspath(r_script)} {data_file} {out_prefix}")
