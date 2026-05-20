"""Group-level PAV summaries for plant pan-genome workflows."""

from __future__ import annotations

from collections import defaultdict
from itertools import combinations
from math import comb, isinf
import os
from typing import Dict, Iterable, List, Tuple

from pgkit.core.parser import parse_frequency_table, parse_pav_matrix
from pgkit.core.utils import check_file, ensure_dir, log, write_tsv


def register(subparsers):
    """Register the group subcommand."""
    p = subparsers.add_parser("group", help="Compare PAV frequencies across sample groups")
    p.add_argument("pav", help="PAV matrix file (pav_matrix.tsv)")
    p.add_argument("metadata", help="Sample metadata TSV with sample and group columns")
    p.add_argument("-o", "--output", default="pgkit_output", help="Output directory")
    p.add_argument("-f", "--frequency", help="Optional frequency_table.tsv for category annotation")
    p.add_argument("--sample-col", default=None, help="Sample column name (default: auto/first)")
    p.add_argument("--group-col", default=None, help="Group column name (default: auto/second)")
    p.add_argument(
        "--specific-min",
        type=float,
        default=1.0,
        help="Minimum within-group frequency for group-specific calls (default: 1.0)",
    )
    p.add_argument(
        "--outside-max",
        type=float,
        default=0.0,
        help="Maximum outside-group frequency for group-specific calls (default: 0.0)",
    )
    p.set_defaults(func=run)


def parse_metadata(filepath: str, sample_col=None, group_col=None) -> Dict[str, str]:
    """Parse sample-to-group metadata from a TSV file."""
    with open(filepath, "r", encoding="utf-8") as f:
        rows = [line.rstrip("\n").split("\t") for line in f if line.strip()]

    if not rows:
        raise ValueError(f"{filepath} is empty")

    first = rows[0]
    lower = [x.strip().lower() for x in first]
    has_header = (
        sample_col in first
        or group_col in first
        or "sample" in lower
        or "species" in lower
        or "group" in lower
        or "population" in lower
    )

    if has_header:
        header = first
        data_rows = rows[1:]
        sample_idx = _column_index(header, sample_col, ["Sample", "Species", "sample", "species"], 0)
        group_idx = _column_index(header, group_col, ["Group", "Population", "group", "population"], 1)
    else:
        data_rows = rows
        sample_idx = 0
        group_idx = 1

    mapping = {}
    for line_number, row in enumerate(data_rows, start=2 if has_header else 1):
        if len(row) <= max(sample_idx, group_idx):
            raise ValueError(f"{filepath}:{line_number} does not contain both sample and group")
        sample = row[sample_idx].strip()
        group = row[group_idx].strip()
        if sample and group:
            mapping[sample] = group
    return mapping


def _column_index(header: List[str], requested, candidates: List[str], default: int) -> int:
    if requested:
        if requested not in header:
            raise ValueError(f"Column {requested!r} not found in metadata header")
        return header.index(requested)
    for name in candidates:
        if name in header:
            return header.index(name)
    return default


def build_category_map(frequency_file: str | None) -> Dict[str, str]:
    """Build optional orthogroup-to-category mapping from frequency_table.tsv."""
    if not frequency_file:
        return {}
    return {row["Orthogroup"]: row["Category"] for row in parse_frequency_table(frequency_file)}


def group_samples(species_list: Iterable[str], metadata: Dict[str, str]) -> Dict[str, List[str]]:
    """Return PAV samples grouped by metadata group."""
    grouped = defaultdict(list)
    for sample in species_list:
        group = metadata.get(sample)
        if group:
            grouped[group].append(sample)
    return dict(sorted(grouped.items()))


def frequency_by_group(pav_dict, grouped_samples, category_map):
    """Calculate per-group presence frequency for each orthogroup."""
    rows = []
    for og_id in sorted(pav_dict):
        for group, samples in grouped_samples.items():
            present = sum(pav_dict[og_id].get(sample, 0) for sample in samples)
            total = len(samples)
            freq = present / total if total else 0.0
            rows.append(
                [
                    og_id,
                    category_map.get(og_id, ""),
                    group,
                    present,
                    total,
                    f"{freq:.6f}",
                ]
            )
    return rows


def group_specific_rows(pav_dict, grouped_samples, category_map, specific_min, outside_max):
    """Find orthogroups enriched as exact or thresholded group-specific calls."""
    rows = []
    all_samples = [sample for samples in grouped_samples.values() for sample in samples]
    for og_id in sorted(pav_dict):
        for group, samples in grouped_samples.items():
            outside = [sample for sample in all_samples if sample not in samples]
            present = sum(pav_dict[og_id].get(sample, 0) for sample in samples)
            outside_present = sum(pav_dict[og_id].get(sample, 0) for sample in outside)
            freq = present / len(samples) if samples else 0.0
            outside_freq = outside_present / len(outside) if outside else 0.0
            if freq >= specific_min and outside_freq <= outside_max:
                rows.append(
                    [
                        og_id,
                        category_map.get(og_id, ""),
                        group,
                        present,
                        len(samples),
                        outside_present,
                        len(outside),
                        f"{freq:.6f}",
                        f"{outside_freq:.6f}",
                    ]
                )
    return rows


def pairwise_rows(pav_dict, grouped_samples, category_map):
    """Calculate pairwise group frequency differences with Fisher exact tests."""
    raw_rows = []
    p_values = []
    for og_id in sorted(pav_dict):
        for group_a, group_b in combinations(grouped_samples, 2):
            samples_a = grouped_samples[group_a]
            samples_b = grouped_samples[group_b]
            a_present = sum(pav_dict[og_id].get(sample, 0) for sample in samples_a)
            b_present = sum(pav_dict[og_id].get(sample, 0) for sample in samples_b)
            a_absent = len(samples_a) - a_present
            b_absent = len(samples_b) - b_present
            p_value = fisher_exact_two_sided(a_present, a_absent, b_present, b_absent)
            odds = odds_ratio(a_present, a_absent, b_present, b_absent)
            a_freq = a_present / len(samples_a) if samples_a else 0.0
            b_freq = b_present / len(samples_b) if samples_b else 0.0
            direction = group_a if a_freq > b_freq else group_b if b_freq > a_freq else "Tie"
            raw_rows.append(
                [
                    og_id,
                    category_map.get(og_id, ""),
                    group_a,
                    group_b,
                    a_present,
                    len(samples_a),
                    b_present,
                    len(samples_b),
                    f"{a_freq:.6f}",
                    f"{b_freq:.6f}",
                    _format_float(odds),
                    f"{p_value:.6g}",
                    direction,
                ]
            )
            p_values.append(p_value)

    fdr_values = bh_fdr(p_values)
    rows = []
    for row, fdr in zip(raw_rows, fdr_values):
        rows.append(row[:-1] + [f"{fdr:.6g}", row[-1]])
    return rows


def odds_ratio(a_present: int, a_absent: int, b_present: int, b_absent: int) -> float:
    """Calculate odds ratio with explicit zero-cell handling."""
    if a_absent == 0 or b_present == 0:
        if a_present == 0 or b_absent == 0:
            return 1.0
        return float("inf")
    return (a_present * b_absent) / (a_absent * b_present)


def fisher_exact_two_sided(a_present: int, a_absent: int, b_present: int, b_absent: int) -> float:
    """Two-sided Fisher exact test for a 2x2 table."""
    row1 = a_present + a_absent
    row2 = b_present + b_absent
    col1 = a_present + b_present
    total = row1 + row2
    observed = _hypergeom_prob(a_present, row1, row2, col1, total)
    low = max(0, col1 - row2)
    high = min(row1, col1)
    p_value = 0.0
    for x in range(low, high + 1):
        p = _hypergeom_prob(x, row1, row2, col1, total)
        if p <= observed + 1e-12:
            p_value += p
    return min(1.0, p_value)


def _hypergeom_prob(x: int, row1: int, row2: int, col1: int, total: int) -> float:
    return comb(row1, x) * comb(row2, col1 - x) / comb(total, col1)


def bh_fdr(p_values: List[float]) -> List[float]:
    """Benjamini-Hochberg FDR correction preserving original order."""
    n = len(p_values)
    adjusted = [1.0] * n
    prev = 1.0
    for rank, idx in enumerate(sorted(range(n), key=lambda i: p_values[i]), start=1):
        value = p_values[idx] * n / rank
        adjusted[idx] = value
    for rank, idx in reversed(
        list(enumerate(sorted(range(n), key=lambda i: p_values[i]), start=1))
    ):
        prev = min(prev, adjusted[idx])
        adjusted[idx] = min(prev, 1.0)
    return adjusted


def _format_float(value: float) -> str:
    if isinf(value):
        return "Inf"
    return f"{value:.6g}"


def run(args):
    """Run group-level PAV summaries."""
    ensure_dir(args.output)
    check_file(args.pav, "PAV matrix")
    check_file(args.metadata, "Metadata")
    if args.frequency:
        check_file(args.frequency, "Frequency table")

    pav_dict, species_list = parse_pav_matrix(args.pav)
    metadata = parse_metadata(args.metadata, args.sample_col, args.group_col)
    grouped = group_samples(species_list, metadata)
    category_map = build_category_map(args.frequency)

    if len(grouped) < 1:
        raise ValueError("No metadata samples match the PAV matrix")

    missing = [sample for sample in species_list if sample not in metadata]
    if missing:
        log(f"Warning: {len(missing)} PAV samples have no metadata group and were skipped")

    log(f"Groups: {', '.join(f'{group}={len(samples)}' for group, samples in grouped.items())}")

    group_freq_file = os.path.join(args.output, "group_frequency.tsv")
    write_tsv(
        group_freq_file,
        ["Orthogroup", "Category", "Group", "Present", "Total", "Frequency"],
        frequency_by_group(pav_dict, grouped, category_map),
    )

    specific_file = os.path.join(args.output, "group_specific_orthogroups.tsv")
    write_tsv(
        specific_file,
        [
            "Orthogroup",
            "Category",
            "Group",
            "Present",
            "Total",
            "Outside_Present",
            "Outside_Total",
            "Frequency",
            "Outside_Frequency",
        ],
        group_specific_rows(
            pav_dict,
            grouped,
            category_map,
            args.specific_min,
            args.outside_max,
        ),
    )

    pairwise_file = os.path.join(args.output, "group_pairwise.tsv")
    write_tsv(
        pairwise_file,
        [
            "Orthogroup",
            "Category",
            "Group_A",
            "Group_B",
            "A_Present",
            "A_Total",
            "B_Present",
            "B_Total",
            "A_Frequency",
            "B_Frequency",
            "Odds_Ratio",
            "P_Value",
            "FDR",
            "Direction",
        ],
        pairwise_rows(pav_dict, grouped, category_map),
    )

    log("Output files:")
    log(f"  Group frequencies: {group_freq_file}")
    log(f"  Group-specific orthogroups: {specific_file}")
    log(f"  Pairwise group tests: {pairwise_file}")
    log("Done!")
