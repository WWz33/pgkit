from pathlib import Path

from pgkit.commands.group import (
    fisher_exact_two_sided,
    group_samples,
    parse_metadata,
    run,
)


class Args:
    pass


def test_parse_metadata_accepts_header_and_auto_columns(tmp_path: Path) -> None:
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "Sample\tGroup\n"
        "A\tWild\n"
        "B\tCultivar\n",
        encoding="utf-8",
    )

    assert parse_metadata(str(metadata)) == {"A": "Wild", "B": "Cultivar"}


def test_group_samples_skips_samples_without_metadata() -> None:
    grouped = group_samples(["A", "B", "C"], {"A": "Wild", "B": "Cultivar"})

    assert grouped == {"Cultivar": ["B"], "Wild": ["A"]}


def test_fisher_exact_two_sided_detects_extreme_split() -> None:
    assert fisher_exact_two_sided(2, 0, 0, 2) == 1 / 3


def test_group_command_writes_frequency_specific_and_pairwise_tables(tmp_path: Path) -> None:
    pav = tmp_path / "pav_matrix.tsv"
    pav.write_text(
        "Orthogroup\tA\tB\tC\tD\n"
        "OG1\t1\t1\t0\t0\n"
        "OG2\t1\t1\t1\t1\n",
        encoding="utf-8",
    )
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text(
        "Sample\tGroup\n"
        "A\tWild\n"
        "B\tWild\n"
        "C\tCultivar\n"
        "D\tCultivar\n",
        encoding="utf-8",
    )
    frequency = tmp_path / "frequency_table.tsv"
    frequency.write_text(
        "Orthogroup\tFrequency\tCategory\tSpecies_Count\n"
        "OG1\t0.5000\tDispensable\t2\n"
        "OG2\t1.0000\tCore\t4\n",
        encoding="utf-8",
    )

    args = Args()
    args.pav = str(pav)
    args.metadata = str(metadata)
    args.frequency = str(frequency)
    args.output = str(tmp_path)
    args.sample_col = None
    args.group_col = None
    args.specific_min = 1.0
    args.outside_max = 0.0

    run(args)

    assert (tmp_path / "group_frequency.tsv").exists()
    pairwise = (tmp_path / "group_pairwise.tsv").read_text(encoding="utf-8")
    specific = (tmp_path / "group_specific_orthogroups.tsv").read_text(encoding="utf-8")
    assert "OG1\tDispensable" in pairwise
    assert "OG1\tDispensable\tWild" in specific
