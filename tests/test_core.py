from pathlib import Path

from pgkit.core.classify import (
    build_og_to_category,
    classify_orthogroups,
    count_genes_per_species_per_category,
)
from pgkit.core.parser import parse_frequency_table, parse_orthogroups_tsv, parse_pav_matrix


def test_parse_orthogroups_preserves_original_gene_strings(tmp_path: Path) -> None:
    orthogroups_file = tmp_path / "Orthogroups.tsv"
    orthogroups_file.write_text(
        "Orthogroup\tA\tB\n"
        "OG0001\tgene1, gene2\t-\n"
        "OG0002\t\tgene3\n",
        encoding="utf-8",
    )

    orthogroups, species = parse_orthogroups_tsv(str(orthogroups_file))

    assert species == ["A", "B"]
    assert orthogroups == {
        "OG0001": {"A": "gene1, gene2", "B": ""},
        "OG0002": {"A": "", "B": "gene3"},
    }


def test_classify_and_gene_counts_keep_category_semantics() -> None:
    orthogroups = {
        "core": {"A": "a1", "B": "b1", "C": "c1"},
        "soft": {"A": "a2", "B": "b2", "C": ""},
        "disp": {"A": "a3", "B": "", "C": "c3"},
        "priv": {"A": "", "B": "b4, b5", "C": ""},
    }
    species = ["A", "B", "C"]

    categories, species_counts = classify_orthogroups(orthogroups, species, 0.66)
    mapping = build_og_to_category(categories)
    counts = count_genes_per_species_per_category(orthogroups, categories, species)

    assert categories == {
        "Core": ["core"],
        "Softcore": ["soft", "disp"],
        "Dispensable": [],
        "Private": ["priv"],
    }
    assert species_counts == {"core": 3, "soft": 2, "disp": 2, "priv": 1}
    assert mapping["core"] == "Core"
    assert counts["B"]["Private"] == 2


def test_parse_pav_and_frequency_tables(tmp_path: Path) -> None:
    pav_file = tmp_path / "pav_matrix.tsv"
    pav_file.write_text("Orthogroup\tA\tB\nOG0001\t1\t0\n", encoding="utf-8")
    freq_file = tmp_path / "frequency_table.tsv"
    freq_file.write_text(
        "Orthogroup\tFrequency\tCategory\tSpecies_Count\n"
        "OG0001\t0.5000\tDispensable\t1\n",
        encoding="utf-8",
    )

    pav, species = parse_pav_matrix(str(pav_file))
    freq = parse_frequency_table(str(freq_file))

    assert species == ["A", "B"]
    assert pav == {"OG0001": {"A": 1, "B": 0}}
    assert freq == [
        {
            "Orthogroup": "OG0001",
            "Frequency": 0.5,
            "Category": "Dispensable",
            "Species_Count": 1,
        }
    ]
