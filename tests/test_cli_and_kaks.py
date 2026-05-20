from pgkit.cli import build_parser
from pgkit.kaks.calculator import estimate_kaks_ng
from pgkit.kaks.sampling import random_sample_pairs


def test_cli_exposes_existing_subcommands() -> None:
    parser = build_parser()
    subcommands = parser._subparsers._group_actions[0].choices

    assert {"run", "pav", "curve", "pie", "bar", "heatmap", "stats", "kaks"} <= set(
        subcommands
    )


def test_kaks_ng_identical_sequences_keep_original_zero_semantics() -> None:
    ka, ks, kaks = estimate_kaks_ng("ATGAAA", "ATGAAA")

    assert ka == 0.0
    assert ks == 0.0
    assert kaks == 0.0


def test_random_sample_pairs_keeps_category_sampling_shape() -> None:
    categories = {"Core": ["OG1"], "Softcore": [], "Dispensable": [], "Private": []}
    og_sequences = {"OG1": {"g1": "M", "g2": "M", "g3": "M"}}

    sampled = random_sample_pairs({}, categories, og_sequences, n_genes=1, n_pairs=2)

    assert len(sampled["Core"]) == 2
    assert {pair[0] for pair in sampled["Core"]} == {"OG1"}
