import itertools
import random

from pgkit.commands.curve import calc_core_pan, iter_sample_sets


def test_calc_core_pan_matches_presence_semantics() -> None:
    pav = {
        "core": {"A": 1, "B": 1},
        "private": {"A": 1, "B": 0},
        "absent": {"A": 0, "B": 0},
    }

    assert calc_core_pan(pav, ("A", "B")) == (1, 2)


def test_iter_sample_sets_exact_mode_does_not_change_small_combinations() -> None:
    species = ["A", "B", "C"]

    assert list(iter_sample_sets(species, 2, simulations=1)) == list(
        itertools.combinations(species, 2)
    )


def test_iter_sample_sets_large_exact_mode_samples_without_materializing_all() -> None:
    species = [f"S{i}" for i in range(20)]
    random.seed(1)

    samples = list(iter_sample_sets(species, 10, simulations=1, max_exact=100))

    assert len(samples) == 200
    assert all(len(sample) == 10 for sample in samples)
