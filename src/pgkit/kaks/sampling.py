"""Sampling helpers for pan-genome Ka/Ks calculations."""

from collections import defaultdict
import random

# ============================================================
# Random sampling
# ============================================================
def random_sample_pairs(orthogroups, categories, og_sequences, n_genes, n_pairs):
    """Random sample gene pairs from each category"""
    sampled = defaultdict(list)

    for category, og_list in categories.items():
        valid_ogs = [og for og in og_list if og in og_sequences and len(og_sequences[og]) >= 2]
        if not valid_ogs:
            continue

        n_sample = min(n_genes, len(valid_ogs))
        for og_id in random.sample(valid_ogs, n_sample):
            genes = list(og_sequences[og_id].keys())
            all_pairs = [(genes[i], genes[j]) for i in range(len(genes)) for j in range(i+1, len(genes))]
            if not all_pairs:
                continue
            n_sample_pairs = min(n_pairs, len(all_pairs))
            for g1, g2 in random.sample(all_pairs, n_sample_pairs):
                sampled[category].append((og_id, g1, g2))

    return sampled
