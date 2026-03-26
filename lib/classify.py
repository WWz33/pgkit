"""
Gene family classification
"""
import math


def classify_orthogroups(orthogroups, species_list, soft_core_threshold=0.9):
    """
    Classify orthogroups by presence frequency
    
    Categories:
        - Core: 100% samples present
        - Softcore: >=threshold (default 90%) samples present
        - Dispensable: >1 sample but <threshold
        - Specific: only 1 sample present
    
    Returns:
        categories: dict {category: [og_ids]}
        species_counts: dict {og_id: count}
    """
    total = len(species_list)
    soft_core_min = math.ceil(total * soft_core_threshold)

    categories = {'Core': [], 'Softcore': [], 'Dispensable': [], 'Specific': []}
    species_counts = {}

    for og_id, og_data in orthogroups.items():
        count = sum(1 for v in og_data.values() if v.strip())
        species_counts[og_id] = count

        if count == total:
            categories['Core'].append(og_id)
        elif count >= soft_core_min:
            categories['Softcore'].append(og_id)
        elif count == 1:
            categories['Specific'].append(og_id)
        else:
            categories['Dispensable'].append(og_id)

    return categories, species_counts


def build_og_to_category(categories):
    """构建orthogroup到分类的映射"""
    mapping = {}
    for cat, ogs in categories.items():
        for og in ogs:
            mapping[og] = cat
    return mapping


def count_genes_per_species_per_category(orthogroups, categories, species_list):
    """
    统计每个物种在每个类别中的基因数量
    
    返回:
        dict {species: {category: gene_count}}
    """
    result = {sp: {cat: 0 for cat in categories} for sp in species_list}

    for cat, ogs in categories.items():
        for og_id in ogs:
            if og_id in orthogroups:
                og_data = orthogroups[og_id]
                for sp in species_list:
                    genes_str = og_data.get(sp, '')
                    if genes_str and genes_str.strip():
                        genes = [g.strip() for g in genes_str.split(',')]
                        result[sp][cat] += len(genes)

    return result
