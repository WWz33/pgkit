"""
解析OrthoFinder输出文件
"""


def parse_orthogroups_tsv(filepath):
    """
    解析Orthogroups.tsv文件
    
    返回:
        orthogroups: dict {og_id: {species: genes_str}}
        species_list: list of species names
    """
    orthogroups = {}
    species_list = []

    with open(filepath, 'r', encoding='utf-8') as f:
        header = f.readline().strip().split('\t')
        species_list = header[1:]

        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            og_id = parts[0]
            og_data = {}
            for i, sp in enumerate(species_list):
                val = parts[i + 1].strip() if i + 1 < len(parts) else ''
                og_data[sp] = val if val and val != '-' else ''
            orthogroups[og_id] = og_data

    return orthogroups, species_list


def parse_unassigned_genes(filepath):
    """
    解析Orthogroups_UnassignedGenes.tsv文件
    
    返回:
        dict {og_id: {species: genes_str}}
    """
    unassigned = {}
    species_list = []

    with open(filepath, 'r', encoding='utf-8') as f:
        header = f.readline().strip().split('\t')
        species_list = header[1:]

        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            og_id = parts[0]
            og_data = {}
            for i, sp in enumerate(species_list):
                val = parts[i + 1].strip() if i + 1 < len(parts) else ''
                og_data[sp] = val if val and val != '-' else ''
            unassigned[og_id] = og_data

    return unassigned, species_list


def parse_pav_matrix(filepath):
    """
    解析PAV矩阵文件
    
    返回:
        pav: dict {og_id: {species: 0/1}}
        species_list: list
    """
    pav = {}
    species_list = []

    with open(filepath, 'r', encoding='utf-8') as f:
        header = f.readline().strip().split('\t')
        species_list = header[1:]

        for line in f:
            parts = line.strip().split('\t')
            og_id = parts[0]
            pav[og_id] = {}
            for i, sp in enumerate(species_list):
                pav[og_id][sp] = int(parts[i + 1]) if i + 1 < len(parts) else 0

    return pav, species_list


def parse_frequency_table(filepath):
    """
    解析频率表文件
    
    返回:
        list of dict
    """
    data = []
    with open(filepath, 'r', encoding='utf-8') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            parts = line.strip().split('\t')
            row = {}
            for i, col in enumerate(header):
                if col == 'Frequency':
                    row[col] = float(parts[i])
                elif col == 'Species_Count':
                    row[col] = int(parts[i])
                else:
                    row[col] = parts[i]
            data.append(row)
    return data
