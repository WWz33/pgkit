"""Sequence loading helpers for Ka/Ks calculations."""

import glob
import os

from Bio import SeqIO

from pgkit.kaks.constants import CODON_TABLES

# ============================================================
# FASTA parsing
# ============================================================
def parse_fasta(filepath):
    """Parse FASTA file"""
    sequences = {}
    for record in SeqIO.parse(filepath, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences


# ============================================================
# CDS validation (length % 3 == 0)
# ============================================================
def load_cds_file(cds_file, genetic_code=1):
    """
    Load and validate CDS sequences

    Validation rules:
    1. Length must be divisible by 3
    2. No internal stop codons
    3. Minimum 30 bp
    4. Valid nucleotides only (A, C, G, T, N)
    """
    cds_dict = {}
    invalid = []
    codon_table = CODON_TABLES.get(genetic_code, CODON_TABLES[1])

    for record in SeqIO.parse(cds_file, "fasta"):
        gene_id = record.id
        seq = str(record.seq).upper()

        # Rule 1: Length divisible by 3
        if len(seq) % 3 != 0:
            invalid.append((gene_id, f"Length {len(seq)} not divisible by 3"))
            continue

        # Rule 2: Minimum length (30 bp = 10 codons)
        if len(seq) < 30:
            invalid.append((gene_id, f"Too short ({len(seq)} bp)"))
            continue

        # Rule 3: No internal stop codons
        has_internal_stop = False
        for i in range(0, len(seq) - 3, 3):
            codon = seq[i:i+3]
            aa = codon_table.get(codon, 'X')
            if aa == '*':
                has_internal_stop = True
                break

        if has_internal_stop:
            invalid.append((gene_id, "Internal stop codon"))
            continue

        # Rule 4: Valid nucleotides
        invalid_chars = set(seq) - set('ACGTN')
        if invalid_chars:
            invalid.append((gene_id, f"Invalid characters: {invalid_chars}"))
            continue

        cds_dict[gene_id] = seq

    return cds_dict, invalid


# ============================================================
# Protein sequence loading
# ============================================================
def load_orthogroup_sequences(seq_dir):
    """Load protein sequences from Orthogroup_Sequences directory"""
    protein_dict = {}
    for fa_file in glob.glob(os.path.join(seq_dir, "*.fa")):
        for record in SeqIO.parse(fa_file, "fasta"):
            protein_dict[record.id] = str(record.seq)
    return protein_dict


def parse_orthogroup_sequences_with_og(seq_dir):
    """Load protein sequences with orthogroup mapping"""
    og_sequences = {}
    for fa_file in glob.glob(os.path.join(seq_dir, "*.fa")):
        og_id = os.path.splitext(os.path.basename(fa_file))[0]
        genes = {}
        for record in SeqIO.parse(fa_file, "fasta"):
            genes[record.id] = str(record.seq)
        og_sequences[og_id] = genes
    return og_sequences
