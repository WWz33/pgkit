"""Alignment helpers for Ka/Ks calculations."""

import os
import subprocess

from Bio import SeqIO

# Protein alignment (MUSCLE / MAFFT / fallback)
# ============================================================
def align_proteins(gene1, gene2, prot1, prot2, tmp_dir):
    """Align two protein sequences"""
    pair_dict = {gene1: prot1, gene2: prot2}
    tmp_file = os.path.join(tmp_dir, f"{gene1}_{gene2}.fa")
    out_file = os.path.join(tmp_dir, f"{gene1}_{gene2}_aln.fa")

    with open(tmp_file, 'w') as f:
        f.write(f">{gene1}\n{prot1}\n>{gene2}\n{prot2}\n")

    aligned = None

    for tool, cmd in [
        ('muscle', ['muscle', '-in', tmp_file, '-out', out_file]),
        ('mafft', ['mafft', '--auto', '--quiet', tmp_file]),
    ]:
        try:
            if tool == 'mafft':
                result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
                if result.returncode == 0:
                    with open(out_file, 'w') as f:
                        f.write(result.stdout)
            else:
                subprocess.run(cmd, capture_output=True, text=True, timeout=30)

            if os.path.exists(out_file):
                aligned = {}
                for record in SeqIO.parse(out_file, "fasta"):
                    aligned[record.id] = str(record.seq)
                break
        except (FileNotFoundError, subprocess.TimeoutExpired):
            continue

    # Cleanup
    for f in [tmp_file, out_file]:
        if os.path.exists(f):
            os.remove(f)

    # Fallback: simple padding alignment
    if aligned is None:
        len1, len2 = len(prot1), len(prot2)
        if len1 < len2:
            prot1 = prot1 + '-' * (len2 - len1)
        elif len2 < len1:
            prot2 = prot2 + '-' * (len1 - len2)
        aligned = {gene1: prot1, gene2: prot2}

    return aligned


# ============================================================
# Back-translation
# ============================================================
def back_translate(protein_alignment, cds_dict):
    """Back-translate protein alignment to CDS alignment"""
    cds_alignment = {}

    for gene_id, prot_seq in protein_alignment.items():
        if gene_id not in cds_dict:
            continue

        cds_seq = cds_dict[gene_id]
        aligned_cds = []
        cds_pos = 0

        for aa in prot_seq:
            if aa == '-':
                aligned_cds.append('---')
            elif cds_pos + 3 <= len(cds_seq):
                aligned_cds.append(cds_seq[cds_pos:cds_pos+3])
                cds_pos += 3
            else:
                aligned_cds.append('---')

        cds_alignment[gene_id] = ''.join(aligned_cds)

    return cds_alignment


# ============================================================
# Generate AXT format (for KaKs_Calculator)
# ============================================================
def generate_axt(gene1, gene2, cds1, cds2, output_file):
    """Generate AXT format file for KaKs_Calculator"""
    with open(output_file, 'w') as f:
        f.write(f"{gene1} {gene2}\n")
        f.write(f"{cds1}\n")
        f.write(f"{cds2}\n\n")
    return output_file
