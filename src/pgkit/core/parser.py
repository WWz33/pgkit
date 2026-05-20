"""Parsers for OrthoFinder and pgkit tabular outputs."""

from __future__ import annotations

from typing import List


def _read_header(handle, filepath: str) -> List[str]:
    header = handle.readline().rstrip("\n").split("\t")
    if len(header) < 2:
        raise ValueError(f"{filepath} must contain an ID column and at least one sample column")
    return header


def parse_orthogroups_tsv(filepath: str):
    """Parse an OrthoFinder ``Orthogroups.tsv`` file."""
    orthogroups = {}

    with open(filepath, "r", encoding="utf-8") as f:
        header = _read_header(f, filepath)
        species_list = header[1:]

        for line_number, line in enumerate(f, start=2):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if not parts[0]:
                raise ValueError(f"{filepath}:{line_number} has an empty orthogroup ID")

            og_data = {}
            for i, sp in enumerate(species_list):
                val = parts[i + 1].strip() if i + 1 < len(parts) else ""
                og_data[sp] = val if val and val != "-" else ""
            orthogroups[parts[0]] = og_data

    return orthogroups, species_list


def parse_unassigned_genes(filepath: str):
    """Parse an OrthoFinder ``Orthogroups_UnassignedGenes.tsv`` file."""
    return parse_orthogroups_tsv(filepath)


def parse_pav_matrix(filepath: str):
    """Parse a pgkit PAV matrix."""
    pav = {}

    with open(filepath, "r", encoding="utf-8") as f:
        header = _read_header(f, filepath)
        species_list = header[1:]

        for line_number, line in enumerate(f, start=2):
            parts = line.rstrip("\n").split("\t")
            if not parts or not parts[0]:
                continue
            pav[parts[0]] = {}
            for i, sp in enumerate(species_list):
                raw_value = parts[i + 1] if i + 1 < len(parts) else "0"
                try:
                    pav[parts[0]][sp] = int(raw_value)
                except ValueError as exc:
                    raise ValueError(
                        f"{filepath}:{line_number} contains non-integer PAV value {raw_value!r}"
                    ) from exc

    return pav, species_list


def parse_frequency_table(filepath: str):
    """Parse a pgkit frequency table."""
    data = []
    with open(filepath, "r", encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        required = {"Orthogroup", "Frequency", "Category", "Species_Count"}
        missing = required - set(header)
        if missing:
            raise ValueError(f"{filepath} is missing required columns: {sorted(missing)}")

        for line_number, line in enumerate(f, start=2):
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            row = {}
            for i, col in enumerate(header):
                value = parts[i] if i < len(parts) else ""
                if col == "Frequency":
                    row[col] = float(value)
                elif col == "Species_Count":
                    row[col] = int(value)
                else:
                    row[col] = value
            data.append(row)
    return data
