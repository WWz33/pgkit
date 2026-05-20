"""Shared utilities for pgkit commands."""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
import os
import sys
from typing import List, Union


PathInput = Union[str, os.PathLike]


def ensure_dir(path: PathInput) -> str:
    """Create a directory if needed and return it as a string path."""
    os.makedirs(path, exist_ok=True)
    return str(path)


def log(msg: str) -> None:
    """Print a timestamped status message."""
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] {msg}")


def check_file(filepath: PathInput, desc: str = "file") -> None:
    """Exit with a clear error if a required file does not exist."""
    if not os.path.exists(filepath):
        print(f"Error: {desc} does not exist: {filepath}")
        sys.exit(1)


def write_tsv(filepath, header, rows) -> None:
    """Write rows to a UTF-8 TSV file."""
    with open(filepath, "w", encoding="utf-8") as f:
        f.write("\t".join(str(x) for x in header) + "\n")
        for row in rows:
            f.write("\t".join(str(x) for x in row) + "\n")


def read_lines(filepath) -> List[str]:
    """Read non-empty stripped lines from a UTF-8 text file."""
    with open(filepath, "r", encoding="utf-8") as f:
        return [line.strip() for line in f if line.strip()]


def script_path(script_name: str) -> str:
    """Return the absolute path to a bundled R script."""
    return str(Path(__file__).resolve().parents[1] / "scripts" / script_name)
