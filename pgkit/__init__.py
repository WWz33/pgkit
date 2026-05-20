"""Source checkout shim for ``python -m pgkit``.

The installable package lives in ``src/pgkit``. This shim keeps direct
checkout execution working without changing runtime package metadata.
"""

from __future__ import annotations

from pathlib import Path

_SRC_PKG = Path(__file__).resolve().parents[1] / "src" / "pgkit"
if _SRC_PKG.exists():
    __path__.append(str(_SRC_PKG))

__all__ = ["__version__"]
__version__ = "0.1.0"
