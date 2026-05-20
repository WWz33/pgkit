# Plan

## Success criteria

- `python -m pgkit --help` works as the primary package entrypoint.
- Existing `python pgkit.py --help` remains compatible.
- Command modules live under a real `pgkit` package with explicit internal imports.
- `kaks.py` is split by responsibility without changing the public `pgkit kaks` CLI.
- Curve sampling avoids materializing huge combination lists and exposes deterministic seeding.
- Parser and utility modules have readable English docstrings and fail loudly on malformed inputs.
- Minimal pytest coverage verifies parser/classifier behavior, curve sampling, and CLI compatibility.

## Implementation steps

1. Create package layout:
   - `pgkit/__main__.py`, `pgkit/cli.py`
   - `pgkit/commands/`
   - `pgkit/core/`
   - `pgkit/kaks/`
   - `pgkit/scripts/`
2. Move existing command modules into `pgkit/commands` and helpers into `pgkit/core`.
3. Split Ka/Ks functions into modules:
   - `pgkit/kaks/io.py`
   - `pgkit/kaks/alignment.py`
   - `pgkit/kaks/calculator.py`
   - `pgkit/kaks/sampling.py`
   - `pgkit/kaks/report.py`
   - `pgkit/commands/kaks.py`
4. Update imports and script paths to package-relative paths.
5. Add `pyproject.toml`, pytest tests, and fix README command examples.
6. Run CLI smoke tests and pytest.
