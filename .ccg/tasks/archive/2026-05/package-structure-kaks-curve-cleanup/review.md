# Review

## Verification

- `pytest -q --basetemp=F:\tmp\pgkit_pytest`: 9 passed.
- `python pgkit.py --help`: passed.
- `python -m pgkit --help`: passed.
- `python pgkit.py kaks --help`: passed.
- `python -m pgkit pav --help`: passed.
- Small `pav --no-plot` smoke test on a synthetic Orthogroups.tsv produced expected Core/Dispensable/Private outputs.

## Scope and compatibility

- Valid-input command semantics, defaults, output filenames, and output columns are preserved.
- Root `pgkit.py`, legacy `src/*.py`, and legacy `lib/*.py` remain as compatibility wrappers.
- Ka/Ks computation code was mechanically split by responsibility; formulas and sampling defaults were not changed.
- Curve exact-combination mode now avoids materializing huge combinations before fallback sampling while preserving the original exact-or-random policy.

## Residual notes

- `.omx/` was already untracked before this task and was not touched.
- A read-only review agent was started but did not return before local validation completed.
