# Script Style Guide

The seven rules every script in `skills/<name>/scripts/` must satisfy.

---

## 1. Module docstring

A module-level docstring at the top of every file stating:
- What the script does (one sentence).
- Inputs (CLI args or function args).
- Outputs (files produced).
- An example invocation.

Example:

```python
"""Validate a tabular dataset before analysis.

Inputs:
    --input PATH        Path to a CSV/TSV/parquet file.
    --id-col COL        Name of the ID column.
    --output-dir DIR    Where to write the validation report.

Outputs:
    <output_dir>/validation_report.json
    <output_dir>/validation_summary.md

Example:
    python validate_data.py --input data.csv --id-col sample_id \
        --output-dir out/
"""
```

Why: the agent reads this when picking which script to run.

---

## 2. Type hints on every public function

```python
def adjust_pvalues(
    pvalues: pd.Series | np.ndarray,
    method: str = "bh",
    alpha: float = 0.05,
) -> tuple[np.ndarray, np.ndarray]:
    ...
```

Why: type hints document the contract, enable static checking, and reduce ambiguity for the next reader.

---

## 3. `argparse` with `--help` for CLI scripts

Every script that's intended to be invoked from the command line must have:

```python
def main() -> None:
    parser = argparse.ArgumentParser(description="What this does.")
    parser.add_argument("--input", required=True, help="Input file path.")
    parser.add_argument("--output", default="out.csv", help="Output path.")
    args = parser.parse_args()
    ...

if __name__ == "__main__":
    main()
```

Why: `--help` is the discoverable entry point. No CLI script should require reading the source to figure out how to invoke it.

For pure library scripts (R or Python with no CLI mode), provide a `main()` function or top-level functions with clear docstrings.

---

## 4. Success verification line

After completing the main work, print a line that includes:
- A ✓ marker.
- What the script did.
- A quantitative outcome (rows, paths, counts).

```python
print(f"✓ Validation completed: {len(df)} rows / {output_path}")
print(f"✓ Multiple testing completed: {n_sig}/{len(df)} significant")
```

Why: the calling agent uses this line to verify the step succeeded.

---

## 5. Graceful ImportError

If a script needs a third-party package that isn't in the standard library, catch the import error and exit with an actionable message:

```python
try:
    import pandas as pd
except ImportError:
    sys.exit("Install pandas with: pip install pandas")
```

Why: opaque `ModuleNotFoundError` traces are user-hostile. The agent should know exactly what to install.

---

## 6. `if __name__ == "__main__":` guard

Always wrap the CLI entry in the guard:

```python
if __name__ == "__main__":
    main()
```

Why: lets the file be imported as a module without running the CLI.

---

## 7. No hardcoded absolute paths

Bad:

```python
df = pd.read_csv("/Users/me/data.csv")
df.to_csv("/Users/me/out/result.csv")
```

Good:

```python
def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    df = pd.read_csv(args.input)
    df.to_csv(args.output, index=False)
```

Why: hardcoded paths break the moment someone else clones the repo.

---

## R-specific notes

R scripts follow the spirit of the rules with R idioms:

- Replace module docstring with a `#` comment block at the top.
- Use `commandArgs(trailingOnly = TRUE)` or `optparse::OptionParser` for CLI.
- `tryCatch(library(pkg), error = function(e) stop("Install pkg with: ..."))`.
- `if (sys.nframe() == 0L) main()` is the equivalent of the `__main__` guard.

Example R header:

```r
# script_name.R
# What this does.
#
# Inputs:
#   - arg1: description
# Outputs:
#   - file1: description
# Example:
#   Rscript script_name.R --input data.csv --output out.csv
```

---

## Anti-patterns to avoid

- Reading a config file from a fixed location.
- Implicit dependence on the current working directory.
- Catching all exceptions and printing "something went wrong" — be specific.
- Print statements that don't include a status marker (✓, ⚠, ✗).
- Skipping the `--help` because "it's just a small script".
- One giant `main()` function — break into named helpers.
- Mutable global state.

---

## When a rule conflicts with the task

The rules describe the default. If a script genuinely needs to deviate (e.g., a pure helper module with no CLI), follow the spirit:

- Document what the script does.
- Make it usable without reading the source.
- Be explicit about contracts.

If you find yourself repeatedly violating a rule for a class of script, propose a refinement to this style guide rather than silently breaking it.
