# Reproducibility — Seeds, Versions, Environments

What "reproducible" actually means in practice: given the same data and code, you get the same numerical result, today and in two years.

---

## The three failure modes

1. **Non-deterministic code.** A random-number step without a seed. Even one unseeded `random.choice` can change everything downstream.
2. **Version drift.** A package update changes default parameters or numerical algorithms. Yesterday's `sklearn==1.2` is not today's `sklearn==1.3`.
3. **Environment differences.** macOS vs Linux, x86 vs ARM, BLAS implementation — subtle but real.

Mitigate all three.

---

## Random seeds

Set every relevant generator at the top of every entry point:

```python
import os, random, numpy as np
SEED = 42
os.environ["PYTHONHASHSEED"] = str(SEED)
random.seed(SEED)
np.random.seed(SEED)

# If using PyTorch
import torch
torch.manual_seed(SEED)
torch.cuda.manual_seed_all(SEED)
torch.use_deterministic_algorithms(True, warn_only=True)

# If using TensorFlow
import tensorflow as tf
tf.random.set_seed(SEED)
tf.keras.utils.set_random_seed(SEED)

# If using scikit-learn estimators that accept random_state, pass it
```

```r
set.seed(42)
```

Common omissions:
- Cross-validation splits (`random_state=` argument).
- Train-test splits.
- Bootstrap resampling.
- Dropout / data augmentation in neural networks.
- KMeans initialization.
- Imputation methods with stochastic components.

Save seeds to disk:

```python
import json
with open("seeds.json", "w") as fh:
    json.dump({"seed": SEED}, fh, indent=2)
```

---

## Version pinning

### Python — requirements.txt and lockfiles

After installing what you need:

```bash
pip freeze > requirements.txt
```

This pins every package, including transitive dependencies. Reading a `requirements.txt` later and `pip install -r` it gives the same environment (modulo OS / Python version).

For a stricter contract, use a lockfile (`pip-tools`, `poetry`, `uv`):

```bash
pip-compile pyproject.toml > requirements.lock
```

### Python version

Pin the interpreter too:

```
# .python-version (used by pyenv)
3.11.6
```

Or in `pyproject.toml`:

```toml
[project]
requires-python = ">=3.11,<3.12"
```

### Conda / mamba environments

```bash
conda env export --no-builds > environment.yml
```

`--no-builds` strips the build hash, making the file portable across platforms.

### R — renv

```r
renv::init()       # at project creation
renv::snapshot()   # save current package versions
renv::restore()    # restore from the snapshot
```

Commit `renv.lock` to version control.

---

## Containers

For multi-language workflows or strict OS-level reproducibility, package the analysis in a container.

```dockerfile
FROM python:3.11-slim
COPY requirements.txt /tmp/
RUN pip install --no-deps -r /tmp/requirements.txt
COPY src /app/src
WORKDIR /app
ENTRYPOINT ["python", "-m", "src.main"]
```

Pin the base image with a content hash (`python:3.11-slim@sha256:...`) for the strongest reproducibility guarantee.

---

## Data versioning

Code reproducibility is half the story. Data reproducibility is the other half.

- Hash input files (`sha256sum data/raw/*.csv > data/raw/hashes.txt`).
- Use DVC, git-lfs, or an immutable artifact store for large files.
- Never modify raw data in place — keep raw separate from processed.

---

## Workflow checklist

- [ ] Seeds set at every stochastic step.
- [ ] `requirements.txt` (or `environment.yml` or `renv.lock`) committed.
- [ ] Python (or R) version pinned.
- [ ] Raw data hashed; processed data has a clear lineage from raw.
- [ ] Analysis runs end-to-end from a clean clone (a "fresh git clone" smoke test).
- [ ] Results saved with a manifest of (code commit, data hash, seed, environment).

---

## Common pitfalls

- **`numpy.random.seed` is global state** — for parallel code, use `numpy.random.default_rng(SEED)` per-worker.
- **`PYTHONHASHSEED` matters** for dictionary iteration order in older Python; setting it gives stable iteration even if you accidentally rely on it.
- **GPU non-determinism** — some CUDA operations are non-deterministic by design. `torch.use_deterministic_algorithms(True)` reveals which.
- **Locale and encoding** — `LC_ALL=C` for stable sorting; explicit `encoding="utf-8"` on file reads.
- **Floating-point summation order** — parallel reductions may yield bit-different sums. Usually negligible, occasionally important.

---

## References

- Sandve GK, Nekrutenko A, Taylor J, Hovig E. *Ten simple rules for reproducible computational research.* PLoS Comp Biol 2013;9:e1003285.
- Wilson G, Bryan J, Cranston K, Kitzes J, Nederbragt L, Teal TK. *Good enough practices in scientific computing.* PLoS Comp Biol 2017;13:e1005510.
- Stodden V, McNutt M, Bailey DH, et al. *Enhancing reproducibility for computational methods.* Science 2016;354:1240–1241.
