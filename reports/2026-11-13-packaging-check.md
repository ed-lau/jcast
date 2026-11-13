# Packaging check: Python versions and dependency bounds (PRs 10, 11)

| | |
|---|---|
| **Date** | 2026-11-13 |
| **Status** | Decisions taken 2026-11-13 (§5). **Implemented in PR 10** (`4a039a9`, §2); **CI in PR 11** (`dee41b6`, §6). No output changes: the code under `jcast/` is unchanged apart from its `__init__.py`. |
| **Question** | Does JCAST build, install and pass its tests on each Python version it should support, at the newest dependency versions and at the lowest that its bounds allow? Which bounds are true? |
| **Reproduce** | `python reports/packaging/check_install.py --uv <uv> --python <python 3.11> … --python <python 3.15> --resolution newest --jobs 5` (≈ 3 min) · `… --python <python 3.11> --resolution lowest` (≈ 2 min). The interpreters were uv's builds (`uv python install`), arm64. |
| **Evidence** (at `4a039a9`) | [packaging/check_install/check_newest.txt](packaging/check_install/check_newest.txt) · [check_lowest.txt](packaging/check_install/check_lowest.txt) |

## 1. Summary

- **JCAST works on Python 3.11 to 3.15** with the newest dependencies:
  pandas 3.1, numpy 2.4 (3.11) or 2.5, scipy 1.17 or 1.18, biopython 1.88,
  matplotlib 3.11. All 145 tests that need no network pass on each.
- **The bounds in `setup.py` were not all true.** On Python 3.11, the lowest
  versions they allowed either have no wheel (numpy 1.20, scipy 1.7,
  matplotlib 3.4) or fail to import (`requests` 2.0). The new lower bounds
  are the first releases with Python 3.11 wheels, or the lowest tested
  release (§3). With them, the tests pass.
- **One test depended on the pytest version**, through the order the tests
  run in (§4.1). Fixed in the test.
- **The packages:**
  - The wheel holds only the `jcast` package (62 kB). Before, `tests` was
    installed as a top-level package, and the chr15 genome, GTF and rMATS
    files (65 MB) into `<prefix>/tests/`.
  - The sdist (60 MB) holds the tests and all their data, so the tests run
    from it. The 0.3.5 sdist was 71 MB, without the UniProt and model test
    data.

## 2. What PR 10 changes

- **`pyproject.toml` holds the metadata** (`[project]`), and `setup.py` is
  removed. The version is set once, in `jcast/__init__.py`
  (`__version__ = '0.4.0'`), which `pyproject.toml` reads.
  `__version_info__` is removed.
- **Python ≥ 3.11**, with classifiers for 3.11 to 3.15. `setup.py` said
  ≥ 3.7, with classifiers for 3.6 to 3.9; the GUI already needed 3.10.
- **Only the `jcast` package is installed**
  (`[tool.setuptools.packages.find]`); no `data_files`.
- **The sdist** (`MANIFEST.in`) holds the tests and `tests/data/`, without
  the GTF caches (`*.cached`) and bytecode.
- **Dependencies** (§3). `requirements.txt` is removed: it repeated the
  dependencies, and listed scikit-learn, which `setup.py` did not.
  scikit-learn is now in the `test` extra, with pytest:
  `pip install -e '.[test]'`. It is used only by the Gaussian mixture in
  `model.py`, which the command line does not use, and its test, which is
  skipped without it.
- **`tox.ini`** runs the tests on Python 3.11 to 3.15 (it ran 3.7 to 3.9).
- **`jcast/__init__.py`** no longer builds `__all__` from the module files
  ([REVIEW.md §4](../REVIEW.md#low) LOW), which made
  `from jcast import *` import every module.
- **Classifiers:** "Science/Research" and "Bio-Informatics", instead of
  "Developers" and "Software Development :: Build Tools".

## 3. Dependency bounds

The lowest check installs, for each dependency, the lowest version within
its bounds (`uv pip install --resolution lowest-direct`). For the packages
that need compiling, it takes only wheels, so each is the oldest release
with a wheel for that Python.

| Dependency | Before | After | Why |
|---|---|---|---|
| numpy | ≥ 1.20 | ≥ 1.23.2 | first release with Python 3.11 wheels |
| scipy | ≥ 1.7 | ≥ 1.9.2 | the same |
| matplotlib | ≥ 3.4 | ≥ 3.6 | the same |
| requests | ≥ 2 | ≥ 2.20 | 2.0 fails to import on Python 3.11 (`collections.MutableMapping`), and 2.18 allows only urllib3 releases that fail the same way; 2.20 works |
| pandas | ≥ 1.5 | unchanged | 1.5.0 has Python 3.11 wheels and passes |
| biopython | ≥ 1.78 | unchanged | runs the whole test suite on Python 3.11 in the maintainer's conda environment; on Linux, 1.79 has Python 3.11 wheels; on macOS arm64, 1.83 is the first with them (the check's lowest) |
| gtfparse | ≥ 1.3, < 2 | unchanged | 1.3.0 is the only release in the bounds; pure Python, so it installs on 3.15. The 2.x and 3.x releases use polars (deferred to v0.5) |
| platformdirs, tqdm | ≥ 2.5, ≥ 4 | unchanged | 2.5.0 and 4.1.0 pass |
| scikit-learn (`test`) | ≥ 1.0 (`requirements.txt`) | ≥ 1.1.3 | first release with Python 3.11 wheels |
| pytest (`test`) | none | ≥ 8 | 8.0.0 passes |
| PySide6 (`gui`) | ≥ 6.6 | unchanged | not tested: no test imports the GUI |

The upper bounds are unchanged: the major version after the newest
tested (`pandas<4`, `numpy<3`, …).

## 4. Results

| Python (arm64) | Newest dependencies | Lowest dependencies |
|---|---|---|
| 3.11.18 | 145 passed (pandas 3.1.1, numpy 2.4.6, scipy 1.17.1) | 145 passed (pandas 1.5.0, numpy 1.23.2, scipy 1.9.2, matplotlib 3.6.0, biopython 1.83, requests 2.20.0, pytest 8.0.0) |
| 3.12.16 | 145 passed (numpy 2.5.4, scipy 1.18.1) | |
| 3.13.18 | 145 passed | |
| 3.14.10 | 145 passed | |
| 3.15.2 | 145 passed | |

- **The three network tests** were deselected (`-m 'not network'`), as
  in CI.
- **The command line:** in each environment, `jcast --help` and
  `python -m jcast --help` exit 0, and the installed version is 0.4.0.
- **The tests ran against the wheel:** from the sdist with its `jcast/`
  removed, so they imported the installed package.
- **The maintainer's conda environment** (Python 3.11, x86_64; pandas 2.3,
  numpy 1.26, biopython 1.78): the whole suite passes, network tests
  included (148 passed).

### 4.1 A test that depended on the pytest version

`test_unusable_cache_directory` checks that a cache directory that cannot
be created is logged. The warning goes to the `jcast` logger, which every
earlier `run_jcast` in the same pytest process leaves not propagating
(`get_logger` sets `propagate = False`). pytest 9.1 and later capture such
records anyway; pytest 8 does not, so the test failed whenever it ran after
a full run, as in the whole suite. The test now sets `propagate` itself.

### 4.2 Notes

- **Forking from a multi-threaded process.** On Python 3.12 and later, the
  two tests that run `-p` warn that the process is multi-threaded when it
  forks: earlier tests in the same process leave a thread running. A
  command-line run with `-p 4` on Python 3.15 gives no warning, because
  JCAST starts its workers before any thread exists
  ([parallel audit](2026-11-12-parallel-audit.md)). This matters for the
  GUI (PR 13), which runs in a Qt process with threads (REVIEW.md §10).
- **PyPI.** It has sdists only, up to 0.3.5 (71 MB each); 0.3.5a was not
  uploaded. A 0.4.0 wheel would make `pip install jcast` download 62 kB
  instead of the sdist.
- **The website** (`docs/install.html`) still says Python 3.7+ and `tox`.
  It is generated from sources not in the repository.

## 5. Decisions (2026-11-13)

1. Python ≥ 3.11 (3.10 reached its end of life in October 2026); test 3.11
   to 3.15, and check the dependencies.
2. Order: PR 10, then PR 11 (CI), then PR 12.
3. CI on Linux and macOS (the maintainer runs JCAST on macOS).
4. Version 0.4.0, to be released on GitHub after PRs 10 to 12.

## 6. CI (PR 11)

[.github/workflows/tests.yml](../.github/workflows/tests.yml) runs on every
push to a branch and on pull requests, without the three network tests:

| Job | Runs on | Steps |
|---|---|---|
| `test` (10 jobs) | `ubuntu-latest` and `macos-latest` (arm64), Python 3.11 to 3.15 | `pip install '.[test]'`; `pytest -m 'not network'` |
| `lowest` | `ubuntu-latest`, Python 3.11 | the lowest version of each dependency, as in §3, in a virtual environment; the tests |
| `package` | `ubuntu-latest`, Python 3.13 | `python -m build`; `twine check --strict`; the wheel holds only `jcast`; `jcast --help` and the tests, from the sdist against the installed wheel |

**Checked locally** on macOS arm64. Each job's `run` steps were taken from
the workflow file and run under `bash -e`, as GitHub runs them, on a fresh
checkout of `c31bf68` (PR 10 and its hash), in a fresh virtual environment
first on `PATH`:

| Job | Python | Result |
|---|---|---|
| `test` | 3.15.2 | 145 passed |
| `lowest` | 3.11.18 | 145 passed (pandas 1.5.0, numpy 1.23.2, scipy 1.9.2, matplotlib 3.6.0, biopython 1.83, requests 2.20.0) |
| `package` | 3.13.18 | both packages pass `twine check --strict`; 145 passed against the installed wheel |

- **The first local run of `lowest` installed into the wrong
  environment.** It used `uv pip install --system`, which found the
  machine's conda base environment through `CONDA_PREFIX`, not the job's
  Python. It added five packages there (removed again). The job now
  creates `.venv` and names it on every `uv` command, so it installs into
  the same place on any machine.
- **Not checked locally:**
  - Linux. There, the `lowest` job will take biopython 1.79, the oldest
    with Python 3.11 wheels.
  - The workflow on GitHub itself: it first runs when the branch is pushed.
    `actionlint` finds no problems in it.
