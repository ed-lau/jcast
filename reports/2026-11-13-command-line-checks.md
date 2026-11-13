# Command-line checks: the rMATS files read, and `-m` (PR 12)

| | |
|---|---|
| **Date** | 2026-11-13 |
| **Status** | Decisions taken 2026-11-13 (§3). **Implemented in PR 12** (`d84f49d`, §4). No output file changes for any run that worked before (§5). |
| **Question** | `-s` selects the splice types to translate, but JCAST read all five rMATS files whatever `-s` said ([REVIEW.md](../REVIEW.md) MED-8). A folder without one of them failed. Which files should be read, and what should `-m` fit on? |
| **Reproduce** | `python reports/cli_checks/verify_12.py --before <checkout> --after <checkout>` (≈ 3 min) |
| **Evidence** (at `d84f49d`) | [cli_checks/verify_12/verify_12.txt](cli_checks/verify_12/verify_12.txt) |

## 1. Summary

- **Before PR 12, a folder had to hold all five rMATS files.**
  - The command line checked only the files of the `-s` splice types, and
    with an `assert`, which `python -O` skips.
  - `RmatsResults` then read all five. A folder holding only
    `SE.MATS.JC.txt` passed the check, and the run failed with
    `FileNotFoundError`.
- **The `-m` model fitted on every file read, so on all five, whatever
  `-s` said.** Fitting on the selected splice types alone changes the
  cutoff a lot on chr15 (§2): 26 → 6 with `-s SE`, and no fit with
  `-s A5SS`.
- **Decided (§3):**
  - Only the selected files must exist.
  - With `-m`, the other files are read too when present, so the model
    still fits on every junction in the folder.
  - Every run that worked before gives the same output.
- **Also in PR 12:**
  - The missing-file check is a usage error.
  - `--g_or_ln` is removed. It was read by nothing.
  - So are the unused `get_rest_nuc` and the scikit-learn Gaussian
    mixture.

## 2. What `-m` fits on

The cutoff of the Gamma + Normal model (`model.general_mixture_model`) on
the bundled chr15 rMATS files, fitted on each splice type's junctions and
on all of them:

| Fitted on | Junctions | Cutoff (reads) |
|---|---|---|
| all five files | 1,105 | 26 |
| MXE | 90 | 57 |
| SE | 720 | 6 |
| RI | 88 | 13 |
| A5SS | 71 | no fit (falls back to `-r`) |
| A3SS | 136 | 4 |
| SE + MXE | 810 | 12 |

The model needs many junctions, and its cutoff reflects the sample's
sequencing depth, which all splice types share. Fitting on a subset would
make the cutoff of `-m -s SE` differ from that of `-m` for the same SE
junctions.

## 3. Decisions (2026-11-13)

1. **Require only the files of the `-s` splice types. With `-m`, also read
   the others that are present,** so that the model fits on every junction
   in the folder. The log records every rMATS file read, and says which
   were read only for the model. Rejected: reading only the selected files,
   which changes the `-m` cutoff when `-s` selects fewer than five (§2).
2. **Remove the scikit-learn Gaussian mixture** (`model.gaussian_mixture`,
   `model.plot_model`, its commented-out call in `main.py`, its test), and
   scikit-learn from the `test` extra. The command line never called it.
3. From the [packaging decisions](2026-11-13-packaging-check.md#5-decisions-2026-11-13):
   remove `--g_or_ln` and `helpers.get_rest_nuc`; report a missing rMATS
   file with `parser.error`.

## 4. What PR 12 changes

- **`RmatsResults(logger, rmats_dir, splice_types, others_if_present)`.**
  - Reads the files of `splice_types` (by default all five; each must
    exist).
  - With `others_if_present`, also reads the other files that exist.
    `main` passes `-s` and `-m`.
  - Its tables are in `tables`, by splice type. `rmats_mxe`, `rmats_se`,
    … are `None` for a file not read.
- **The log:**
  - One line per rMATS file read, with its sha256.
  - A file read only for the model is marked
    `(not selected by -s; read for the -m model)`.
  - Before, the log listed only the selected files, though `-m` read all
    five.
- **The command line:**
  - A missing rMATS folder or selected file is a usage error (exit status
    2), naming the files: `rMATS file not found: …; select the splice
    types the folder has with -s`.
  - `--g_or_ln` is rejected as an unrecognized argument.
  - The `-s` and `-m` help say which files are read.
- **Removed:**
  - `helpers.get_rest_nuc`, which fetched sequence from the Ensembl REST
    API with `print` logging and bare `except:`, and was called by
    nothing;
  - the scikit-learn path (§3).
- **Tests:** `tests/test_command_line.py`, 8 tests, all failing on the code
  before (`af187e8`).

## 5. Verification

`verify_12.py` ran the chr15 test set (`-c --uniprot`) with the code before
(`af187e8`) and after (frozen in a worktree):

| Run | Output files | Cutoff | Log |
|---|---|---|---|
| bundled folder, default | 14, byte-identical | 1 / 1 | summary identical; `Arguments` without `g_or_ln` |
| bundled folder, `-m` | 11, byte-identical | 26 / 26 | the same |
| bundled folder, `-s SE -m` | 11, byte-identical | 26 / 26 | the same, plus one line for each of the four files read only for the model |
| bundled folder, `-s SE RI` | 14, byte-identical | 1 / 1 | summary identical |
| folder with only `SE.MATS.JC.txt`, `-s SE` | before: `FileNotFoundError` for `MXE.MATS.JC.txt`; after: identical to the bundled folder with `-s SE` | | |
| folder with only `SE.MATS.JC.txt`, `-s SE -m` | after only | 6 | the model on the SE junctions alone |

The last row is the one case where the folder's contents set the cutoff:
with `-m`, the cutoff comes from every rMATS file present, and here only
SE is present. Before PR 12 this run failed.

The whole test suite passes: 153 tests (the 148 before, less the 3 of the
scikit-learn path, plus the 8 new).
