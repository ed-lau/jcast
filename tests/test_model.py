# -*- coding: utf-8 -*-

""" Tests for read count model related functionalities """

import csv
import math
import shutil

import numpy as np
import pandas as pd
import pytest

from jcast import main, model
from jcast.junctions import Junction, RmatsResults, read_count

from tests.conftest import DATA, HUMAN_GENOME, HUMAN_GTF, RMATS

DATASETS = sorted(p.name for p in (DATA / 'model').iterdir() if p.is_dir())


def _junction_counts(logger, dataset):
    rmats_results = RmatsResults(logger=logger, rmats_dir=DATA / 'model' / dataset)
    rmats_results.get_junction_count_array()
    return rmats_results.min_jc_array


@pytest.mark.parametrize('dataset, expected_min_count', [
    ('encode_human_heart', 14),
    ('encode_mouse_heart', 6),
    ('ipsc_human_d14', 15),
])
def test_general_mixture_model(logger, tmp_path, dataset, expected_min_count):
    """
    test the Gamma + Normal read count cutoff against values recorded when replacing pomegranate. 2026-11-12 the
    cutoff is in reads; it was on the count + 1 scale, one more (15, 7, 16); see
    reports/2026-11-12-thresholds-runs-canonical-audit.md §2.3
    """

    ln_sjc, params, min_count = model.general_mixture_model(min_jc_array=_junction_counts(logger, dataset))

    assert min_count == expected_min_count
    shape, rate, mu, sd, w0, w1 = params
    assert mu > (shape - 1) / rate     # normal (expressed) component sits above the gamma (trace) mode
    assert w0 + w1 == pytest.approx(1)

    model.plot_general_mixture_model(ln_sjc, params, min_count, write_dir=str(tmp_path), filename=dataset)
    assert (tmp_path / f'{dataset}.pdf').exists()


def test_general_mixture_model_degenerate(logger):
    """ stomach has no stable two-component fit; the model should report that rather than crash """

    _, params, min_count = model.general_mixture_model(min_jc_array=_junction_counts(logger, 'encode_human_stomach'))

    assert params is None
    assert min_count is None


def test_cutoff_is_the_models_boundary_in_reads(logger):
    """ the cutoff c is the smallest read count whose ln(count + 1) lies at or above the components' crossing """
    import math
    import numpy as np
    from scipy import stats

    counts = _junction_counts(logger, 'encode_human_heart')
    _, (shape, rate, mu, sd, w0, w1), cutoff = model.general_mixture_model(min_jc_array=counts)

    def expressed(reads):
        x = math.log(reads + 1)
        return w1 * stats.norm.pdf(x, mu, sd) >= w0 * stats.gamma.pdf(x, shape, scale=1 / rate)
    assert expressed(cutoff) and not expressed(cutoff - 1)


# ---- PR 19: what the GUI needs (reports/2026-11-13-gui-cli-scope.md §8) ----

@pytest.mark.parametrize('folder', [RMATS] + [DATA / 'model' / d for d in DATASETS])
def test_read_counts_as_from_each_junction(logger, folder):
    """ the read counts from the tables' columns are those of a Junction made for each row, as before PR 19 """
    results = RmatsResults(logger=logger, rmats_dir=folder)
    results.get_junction_count_array()
    tot = pd.concat(results.rmats_list, ignore_index=True)
    junctions = [Junction(**tot.iloc[i].to_dict()) for i in range(len(tot))]

    assert list(results.counts.min_jc) == [j.min_jc for j in junctions]
    assert list(results.counts.junction_type) == [j.junction_type for j in junctions]
    assert list(results.counts.fdr) == [j.fdr for j in junctions]
    expected = np.array([[j.min_jc + 1] for j in junctions])
    assert results.min_jc_array.dtype == expected.dtype and np.array_equal(results.min_jc_array, expected)


@pytest.mark.parametrize('counts, expected', [
    ((3, 5, 4, 1), 6),                          # min(ijc 3 + 4, sjc 5 + 1)
    (('1,2,3', '4,5', '0,0,1', '6'), 7),        # replicates, comma-separated
    ((float('nan'), 5, 4, 1), 0),               # an inclusion count that does not parse: both samples' count 0
    (('5.0', 5, 4, 1), 0),                      # a float, from a column with a missing value, does not parse either
])
def test_read_count(counts, expected):
    assert read_count(*counts) == expected


def test_preview_is_the_runs_model(tmp_path):
    """ the -m model as a run fits it: on every rMATS file in the folder (cutoff 26 on chr15), here also with one
    file only (6, as the run's with only SE.MATS.JC.txt present; reports/2026-11-13-command-line-checks.md §5) """
    counts, ln_sjc, params, cutoff, no_fit = model.preview(RMATS)
    assert (len(counts), cutoff, no_fit) == (1105, 26, None) and len(ln_sjc) == 1105 and params is not None

    shutil.copy(RMATS / 'SE.MATS.JC.txt', tmp_path)
    assert model.preview(tmp_path)[3] == 6

    with pytest.raises(FileNotFoundError):
        model.preview(tmp_path / 'none')


@pytest.mark.slow
def test_filter_counts_are_the_runs(tmp_path, cache_dir):
    """ what the preview says the filters keep and skip, per splice type, is what a run's junctions.tsv records """
    out = tmp_path / 'out'
    main.main([str(a) for a in [RMATS, HUMAN_GTF, HUMAN_GENOME, '-s', 'RI', 'A5SS', '-r', 5, '-q', 0, 0.05, '-o', out,
                                '--cache_dir', cache_dir]])
    with open(out / 'junctions.tsv') as fh:
        rows = list(csv.DictReader(fh, delimiter='\t'))

    table = model.filter_counts(model.preview(RMATS)[0], 5, (0, 0.05))
    for splice_type in ['RI', 'A5SS']:
        fates = [r['fate'] for r in rows if r['junction_type'] == splice_type]
        low, fdr = fates.count('skipped_low_count'), fates.count('skipped_fdr')
        assert table.loc[splice_type].tolist() == [len(fates), low, fdr, len(fates) - low - fdr]
        assert low and fdr
    assert list(table.index) == ['MXE', 'SE', 'RI', 'A5SS', 'A3SS']


def test_draw_model_on_a_callers_axes():
    """ the GUI draws the model on an Axes of its own """
    from matplotlib.figure import Figure

    _, ln_sjc, params, cutoff, _ = model.preview(RMATS)
    fig = Figure()
    ax = fig.subplots()
    model.draw_model(ax, ln_sjc, params, cutoff)

    assert len(fig.axes) == 2       # and the twin axis in reads
    assert len(ax.lines) == 4 and ax.lines[3].get_xdata()[0] == pytest.approx(math.log(cutoff + 1))
    assert ax.texts[0].get_text() == 'cutoff = 26 reads'


# ---- PR 20: when the -m model cannot be fitted, the run falls back to -r ----

EMPTY = DATA / 'model' / 'encode_mouse_fetalheart'     # the five rMATS files, each its header only


def test_fit_without_junctions():
    """ 2026-11-13 general_mixture_model raised ZeroDivisionError, and the run failed """
    ln_sjc, params, cutoff, no_fit = model.fit(np.array([]))

    assert (len(ln_sjc), params, cutoff) == (0, None, None)
    assert no_fit == 'failed: the rMATS files read hold no junctions'


def test_fit_that_fails(monkeypatch, logger):
    """ any failure of the fit is no fit, with its reason """
    def fails(min_jc_array):
        raise FloatingPointError('overflow')
    monkeypatch.setattr(model, 'general_mixture_model', fails)

    assert model.fit(_junction_counts(logger, 'encode_human_heart'))[1:] == \
           (None, None, 'failed (FloatingPointError: overflow)')


def test_fit_degenerate_and_fitted(logger):
    assert model.fit(_junction_counts(logger, 'encode_human_stomach'))[1:] == (None, None, 'was degenerate')
    assert model.fit(_junction_counts(logger, 'encode_human_heart'))[2:] == (14, None)


def test_preview_without_junctions():
    assert model.preview(EMPTY)[3:] == (None, 'failed: the rMATS files read hold no junctions')


@pytest.mark.parametrize('extra, read', [([], 1), (['-r', 3], 3)])
def test_a_run_without_junctions_falls_back_to_r(tmp_path, cache_dir, capsys, extra, read):
    """ -m on rMATS files without junctions: the run warns, uses -r (1 by default) and completes. 2026-11-13 it failed
    with ZeroDivisionError """
    out = tmp_path / 'out'
    main.main([str(a) for a in [EMPTY, HUMAN_GTF, HUMAN_GENOME, '-m', '-o', out, '--cache_dir', cache_dir] + extra])

    warning = f'Mixture model fit failed: the rMATS files read hold no junctions; falling back to -r {read}.'
    assert warning in capsys.readouterr().err
    log = (out / 'logfile.log').read_text()
    assert f'~WARNING~{warning}' in log and f'read count of at least {read} (-r)' in log
    assert (out / 'junctions.tsv').read_text().count('\n') == 1 and not (out / 'model.pdf').exists()


def test_a_run_whose_fit_fails_falls_back_to_r(tmp_path, cache_dir, capsys, monkeypatch):
    def fails(min_jc_array):
        raise FloatingPointError('overflow')
    monkeypatch.setattr(model, 'general_mixture_model', fails)
    out = tmp_path / 'out'
    main.main([str(a) for a in [RMATS, HUMAN_GTF, HUMAN_GENOME, '-m', '-s', 'RI', '-r', 100000, '-o', out,
                                '--cache_dir', cache_dir]])

    assert 'Mixture model fit failed (FloatingPointError: overflow); falling back to -r 100000.' in \
           capsys.readouterr().err
    assert 'read count of at least 100000 (-r)' in (out / 'logfile.log').read_text()


@pytest.mark.parametrize('folder, extra, expected', [
    (RMATS, ['-m'], 'at least 26 (-m)'),
    (RMATS, ['-r', 100000], 'at least 100000 (-r)'),
    (DATA / 'model' / 'encode_human_stomach', ['-m', '-r', 100000], 'at least 100000 (-r)'),     # a degenerate fit
])
def test_the_cutoff_names_its_option(tmp_path, cache_dir, folder, extra, expected):
    """ the log names the option the cutoff came from: -m when fitted; -r without -m, or when -m has no fit.
    2026-11-13 it said -m whenever -m was given """
    out = tmp_path / 'out'
    main.main([str(a) for a in [folder, HUMAN_GTF, HUMAN_GENOME, '-s', 'RI', '-o', out, '--cache_dir', cache_dir]
               + extra])
    assert f'Translating junctions with a read count of {expected}:' in (out / 'logfile.log').read_text()
