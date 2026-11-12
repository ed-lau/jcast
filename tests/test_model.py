# -*- coding: utf-8 -*-

""" Tests for read count model related functionalities """

import pytest

from jcast import model
from jcast.junctions import RmatsResults

from tests.conftest import DATA


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


@pytest.mark.parametrize('dataset', ['encode_human_heart', 'encode_mouse_heart', 'ipsc_human_d14'])
def test_gaussian_mixture(logger, tmp_path, dataset):
    """ test the (unused) scikit-learn Box-Cox Gaussian mixture path """

    min_jc_array = _junction_counts(logger, dataset)
    pt, gmm, min_count = model.gaussian_mixture(min_jc_array=min_jc_array)
    model.plot_model(min_jc_array=min_jc_array,
                     pt=pt, gmm=gmm, min_count=min_count,
                     write_dir=str(tmp_path), filename=dataset)

    assert gmm.means_[1] > gmm.means_[0]


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
