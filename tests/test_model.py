# -*- coding: utf-8 -*-

""" Tests for read count model related functionalities """

import pytest

from jcast import model
from jcast.junctions import RmatsResults

from tests.conftest import DATA


def _junction_counts(logger, dataset):
    rmats_results = RmatsResults(logger=logger, rmats_dir=DATA / 'model' / dataset)
    rmats_results.get_junction_count_array()
    return rmats_results.sum_sjc_array


@pytest.mark.parametrize('dataset, expected_min_count', [
    ('encode_human_heart', 15),
    ('encode_mouse_heart', 7),
    ('ipsc_human_d14', 16),
])
def test_general_mixture_model(logger, tmp_path, dataset, expected_min_count):
    """ test the Gamma + Normal read count cutoff against values recorded when replacing pomegranate """

    ln_sjc, params, min_count = model.general_mixture_model(sum_sjc_array=_junction_counts(logger, dataset))

    assert min_count == expected_min_count
    shape, rate, mu, sd, w0, w1 = params
    assert mu > (shape - 1) / rate     # normal (expressed) component sits above the gamma (trace) mode
    assert w0 + w1 == pytest.approx(1)

    model.plot_general_mixture_model(ln_sjc, params, min_count, write_dir=str(tmp_path), filename=dataset)
    assert (tmp_path / f'{dataset}.pdf').exists()


def test_general_mixture_model_degenerate(logger):
    """ stomach has no stable two-component fit; the model should report that rather than crash """

    _, params, min_count = model.general_mixture_model(sum_sjc_array=_junction_counts(logger, 'encode_human_stomach'))

    assert params is None
    assert min_count is None


@pytest.mark.parametrize('dataset', ['encode_human_heart', 'encode_mouse_heart', 'ipsc_human_d14'])
def test_gaussian_mixture(logger, tmp_path, dataset):
    """ test the (unused) scikit-learn Box-Cox Gaussian mixture path """

    sum_sjc_array = _junction_counts(logger, dataset)
    pt, gmm, min_count = model.gaussian_mixture(sum_sjc_array=sum_sjc_array)
    model.plot_model(sum_sjc_array=sum_sjc_array,
                     pt=pt, gmm=gmm, min_count=min_count,
                     write_dir=str(tmp_path), filename=dataset)

    assert gmm.means_[1] > gmm.means_[0]
