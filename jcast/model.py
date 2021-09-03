# -*- coding: utf-8 -*-

""" Mixture model for selecting read count cutoff """

import os
import numpy as np
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import PowerTransformer
import matplotlib.pyplot as plt
import scipy.stats as stats

# TODO: Include a gamma mixture model function, perhaps using pomegranate
def gaussian_mixture(sum_sjc_array):
    """
    Given a numpy array containing all splice junction sum SJC values, return a
    sci-kit learn Gaussian mixture model, along with the power transform object and
    the estimated minimum read count at the decision boundary

    :param sum_sjc_array: numpy array of junction counts
    :return: pt, gmm, min_count: scikit-learn power transform, gaussian mixture model, decision boundary
    """
    # Array of junction sum read counts

    pt = PowerTransformer(method='box-cox')
    pt.fit(sum_sjc_array)
    sum_sjc_array_transformed = pt.transform(sum_sjc_array)

    gmm = GaussianMixture(n_components=2,
                          covariance_type='diag',
                          random_state=1,
                          ).fit(sum_sjc_array_transformed)

    # Sort the classification components so the high-read distribution appears second
    order = gmm.means_.argsort(axis=0)[:, 0]
    gmm.means_ = gmm.means_[order]
    gmm.covariances_ = gmm.covariances_[order]
    gmm.weights_ = gmm.weights_[order]
    gmm.precisions_ = gmm.precisions_[order]
    gmm.precisions_cholesky_ = gmm.precisions_cholesky_[order]

    # Get the decision boundary (minimum count required to be predicted as second predicted class)
    # TODO: catch errors where model may fail and min_count remains 1
    min_count = 1

    # Loop only from 1 count to the mean of the second distribution
    for i in range(1, round(pt.inverse_transform([gmm.means_[1]])[0][0])):
        dist = gmm.predict(pt.transform(np.array([[i]])))
        if dist == [1]:
            min_count = i
            break

    return pt, gmm, min_count


def plot_model(sum_sjc_array,
               pt,
               gmm,
               min_count: int,
               write_dir: str,
               filename: str,
               ):
    """
    :param sum_sjc_array: numpy array of junction counts
    :param pt: scikit-learn power transform
    :param gmm: scikit-learn gaussian mixture model
    :param min_count: int decision boundary count
    :param write_dir: str output directory
    :param filename: str output name
    :return: True
    """

    sum_sjc_array_transformed = pt.transform(sum_sjc_array)

    weights = gmm.weights_
    means = gmm.means_
    covars = gmm.covariances_

    # Plot out the model figure and decision boundary
    fig, ax = plt.subplots(constrained_layout=True)
    ax.hist(sum_sjc_array_transformed, bins=50, histtype='bar', density=True, ec='red', alpha=0.5)
    plt.style.use('seaborn-white')
    ax.set_xlabel('Box-Cox transformed total skipped junction counts')
    ax.set_ylabel('Frequency')
    ax.set_title('Skipped junction count distributions')
    tcks, values = plt.xticks()
    new_tcks = np.int_(np.sort(pt.inverse_transform(np.array([[t] for t in tcks])).ravel()))

    # Second axis for untransformed reads
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(tcks[1:-1])
    ax2.set_xlabel('Untransformed total splice junction read counts')
    ax2.set_xticklabels(new_tcks[1:-1])
    ax2.axvline(linewidth=4, ls='--', color='blue', x=pt.transform(np.array([[min_count]])))

    f_axis = sum_sjc_array_transformed.copy().ravel()
    f_axis.sort()
    ax.plot(f_axis, weights[0] * stats.norm.pdf(f_axis, means[0], np.sqrt(covars[0])).ravel(), c='red')
    ax.plot(f_axis, weights[1] * stats.norm.pdf(f_axis, means[1], np.sqrt(covars[1])).ravel(), c='blue')

    plt.savefig(fname=os.path.join(write_dir, filename + '.png'))

    return True
