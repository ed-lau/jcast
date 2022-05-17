# -*- coding: utf-8 -*-

""" Mixture model for selecting read count cutoff """

import os
import math
import numpy as np
import scipy.stats as stats
import scipy.signal as signal
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import PowerTransformer
from pomegranate import GeneralMixtureModel, GammaDistribution, NormalDistribution, LogNormalDistribution


def general_mixture_model(sum_sjc_array, n_iter=50):
    """Model skip junction histogram with a mixed distribution.

        Args:
            sum_sjc_array(numpy array): rMATS skipped-junction count results matrix.
            n_iter(int): Number of fit attempts to make, resulting distribution parameters will be averaged. 
        Returns:
            ln_sjc(nested list of floats): Natural log of read counts from rMATS, in nested list format for statistical classes.
            mean_mix_model(Pomegranate GeneralMixtureModel): Model object resulting from repeated fit attempts on provided sjc count data.
            min_count(int): Minimum number of reads to consider a sequence to be meaningfully expressed and not trace. 
    """
    passing_params = []
    ln_sjc = [[math.log(xi[0])] for xi in sum_sjc_array]
    for i in range(n_iter):
        mix_model = GeneralMixtureModel.from_samples([GammaDistribution, NormalDistribution], n_components=2, X=ln_sjc)
        dist0 = mix_model.distributions[0]
        dist1 = mix_model.distributions[1]
        all_params = dist0.parameters+dist1.parameters
        if not any(math.isnan(param) for param in all_params):
            out = np.array(all_params + np.exp(mix_model.weights).tolist())
            passing_params.append(out)
    
    passing_params = np.array(passing_params)
    n_params = len(passing_params[0])
    n_samples = len(passing_params[:,0])
    # Takes gaussian KDE of estimated values for each parameter, picks highest peak in KDE as param.
    best_params = []
    x = np.linspace(-5, 20, 500000)
    for i in range(len(passing_params[0])):
        kde = stats.gaussian_kde(passing_params[:,i])
        y = kde.evaluate(x)
        peaks = signal.find_peaks(y)[0]
        peak_heights = [y[j] for j in peaks]
        best_params.append(x[peaks[np.argmax(peak_heights)]])

    best_dist0 = GammaDistribution(best_params[0], best_params[1])
    # Pomegranate GammaDistribution uses shape and rate params, scipy uses shape and scale where scale = 1/rate.
    best_dist1 = NormalDistribution(best_params[2], best_params[3])
    best_mix_model = GeneralMixtureModel([best_dist0, best_dist1], weights=best_params[4:])
    
    x = np.linspace(0, 10, 2000)
    eval_gamma = best_dist0.probability(x) * best_params[4]
    eval_norm = best_dist1.probability(x) * best_params[5]
    
    # Always an intersection near zero, index 1 is intersection of interest.
    intersection_idx = np.argwhere(np.diff(np.sign(eval_gamma - eval_norm))).flatten()
    min_count = math.ceil(np.exp(x[intersection_idx][1]))

    return ln_sjc, best_mix_model, min_count


def plot_general_mixture_model(ln_sjc, 
                               best_mix_model, 
                               min_count, 
                               write_dir: str, 
                               filename: str):
    """Plot Pomegranate GeneralMixtureModel of skip junction read counts and fitted distributions.

        Args:
            ln_sjc(nested list of floats): Natural log of read counts from rMATS, in nested list format for statistical classes.
            best_mix_model(Pomegranate GeneralMixtureModel): Model object resulting from repeated fit attempts on provided sjc count data.
            min_count(int): Minimum number of reads to consider a sequence to be meaningfully expressed and not trace. 
            gamma_or_lognorm(int): Swtich to use Gamma or LogNormal Distribution for modeling low-end of skipped-junction histogram, 0 for Gamma - anything else for LogNormal.
            write_dir(str): Path to directory for plot output.
            filename(str): Name of plot output file.
        Returns:
            None
    """
    best_dist0 = best_mix_model.distributions[0]
    best_dist1 = best_mix_model.distributions[1]
    best_params = best_dist0.parameters+best_dist1.parameters
    best_weights = np.exp(best_mix_model.weights).tolist()

    x = np.linspace(0, 10, 2000)

    eval_dist0 = best_dist0.probability(x) * best_weights[0]
    eval_dist1 = best_dist1.probability(x) * best_weights[1]
    eval_mixture_dist = best_mix_model.probability(x)

    nd_ln = np.asarray(ln_sjc)

    # Plot out the model figure and decision boundary
    fig, ax = plt.subplots(constrained_layout=True, figsize=(10, 4))
    ax.hist(nd_ln, bins=50, histtype='bar', density=True, ec='red', alpha=0.5)
    plt.style.use('seaborn-white')
    ax.set_xlabel('Natural logarithm of total skipped junction counts')
    ax.set_ylabel('Frequency')
    ax.set_title('Skipped junction count distributions')
    ax.set_ylim(0,1)

    ax.plot(x, eval_dist0, color="C0")
    ax.plot(x, eval_dist1, color="C3")
    ax.plot(x, eval_mixture_dist, color="black")
    ax.fill_between(x, 0, eval_dist0, color="C0", alpha=0.2)
    ax.fill_between(x, 0, eval_dist1, color="C3", alpha=0.2)
    ax.axvline(linewidth=2, ls='--', color='black', x=np.log(min_count))
    ax.text(np.log(min_count) + math.log(1.25), 0.75, "min_counts = " + str(min_count))

    tcks, values = plt.xticks()
    new_tcks = np.int_(np.sort(np.asarray([math.exp(t) for t in tcks]).ravel()))

    # Second axis for untransformed reads
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(tcks[1:-1])
    ax2.set_xlabel('Untransformed total splice junction read counts')
    ax2.set_xticklabels(new_tcks[1:-1])

    plt.savefig(fname=os.path.join(write_dir, filename + '.pdf'))


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
