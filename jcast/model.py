# -*- coding: utf-8 -*-

""" Mixture model for selecting read count cutoff """

import os
import math
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import PowerTransformer
from pomegranate import GeneralMixtureModel, GammaDistribution, NormalDistribution, LogNormalDistribution





def general_mixture_model(sum_sjc_array, gamma_or_lognorm=0, n_iter=20):
    """Model skip junction histogram with a mixed distribution.

        Args:
            sum_sjc_array(numpy array): rMATS skipped-junction count results matrix.
            gamma_or_lognorm(int): Swtich to use Gamma or LogNormal Distribution to model low-end of skipped-junction histogram, 0 for Gamma - anything else for LogNormal.
            n_iter(int): Number of fit attempts to make, resulting distribution parameters will be averaged. 
        Returns:
            ln_sjc(nested list of floats): Natural log of read counts from rMATS, in nested list format for statistical classes.
            mean_mix_model(Pomegranate GeneralMixtureModel): Model object resulting from repeated fit attempts on provided sjc count data.
            min_count(int): Minimum number of reads to consider a sequence to be meaningfully expressed and not trace. 
            gamma_or_lognorm(int): Value of swtich to use Gamma or LogNormal Distribution for modeling low-end of skipped-junction histogram, 0 for Gamma - anything else for LogNormal.
    """
    passing_params = []
    ln_sjc = [[math.log(xi[0])] for xi in sum_sjc_array]
    for i in range(n_iter):
        if gamma_or_lognorm == 0:
            mix_model = GeneralMixtureModel.from_samples([GammaDistribution, NormalDistribution], n_components=2, X=ln_sjc)
        else:
            mix_model = GeneralMixtureModel.from_samples([LogNormalDistribution, NormalDistribution], n_components=2, X=ln_sjc)
        dist0 = mix_model.distributions[0]
        dist1 = mix_model.distributions[1]
        all_params = np.array(dist0.parameters+dist1.parameters)
        if not any(math.isnan(param) for param in all_params):
            passing_params.append(all_params)
    passing_params = np.array(passing_params)
    mean_params = [np.mean(passing_params[:,i]) for i in range(len(passing_params[0]))]

    if gamma_or_lognorm == 0:
        mean_dist0 = GammaDistribution(mean_params[0], mean_params[1])
        mean_dist0_ln_95_point = stats.gamma.ppf(0.95, mean_params[0], 0, 1/mean_params[1])
    else:
        mean_dist0 = LogNormalDistribution(mean_params[0], mean_params[1])
        mean_dist0_ln_95_point = stats.lognorm.ppf(0.95, mean_params[1], mean_params[0])
    mean_dist1 = NormalDistribution(mean_params[2], mean_params[3])
    mean_mix_model = GeneralMixtureModel([mean_dist0, mean_dist1])

    mean_dist0_95_point = math.ceil(math.exp(mean_dist0_ln_95_point))
    mean_dist1_ln_5_point = stats.norm.ppf(0.05, loc=mean_params[2], scale=mean_params[3])
    mean_dist1_5_point = math.ceil(math.exp(mean_dist1_ln_5_point))
    
    # If the 95th percentile of lower dist is beyond 5th percentile of higher dist, take 95th percentile as min_counts.
    # Else use 5th percentile as more stringent cutoff (this will probably never happen).
    if mean_dist0_95_point >= mean_dist1_5_point:
        min_count = mean_dist0_95_point
    else:
        min_count = mean_dist1_5_point

    return ln_sjc, mean_mix_model, min_count, gamma_or_lognorm

def plot_general_mixture_model(ln_sjc, 
                               mean_mix_model, 
                               min_count, 
                               gamma_or_lognorm, 
                               write_dir: str, 
                               filename: str):
    """Plot GeneralMixtureModel of skip junction read counts and fitted distributions.

        Args:
            ln_sjc(nested list of floats): Natural log of read counts from rMATS, in nested list format for statistical classes.
            mean_mix_model(Pomegranate GeneralMixtureModel): Model object resulting from repeated fit attempts on provided sjc count data.
            min_count(int): Minimum number of reads to consider a sequence to be meaningfully expressed and not trace. 
            gamma_or_lognorm(int): Swtich to use Gamma or LogNormal Distribution for modeling low-end of skipped-junction histogram, 0 for Gamma - anything else for LogNormal.
            write_dir(str): Path to directory for plot output.
            filename(str): Name of plot output file.
        Returns:
            None
    """
    mean_dist0 = mean_mix_model.distributions[0]
    mean_dist1 = mean_mix_model.distributions[1]
    mean_params = mean_dist0.parameters+mean_dist1.parameters

    x = np.linspace(0, 10, 2000)

    eval_dist0 = mean_dist0.probability(x)
    eval_dist1 = mean_dist1.probability(x)

    if gamma_or_lognorm == 0:
        # Pomegranate GammaDistribution uses shape and rate params, scipy uses shape and scale where scale = 1/rate.
        dist0_ln_95_point = stats.gamma.ppf(0.95, mean_params[0], 0, 1/mean_params[1])
    else:
        # Pomegranate LogNormalDistribution orders parameters [sigma, mu] while scipy.lognorm is [mu, sigma].
        dist0_ln_95_point = stats.lognorm.ppf(0.95, mean_params[1], mean_params[0])
    dist1_ln_5_point = stats.norm.ppf(0.05, loc=mean_params[2], scale=mean_params[3])
    dist0_95_point = math.ceil(math.exp(dist0_ln_95_point))
    dist1_5_point = math.ceil(math.exp(dist1_ln_5_point))

    dist0_ln_where = [True if xi < dist0_ln_95_point else False for xi in x]
    dist1_ln_where = [True if xi > dist1_ln_5_point else False for xi in x]

    nd_ln = np.asarray(ln_sjc)

    # Plot out the model figure and decision boundary
    fig, ax = plt.subplots(constrained_layout=True, figsize=(10, 4))
    ax.hist(nd_ln, bins=50, histtype='bar', density=True, ec='red', alpha=0.5)
    plt.style.use('seaborn-white')
    ax.set_xlabel('Natural logarithm of total skipped junction counts')
    ax.set_ylabel('Frequency')
    ax.set_title('Skipped junction count distributions')

    ax.plot(x, eval_dist0, color="C0")
    ax.plot(x, eval_dist1, color="C3")
    ax.fill_between(x, 0, eval_dist0, where=dist0_ln_where, color="C0", alpha=0.2, hatch="\\")
    ax.fill_between(x, 0, eval_dist1, where=dist1_ln_where, color="C3", alpha=0.2, hatch="/")
    ax.axvline(linewidth=2, ls='--', color='C0', x=dist0_ln_95_point)
    ax.axvline(linewidth=2, ls='--', color='C3', x=dist1_ln_5_point)
    ax.text(dist0_ln_95_point + math.log(1.25), 0.75, "min_counts = " + str(dist0_95_point))


    tcks, values = plt.xticks()
    new_tcks = np.int_(np.sort(np.asarray([math.exp(t) for t in tcks]).ravel()))

    # Second axis for untransformed reads
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(tcks[1:-1])
    ax2.set_xlabel('Untransformed total splice junction read counts')
    ax2.set_xticklabels(new_tcks[1:-1])

    plt.savefig(fname=os.path.join(write_dir, filename + '.png'))

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
