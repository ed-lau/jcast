# -*- coding: utf-8 -*-

""" Mixture model for selecting read count cutoff """

import os
import math
import numpy as np
import scipy.special as special
import scipy.stats as stats


def general_mixture_model(min_jc_array,
                          max_iter=500,
                          tol=1e-6,
                          min_weight=0.05):
    """Model the junction read count histogram with a Gamma (trace) + Normal (expressed) mixture, fitted by EM.

        Args:
            min_jc_array(numpy array): junction read counts (Junction.min_jc) plus 1.
            max_iter(int): Maximum number of EM iterations.
            tol(float): Relative change in log-likelihood at which EM is considered converged.
            min_weight(float): Fit is treated as degenerate if either component weight falls below this.
        Returns:
            ln_sjc(numpy array): Natural log of read counts from rMATS.
            params(tuple): (gamma shape, gamma rate, normal mean, normal sd, gamma weight, normal weight),
                or None if the fit is degenerate.
            min_count(int): The read count cutoff: a junction with at least this many reads is at or above the
                boundary between the two components (meaningfully expressed, not trace). None if the fit is
                degenerate.
    """
    ln_sjc = np.log(np.asarray(min_jc_array, dtype=float).ravel())
    x = np.clip(ln_sjc, 1e-6, None)

    # Initialize the gamma on the lower half of the data and the normal on the upper half.
    lo, hi = x[x <= np.median(x)], x[x > np.median(x)]
    shape, rate = lo.mean() ** 2 / lo.var(), lo.mean() / lo.var()
    mu, sd = hi.mean(), hi.std()
    weights = np.array([0.5, 0.5])

    ll_old = -np.inf
    for _ in range(max_iter):
        # E-step
        log_p = np.vstack([np.log(weights[0]) + stats.gamma.logpdf(x, shape, scale=1 / rate),
                           np.log(weights[1]) + stats.norm.logpdf(x, mu, sd)])
        ll = special.logsumexp(log_p, axis=0)
        resp = np.exp(log_p - ll)

        # M-step
        weights = resp.mean(axis=1)
        if weights.min() < min_weight:
            return ln_sjc, None, None
        mu = np.average(x, weights=resp[1])
        sd = math.sqrt(np.average((x - mu) ** 2, weights=resp[1]))
        # Gamma MLE: solve log(k) - digamma(k) = log(mean) - mean(log x), Minka (2002) start + Newton steps
        mean = np.average(x, weights=resp[0])
        s = math.log(mean) - np.average(np.log(x), weights=resp[0])
        if not s > 0:
            return ln_sjc, None, None
        shape = (3 - s + math.sqrt((s - 3) ** 2 + 24 * s)) / (12 * s)
        for _ in range(5):
            shape -= (math.log(shape) - special.digamma(shape) - s) / (1 / shape - special.polygamma(1, shape))
        rate = shape / mean
        if not all(np.isfinite([shape, rate, mu, sd])) or sd <= 0:
            return ln_sjc, None, None

        if ll.sum() - ll_old < tol * abs(ll.sum()):
            break
        ll_old = ll.sum()

    params = (shape, rate, mu, sd, weights[0], weights[1])

    # Decision boundary b: where the weighted densities cross between the gamma mode and the normal mean. The model is
    # fitted to ln(count + 1), so a junction is at or above b when count + 1 >= e^b, i.e. count >= ceil(e^b) - 1.
    # 2026-11-12 this returned ceil(e^b), on the count + 1 scale, which was then compared to read counts with <=: the
    # cut was two reads above the boundary (reports/2026-11-12-thresholds-runs-canonical-audit.md §2.3)
    grid = np.linspace(max((shape - 1) / rate, 1e-3), mu, 2000)
    diff = (weights[0] * stats.gamma.pdf(grid, shape, scale=1 / rate)
            - weights[1] * stats.norm.pdf(grid, mu, sd))
    crossings = np.flatnonzero(np.diff(np.sign(diff)))
    if len(crossings) == 0:
        return ln_sjc, None, None
    min_count = math.ceil(np.exp(grid[crossings[0] + 1])) - 1

    return ln_sjc, params, min_count


def plot_general_mixture_model(ln_sjc,
                               params,
                               min_count,
                               write_dir: str,
                               filename: str):
    """Plot the Gamma + Normal mixture fitted to skip junction read counts.

        Args:
            ln_sjc(numpy array): Natural log of read counts from rMATS.
            params(tuple): Output of general_mixture_model.
            min_count(int): The read count cutoff from general_mixture_model.
            write_dir(str): Path to directory for plot output.
            filename(str): Name of plot output file.
        Returns:
            None
    """
    import matplotlib.pyplot as plt

    shape, rate, mu, sd, w0, w1 = params
    x = np.linspace(0, 10, 2000)
    eval_dist0 = stats.gamma.pdf(x, shape, scale=1 / rate) * w0
    eval_dist1 = stats.norm.pdf(x, mu, sd) * w1

    fig, ax = plt.subplots(constrained_layout=True, figsize=(10, 4))
    ax.hist(ln_sjc, bins=50, histtype='bar', density=True, ec='red', alpha=0.5)
    ax.set_xlabel('Natural logarithm of read count + 1')
    ax.set_ylabel('Frequency')
    ax.set_title('Junction read counts (smaller of summed inclusion and skipping junction reads)')
    ax.set_ylim(0, 1)

    ax.plot(x, eval_dist0, color="C0")
    ax.plot(x, eval_dist1, color="C3")
    ax.plot(x, eval_dist0 + eval_dist1, color="black")
    ax.fill_between(x, 0, eval_dist0, color="C0", alpha=0.2)
    ax.fill_between(x, 0, eval_dist1, color="C3", alpha=0.2)
    ax.axvline(linewidth=2, ls='--', color='black', x=np.log(min_count + 1))
    ax.text(np.log(min_count + 1) + math.log(1.25), 0.75, "cutoff = " + str(min_count) + " reads")

    # Second axis for untransformed reads
    tcks = ax.get_xticks()
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(tcks[1:-1])
    ax2.set_xlabel('Read count')
    ax2.set_xticklabels((np.exp(tcks[1:-1]) - 1).round().astype(int))

    fig.savefig(fname=os.path.join(write_dir, filename + '.pdf'))
    plt.close(fig)
