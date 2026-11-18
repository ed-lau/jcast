# -*- coding: utf-8 -*-

"""Row 23 audit (reports/2026-11-15-read-count-model-audit.md §8): how well the -m model fits junction read counts,
and how stable its cutoff is.

The -m model (jcast/model.py) is a Gamma ("trace") plus a Normal ("expressed") on ln(count + 1), fitted by EM to every
junction's read count (Junction.min_jc: the smaller of the summed inclusion and skipping junction counts). The cutoff is
where the weighted densities cross. The EM here (em) is general_mixture_model on the distinct counts, each weighted by
how many junctions have it; main checks that it gives the same fit.

Per rMATS folder, into --out:
    fit.tsv         one row per model, on the count scale, where the data are: log-likelihood, BIC and the cutoff the
                    model implies (the lowest count whose posterior for the higher component is at least 0.5). A model
                    on ln(count + 1) gives a count c the mass of its density over [ln(c + 1/2), ln(c + 3/2)]. Every
                    model is truncated at 1: rMATS reports no event with fewer reads in either form. The models:
                        -m as run         general_mixture_model's fit, not refitted (its cutoff is JCAST's)
                        Gamma + Normal    the same family, fitted by maximum likelihood on the counts
                        Normal + Normal   on ln(count + 1), as in the 2019 paper
                        NB + NB           two negative binomials
                        NB + PLN          a negative binomial (trace) and a Poisson-lognormal (expressed)
                        NB x K            K negative binomials (K = 1, 3, 4)
                    The fits are by L-BFGS-B from a few starts (starts_for); a likelihood here is a lower bound on
                    the family's best.
    rootogram.tsv   observed and expected junctions per count bin, per model
    stability.tsv   the -m cutoff by bootstrap, by the EM's start (the lower part is x <= its q quantile; the code
                    uses the median) and stopping tolerance, run to convergence, and with each count spread evenly
                    over its interval (the EM then sees no lattice of ln 2, ln 3, ...)
    fit_summary.txt the check that em gives general_mixture_model's fit, and the bootstrap's range
    subsets.tsv     the -m model fitted to each splice type alone, and to the events rMATS calls annotated or novel
                    (fromGTF.novelJunction / novelSpliceSite.<type>.txt, when present)
    depth.tsv       the -m model fitted to the first k replicates' counts, summed (k = 1 .. all, sample 1's replicates
                    first): how the cutoff grows with the reads summed
    form.tsv        per splice type, which form's count is the smaller, as counted and per junction position (rMATS's
                    IncFormLen and SkipFormLen)
    fit.png         junctions per read count, observed and expected (left); observed / expected per bin (right)

With --bundled <dir>, every rMATS folder in it, into --out/bundled.tsv: the -m cutoff as run and converged, its
bootstrap range, and each two-component model's BIC and cutoff.

Usage:
    python reports/model_audit/model_fit.py --rmats <dir> [--bootstrap 100] --out <dir>
    python reports/model_audit/model_fit.py --bundled tests/data/model --out <dir>
"""

import argparse
import math
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.optimize as optimize
import scipy.special as special
import scipy.stats as stats

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from jcast.model import general_mixture_model  # noqa: E402

SPLICE_TYPES = ['SE', 'MXE', 'A5SS', 'A3SS', 'RI']
EDGES = [1, 2, 3, 4, 5, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024, 2048, 4096, 8192]
TAIL = EDGES[-1]                      # the last bin is >= TAIL

# the reference palette's light steps (dataviz skill, references/palette.md), as in cutoff_recovery.py
INK, INK2, GRID, SURFACE = '#0b0b0b', '#52514e', '#e4e3df', '#fcfcfb'
SERIES = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4']


# --- reading -------------------------------------------------------------------------------------------------------

def replicates(col: pd.Series) -> np.ndarray:
    """ a column of comma-separated replicate counts as a matrix, one column per replicate (not a number: 0) """
    return col.astype(str).str.split(',', expand=True).apply(pd.to_numeric, errors='coerce').fillna(0).to_numpy(int)


def read_events(rmats_dir) -> pd.DataFrame:
    """ per event: splice type, rMATS's call (annotated or novel, None without the fromGTF files), the summed inclusion
    and skipping counts, IncFormLen and SkipFormLen, and the replicate counts (I and S: sample 1's, then sample 2's) """
    out = []
    for st in SPLICE_TYPES:
        f = Path(rmats_dir) / f'{st}.MATS.JC.txt'
        if not f.exists():
            continue
        m = pd.read_csv(f, sep='\t')
        ijc = np.hstack([replicates(m['IJC_SAMPLE_1']), replicates(m['IJC_SAMPLE_2'])])
        sjc = np.hstack([replicates(m['SJC_SAMPLE_1']), replicates(m['SJC_SAMPLE_2'])])
        files = [Path(rmats_dir) / f'fromGTF.{k}.{st}.txt' for k in ['novelJunction', 'novelSpliceSite']]
        if any(p.exists() for p in files):
            novel = set()
            for p in files:
                if p.exists():
                    novel |= set(pd.read_csv(p, sep='\t', usecols=['ID'])['ID'].astype(str))
            call = np.where(m['ID'].astype(str).isin(novel), 'novel', 'annotated')
        else:
            call = None
        out.append(pd.DataFrame({'junction_type': st, 'call': call, 'ijc': ijc.sum(1), 'sjc': sjc.sum(1),
                                 'inc_len': m['IncFormLen'], 'skip_len': m['SkipFormLen'],
                                 'I': list(ijc), 'S': list(sjc)}))
    ev = pd.concat(out, ignore_index=True)
    ev['min_jc'] = np.minimum(ev['ijc'], ev['sjc'])
    return ev


def distinct(counts) -> tuple:
    """ the distinct counts and how many junctions have each """
    c, w = np.unique(np.asarray(counts, dtype=int), return_counts=True)
    return c, w.astype(float)


# --- the -m model, weighted -----------------------------------------------------------------------------------------

def _quantile(x, w, q) -> float:
    """ np.quantile (linear) of the data that x with multiplicities w stand for; x sorted """
    n = w.sum()
    pos = q * (n - 1)
    cum = np.cumsum(w)
    lo, hi = x[np.searchsorted(cum, math.floor(pos) + 1)], x[np.searchsorted(cum, math.ceil(pos) + 1)]
    return lo + (hi - lo) * (pos - math.floor(pos))


def crossing(params) -> int:
    """ general_mixture_model's cutoff for these parameters """
    shape, rate, mu, sd, w0, w1 = params
    grid = np.linspace(max((shape - 1) / rate, 1e-3), mu, 2000)
    diff = w0 * stats.gamma.pdf(grid, shape, scale=1 / rate) - w1 * stats.norm.pdf(grid, mu, sd)
    cross = np.flatnonzero(np.diff(np.sign(diff)))
    return None if len(cross) == 0 else math.ceil(np.exp(grid[cross[0] + 1])) - 1


def em(counts, w, q=0.5, max_iter=500, tol=1e-6, min_weight=0.05) -> dict:
    """ general_mixture_model on distinct counts with multiplicities w, starting from a split at the q quantile """
    x = np.clip(np.log(np.asarray(counts, dtype=float) + 1), 1e-6, None)
    split = _quantile(x, w, q)
    lo, hi = x <= split, x > split
    m_lo = np.average(x[lo], weights=w[lo])
    v_lo = np.average((x[lo] - m_lo) ** 2, weights=w[lo])
    shape, rate = m_lo ** 2 / v_lo, m_lo / v_lo
    mu = np.average(x[hi], weights=w[hi])
    sd = math.sqrt(np.average((x[hi] - mu) ** 2, weights=w[hi]))
    weights = np.array([0.5, 0.5])
    ll_old, it, total = -np.inf, 0, w.sum()
    for it in range(1, max_iter + 1):
        log_p = np.vstack([np.log(weights[0]) + stats.gamma.logpdf(x, shape, scale=1 / rate),
                           np.log(weights[1]) + stats.norm.logpdf(x, mu, sd)])
        ll = special.logsumexp(log_p, axis=0)
        resp = np.exp(log_p - ll) * w
        weights = resp.sum(axis=1) / total
        if weights.min() < min_weight:
            return {'params': None, 'cutoff': None, 'iterations': it}
        mu = np.average(x, weights=resp[1])
        sd = math.sqrt(np.average((x - mu) ** 2, weights=resp[1]))
        mean = np.average(x, weights=resp[0])
        s = math.log(mean) - np.average(np.log(x), weights=resp[0])
        if not s > 0:
            return {'params': None, 'cutoff': None, 'iterations': it}
        shape = (3 - s + math.sqrt((s - 3) ** 2 + 24 * s)) / (12 * s)
        for _ in range(5):
            shape -= (math.log(shape) - special.digamma(shape) - s) / (1 / shape - special.polygamma(1, shape))
        rate = shape / mean
        total_ll = (ll * w).sum()
        if total_ll - ll_old < tol * abs(total_ll):
            break
        ll_old = total_ll
    params = (shape, rate, mu, sd, weights[0], weights[1])
    return {'params': params, 'cutoff': crossing(params), 'iterations': it}


# --- models on the count scale -------------------------------------------------------------------------------------

def _logdiff(a, b):
    """ log(exp(a) - exp(b)) for a >= b """
    with np.errstate(divide='ignore', invalid='ignore'):
        return a + np.log(-np.expm1(np.minimum(b - a, 0)))


def log_interval(dist, lo, hi):
    """ log P(lo < X <= hi) for a frozen scipy distribution, from the nearer tail """
    upper = lo > dist.median()
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(upper, _logdiff(dist.logsf(lo), dist.logsf(hi)), _logdiff(dist.logcdf(hi), dist.logcdf(lo)))


def pln_logpmf(c, mu, sigma, split=300):
    """ the Poisson-lognormal's log probability of each count c: below split, Poisson(c | e^z) integrated over z on a
    grid; from split up, the lognormal's mass over [c - 1/2, c + 1/2] (the Poisson's spread is then far below sigma) """
    c = np.asarray(c, dtype=float)
    out = np.empty(c.shape)
    small = c < split
    if small.any():
        z = np.linspace(mu - 9 * sigma, mu + 9 * sigma, 3001)
        lw = stats.norm.logpdf(z, mu, sigma) + np.log(z[1] - z[0])
        out[small] = special.logsumexp(stats.poisson.logpmf(c[small][:, None], np.exp(z)[None, :]) + lw[None, :],
                                       axis=1)
    if (~small).any():
        out[~small] = log_interval(stats.norm(mu, sigma), np.log(c[~small] - 0.5), np.log(c[~small] + 0.5))
    return out


class Component:
    """ one mixture component on the counts: its log probability of each count, and its location (to order them) """
    n_params = 2

    def __init__(self, kind):
        self.kind = kind

    def logpmf(self, c, a, b):
        if self.kind == 'gamma_ln':
            d = stats.gamma(math.exp(a), scale=math.exp(-b))
            return log_interval(d, np.log(c + 0.5), np.log(c + 1.5))
        if self.kind == 'norm_ln':
            d = stats.norm(a, math.exp(b))
            return log_interval(d, np.log(c + 0.5), np.log(c + 1.5))
        if self.kind == 'nb':
            # size from 1e-6 (near the log-series the zero-truncated NB tends to; smaller, the truncation divides
            # rounding errors) to 1e6 (near Poisson); p = size / (size + mean)
            b = min(max(b, -13.8), 13.8)
            return stats.nbinom.logpmf(c, math.exp(b), special.expit(b - a))
        if self.kind == 'pln':
            return pln_logpmf(c, a, math.exp(b))
        raise ValueError(self.kind)

    def log_zero(self, a, b):
        """ log of the mass below a count of 1 """
        if self.kind in ('gamma_ln', 'norm_ln'):
            d = stats.gamma(math.exp(a), scale=math.exp(-b)) if self.kind == 'gamma_ln' else stats.norm(a, math.exp(b))
            return float(d.logcdf(math.log(1.5)))
        return float(self.logpmf(np.array([0.0]), a, b)[0])

    def location(self, a, b):
        """ the component's centre on the ln-count scale """
        return {'gamma_ln': math.exp(a - b), 'norm_ln': a, 'nb': a, 'pln': a}[self.kind]


class Mixture:
    """ a mixture of Components on the counts, truncated at 1; theta: each component's two parameters (unconstrained),
    then K - 1 weight logits """

    def __init__(self, name, kinds):
        self.name, self.parts = name, [Component(k) for k in kinds]

    @property
    def n_params(self):
        return 3 * len(self.parts) - 1

    def unpack(self, theta):
        k = len(self.parts)
        pairs = [(theta[2 * i], theta[2 * i + 1]) for i in range(k)]
        logits = np.concatenate([[0.0], theta[2 * k:]])
        return pairs, logits - special.logsumexp(logits)

    def component_logpmf(self, c, theta):
        """ (K, len(c)): log weight + log probability, untruncated; and the log of the mixture's mass below 1 """
        pairs, lw = self.unpack(theta)
        lp = np.vstack([lw[i] + p.logpmf(c, *pairs[i]) for i, p in enumerate(self.parts)])
        l0 = special.logsumexp([lw[i] + p.log_zero(*pairs[i]) for i, p in enumerate(self.parts)])
        return lp, l0

    def loglik(self, theta, c, w):
        with np.errstate(all='ignore'):
            lp, l0 = self.component_logpmf(c, theta)
            ll = special.logsumexp(lp, axis=0) - np.log(-np.expm1(l0))
        total = float((ll * w).sum())
        return total if np.isfinite(total) else -1e300

    def fit(self, c, w, starts):
        best = None
        for theta0 in starts:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                r = optimize.minimize(lambda t: -self.loglik(t, c, w), np.asarray(theta0, float), method='L-BFGS-B',
                                      options={'maxiter': 3000, 'maxfun': 30000})
            if best is None or r.fun < best.fun:
                best = r
        return best.x, -best.fun

    def probabilities(self, theta, upto):
        """ each count's probability under the truncated mixture, 1 .. upto - 1, and the posterior of the component
        with the highest location """
        c = np.arange(1, upto, dtype=float)
        with np.errstate(all='ignore'):
            lp, l0 = self.component_logpmf(c, theta)
            total = special.logsumexp(lp, axis=0)
        pairs, _ = self.unpack(theta)
        high = int(np.argmax([p.location(*pairs[i]) for i, p in enumerate(self.parts)]))
        return np.exp(total - np.log(-np.expm1(l0))), np.exp(lp[high] - total)


def em_theta(params) -> np.ndarray:
    """ general_mixture_model's parameters as a Gamma + Normal Mixture's theta """
    shape, rate, mu, sd, w0, w1 = params
    return np.array([math.log(shape), math.log(rate), mu, math.log(sd), math.log(w1 / w0)])


def nb_moments(c, w):
    m = np.average(c, weights=w)
    v = np.average((c - m) ** 2, weights=w)
    size = m ** 2 / (v - m) if v > m * 1.01 else 100.0
    return math.log(m), math.log(size)


def starts_for(kinds, c, w, em_params) -> list:
    """ starting points: a split of the counts at each of a few quantiles, each part's moments for its component """
    starts = []
    k = len(kinds)
    if kinds == ['gamma_ln', 'norm_ln'] and em_params is not None:
        starts.append(em_theta(em_params))
    cum = np.cumsum(w) / w.sum()
    splits = [[0.5], [0.3], [0.7]] if k == 2 else [list(np.linspace(0, 1, k + 1)[1:-1])]
    if k == 1:
        splits = [[]]
    for qs in splits:
        edges = [0] + [int(np.searchsorted(cum, q)) + 1 for q in qs] + [len(c)]
        for i in range(1, k):                        # at least two distinct counts in each part
            edges[i] = max(edges[i], edges[i - 1] + 2)
        theta, ok = [], True
        for kind, (a, b) in zip(kinds, zip(edges[:-1], edges[1:])):
            cc, ww = c[a:b], w[a:b]
            if len(cc) < 2:
                ok = False
                break
            lx = np.log(cc + 1.0)
            m = np.average(lx, weights=ww)
            s = math.sqrt(max(np.average((lx - m) ** 2, weights=ww), 1e-3))
            if kind == 'gamma_ln':
                theta += [math.log(m ** 2 / s ** 2), math.log(m / s ** 2)]
            elif kind == 'norm_ln':
                theta += [m, math.log(s)]
            elif kind == 'nb':
                theta += list(nb_moments(cc, ww))
            else:
                theta += [float(np.average(np.log(cc), weights=ww)), math.log(s)]
        if ok:
            starts.append(theta + [0.0] * (k - 1))
    return starts


MODELS = [('Gamma + Normal', ['gamma_ln', 'norm_ln']), ('Normal + Normal', ['norm_ln', 'norm_ln']),
          ('NB + NB', ['nb', 'nb']), ('NB + PLN', ['nb', 'pln']), ('NB x 1', ['nb']),
          ('NB x 3', ['nb', 'nb', 'nb']), ('NB x 4', ['nb', 'nb', 'nb', 'nb'])]


def first_at_least(posterior, level=0.5):
    """ the lowest count whose posterior is at least level, from the first count up """
    idx = np.flatnonzero(posterior >= level)
    return int(idx[0]) + 1 if len(idx) else None


def fit_models(c, w, em_params, models=MODELS) -> tuple:
    """ fit.tsv's rows, and each model's probability of each count (1 .. TAIL - 1) """
    n = w.sum()
    rows, probs = [], {}
    fitted = [('-m as run', Mixture('-m as run', ['gamma_ln', 'norm_ln']), em_theta(em_params))] \
        if em_params is not None else []
    for name, kinds in models:
        mix = Mixture(name, kinds)
        fitted.append((name, mix, mix.fit(c, w, starts_for(kinds, c, w, em_params))[0]))
    for name, mix, theta in fitted:
        ll = mix.loglik(theta, c, w)
        p, post = mix.probabilities(theta, TAIL)
        probs[name] = p
        pairs, lw = mix.unpack(theta)
        rows.append({'model': name, 'parameters': mix.n_params, 'log-likelihood': round(ll, 1),
                     'BIC': round(-2 * ll + mix.n_params * math.log(n), 1),
                     'cutoff (posterior of the higher component >= 0.5)':
                         first_at_least(post) if len(mix.parts) == 2 else None,
                     'JCAST cutoff (the densities cross)': crossing(em_params) if name == '-m as run' else None,
                     'components': '; '.join(f'{p.kind}({a:.3f}, {b:.3f})' for p, (a, b) in zip(mix.parts, pairs)),
                     'weights': ', '.join(f'{x:.3f}' for x in np.exp(lw))})
    t = pd.DataFrame(rows)
    t['delta BIC'] = (t['BIC'] - t['BIC'].min()).round(1)
    return t, probs


def rootogram(c, w, probs) -> pd.DataFrame:
    """ observed and expected junctions per count bin (EDGES; the last bin is >= TAIL) """
    n = w.sum()
    labels = [f'{a}' if b - a == 1 else f'{a}-{b - 1}' for a, b in zip(EDGES[:-1], EDGES[1:])] + [f'>={TAIL}']
    idx = np.searchsorted(EDGES, c, side='right') - 1
    observed = np.bincount(idx, weights=w, minlength=len(EDGES))
    t = pd.DataFrame({'bin': labels, 'from': EDGES, 'width': list(np.diff(EDGES)) + [np.nan],
                      'observed': observed.astype(int)})
    counts = np.arange(1, TAIL)
    for name, p in probs.items():
        e = np.bincount(np.searchsorted(EDGES, counts, side='right') - 1, weights=p, minlength=len(EDGES))
        e[-1] = max(1 - p.sum(), 0)
        t[name] = (n * e).round(1)
    return t


# --- stability, subsets, depth, form --------------------------------------------------------------------------------

def stability(c, w, n_boot, seed=23) -> pd.DataFrame:
    rows = []
    for q in [0.35, 0.5, 0.65]:
        for tol in [1e-4, 1e-5, 1e-6, 1e-7, 1e-8]:
            r = em(c, w, q=q, tol=tol)
            rows.append({'test': 'start and tolerance', 'start quantile': q, 'tolerance': tol, 'max_iter': 500,
                         'iterations': r['iterations'], 'cutoff': r['cutoff']})
        r = em(c, w, q=q, tol=1e-12, max_iter=200000)
        rows.append({'test': 'converged', 'start quantile': q, 'tolerance': 1e-12, 'max_iter': 200000,
                     'iterations': r['iterations'], 'cutoff': r['cutoff'],
                     'gamma weight': None if r['params'] is None else round(r['params'][4], 4)})
    rng = np.random.default_rng(seed)
    p = w / w.sum()
    for b in range(n_boot):
        wb = rng.multinomial(int(w.sum()), p).astype(float)
        keep = wb > 0
        r = em(c[keep], wb[keep])
        rows.append({'test': 'bootstrap', 'replicate': b, 'start quantile': 0.5, 'tolerance': 1e-6, 'max_iter': 500,
                     'iterations': r['iterations'], 'cutoff': r['cutoff']})
    # the counts spread evenly over the intervals a count model gives them, ln(c + 1/2) to ln(c + 3/2): the EM then
    # sees a continuous variable, not the lattice ln 2, ln 3, ...
    every = np.repeat(c, w.astype(int)).astype(float)
    for d in range(3):
        r = em(every - 0.5 + rng.uniform(size=len(every)), np.ones(len(every)))
        rows.append({'test': 'spread over its interval', 'replicate': d, 'start quantile': 0.5, 'tolerance': 1e-6,
                     'max_iter': 500, 'iterations': r['iterations'], 'cutoff': r['cutoff'],
                     'gamma weight': None if r['params'] is None else round(r['params'][4], 4)})
    return pd.DataFrame(rows)


def kept_share(counts, cutoff):
    return None if cutoff is None else round(float((np.asarray(counts) >= cutoff).mean()), 4)


def subsets(ev) -> pd.DataFrame:
    groups = [('all', ev)] + [(st, ev[ev['junction_type'] == st]) for st in SPLICE_TYPES]
    if ev['call'].notna().any():
        groups += [(f'{k} events', ev[ev['call'] == k]) for k in ['annotated', 'novel']]
    rows = []
    for name, g in groups:
        if len(g) < 50:
            continue
        r = em(*distinct(g['min_jc']))
        rc = em(*distinct(g['min_jc']), tol=1e-12, max_iter=200000)
        rows.append({'subset': name, 'junctions': len(g), 'median count': float(g['min_jc'].median()),
                     'cutoff': r['cutoff'], 'junctions kept': kept_share(g['min_jc'], r['cutoff']),
                     'gamma weight': None if r['params'] is None else round(r['params'][4], 4),
                     'cutoff, converged': rc['cutoff']})
    return pd.DataFrame(rows)


def depth(ev) -> pd.DataFrame:
    i_mat, s_mat = np.vstack(ev['I'].to_numpy()), np.vstack(ev['S'].to_numpy())
    rows = []
    for k in range(1, i_mat.shape[1] + 1):
        count = np.minimum(i_mat[:, :k].sum(1), s_mat[:, :k].sum(1))
        count = count[count >= 1]
        r = em(*distinct(count))
        rows.append({'replicates summed': k, 'junctions with a read in each form': len(count),
                     'junction reads (both forms, all events)': int(i_mat[:, :k].sum() + s_mat[:, :k].sum()),
                     'median count': float(np.median(count)), 'cutoff': r['cutoff'],
                     'junctions kept': kept_share(count, r['cutoff'])})
    return pd.DataFrame(rows)


def form(ev) -> pd.DataFrame:
    rows = []
    for st in ['all'] + SPLICE_TYPES:
        g = ev if st == 'all' else ev[ev['junction_type'] == st]
        if not len(g):
            continue
        raw = g['ijc'] < g['sjc']
        norm = g['ijc'] / g['inc_len'] < g['sjc'] / g['skip_len']
        rows.append({'splice_type': st, 'events': len(g),
                     'median IncFormLen / SkipFormLen': float((g['inc_len'] / g['skip_len']).median()),
                     'smaller: inclusion, as counted': round(float(raw.mean()), 4),
                     'smaller: inclusion, per junction position': round(float(norm.mean()), 4),
                     'the two disagree': round(float((raw != norm).mean()), 4)})
    return pd.DataFrame(rows)


# --- the chart ------------------------------------------------------------------------------------------------------

def _style(ax):
    ax.set_facecolor(SURFACE)
    for s in ['top', 'right']:
        ax.spines[s].set_visible(False)
    for s in ['left', 'bottom']:
        ax.spines[s].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def plot_fit(rg, path, title, shown=('-m as run', 'Gamma + Normal', 'NB + PLN')):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    t = rg.iloc[:-1]                                  # the open last bin has no width
    x = np.sqrt(t['from'] * (t['from'] + t['width'] - 1))   # the bin's geometric centre
    fig, (left, right) = plt.subplots(1, 2, figsize=(10.5, 4.4), facecolor=SURFACE)
    _style(left)
    left.plot(x, t['observed'] / t['width'], color=INK, linestyle='none', marker='o', markersize=5,
              label='observed')
    for name, color in zip(shown, SERIES):
        left.plot(x, t[name] / t['width'], color=color, linewidth=2, label=name)
    left.set_xscale('log')
    left.set_yscale('log')
    left.set_xlabel('read count (min of summed inclusion and skipping counts)', color=INK, fontsize=9)
    left.set_ylabel('junctions per read count', color=INK, fontsize=9)
    left.set_title('Observed and expected', color=INK, fontsize=10, loc='left')
    left.legend(frameon=False, fontsize=8.5, loc='lower left', labelcolor=INK)
    _style(right)
    for name, color in zip(shown, SERIES):
        right.plot(x, t['observed'] / t[name], color=color, linewidth=2, marker='o', markersize=4, label=name)
    right.axhline(1, color=INK2, linewidth=0.8, linestyle=(0, (3, 3)))
    right.set_xscale('log')
    right.set_yscale('log')
    right.set_xlabel('read count', color=INK, fontsize=9)
    right.set_ylabel('observed / expected', color=INK, fontsize=9)
    right.set_title('Observed / expected, per bin', color=INK, fontsize=10, loc='left')
    for ax in (left, right):
        ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f'{v:g}'))
    right.set_ylim(0.5, 2)
    right.set_yticks([0.5, 0.67, 0.8, 1, 1.25, 1.5, 2])
    right.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f'{v:g}'))
    right.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    fig.suptitle(title, color=INK, fontsize=11, x=0.01, ha='left')
    fig.tight_layout()
    fig.savefig(path, dpi=150, facecolor=SURFACE)
    plt.close(fig)


# --- main -----------------------------------------------------------------------------------------------------------

def one_folder(rmats_dir, out, n_boot, name):
    ev = read_events(rmats_dir)
    c, w = distinct(ev['min_jc'])
    ref_ln, ref_params, ref_cutoff = general_mixture_model(np.asarray(ev['min_jc']) + 1)
    r = em(c, w)
    same = ref_params is not None and np.allclose(r['params'], ref_params, rtol=1e-6) and r['cutoff'] == ref_cutoff
    lines = [f'{rmats_dir}: {len(ev)} events, {len(c)} distinct counts',
             f'general_mixture_model: cutoff {ref_cutoff}, parameters {ref_params}',
             f'em (weighted, distinct counts): cutoff {r["cutoff"]}, parameters {r["params"]}; same: {same}',
             f'junctions with a count of 1: {int(w[c == 1].sum())} ({w[c == 1].sum() / w.sum():.1%}); '
             f'kept at the cutoff: {kept_share(ev["min_jc"], ref_cutoff)}']
    fit, probs = fit_models(c, w, ref_params)
    rg = rootogram(c, w, probs)
    st = stability(c, w, n_boot)
    sub, dp, fm = subsets(ev), depth(ev), form(ev)
    for t, f in [(fit, 'fit'), (rg, 'rootogram'), (st, 'stability'), (sub, 'subsets'), (dp, 'depth'), (fm, 'form')]:
        t.to_csv(out / f'{f}.tsv', sep='\t', index=False)
    boot = st.loc[st['test'] == 'bootstrap', 'cutoff']
    lines.append(f'bootstrap ({len(boot)}): cutoffs {boot.min()} to {boot.max()}, '
                 f'2.5-97.5%: {boot.quantile(0.025)} to {boot.quantile(0.975)}')
    (out / 'fit_summary.txt').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))
    with pd.option_context('display.width', 250, 'display.max_columns', 40, 'display.max_colwidth', 80):
        for t in (fit, rg, st[st['test'] != 'bootstrap'], sub, dp, fm):
            print(t.to_string(index=False))
    plot_fit(rg, out / 'fit.png', f'Junction read counts and the models fitted to them: {name}')


def bundled(root, out, n_boot):
    rows = []
    for d in sorted(p for p in Path(root).iterdir() if (p / 'SE.MATS.JC.txt').exists()):
        ev = read_events(d)
        if ev.empty:                                 # encode_mouse_fetalheart: the files hold headers only
            rows.append({'folder': d.name, 'junctions': 0})
            continue
        c, w = distinct(ev['min_jc'])
        r, rc = em(c, w), em(c, w, tol=1e-12, max_iter=200000)
        rng = np.random.default_rng(23)
        boot = []
        for _ in range(n_boot):
            wb = rng.multinomial(int(w.sum()), w / w.sum()).astype(float)
            boot.append(em(c[wb > 0], wb[wb > 0])['cutoff'])
        boot = pd.Series(boot, dtype=float)
        fit, _ = fit_models(c, w, r['params'], [m for m in MODELS if len(m[1]) == 2])
        row = {'folder': d.name, 'junctions': len(ev), 'count 0': int((ev['min_jc'] == 0).sum()),
               'median count': float(ev['min_jc'].median()), '-m cutoff': r['cutoff'],
               '-m cutoff, converged': rc['cutoff'], 'bootstrap 2.5%': boot.quantile(0.025),
               'bootstrap 97.5%': boot.quantile(0.975), 'bootstrap: no fit': int(boot.isna().sum())}
        for _, f in fit.iterrows():
            row[f'{f["model"]}: delta BIC'] = f['delta BIC']
            row[f'{f["model"]}: cutoff'] = f['cutoff (posterior of the higher component >= 0.5)']
        rows.append(row)
        print(row, flush=True)
    t = pd.DataFrame(rows)
    t.to_csv(out / 'bundled.tsv', sep='\t', index=False)
    with pd.option_context('display.width', 250, 'display.max_columns', 40):
        print(t.T.to_string())


def main():
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--rmats', help='an rMATS folder')
    p.add_argument('--bundled', help='a folder of rMATS folders')
    p.add_argument('--name', default=None, help='the data set, for the chart title')
    p.add_argument('--bootstrap', type=int, default=100)
    p.add_argument('--out', required=True)
    a = p.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    if a.bundled:
        bundled(a.bundled, out, a.bootstrap)
    else:
        one_folder(a.rmats, out, a.bootstrap, a.name or Path(a.rmats).name)


if __name__ == '__main__':
    main()
