# -*- coding: utf-8 -*-

""" Mixture model for selecting read count cutoff """

import logging
import os
import math
import numpy as np
import pandas as pd
import scipy.special as special
import scipy.stats as stats

from jcast import fates
from jcast.junctions import RmatsResults
from jcast.outcomes import SPLICE_TYPES

# the read count bins of the frame chart, [lo, next lo), as in reports/2026-11-15-read-count-model-audit.md §4.4
FRAME_BINS = [1, 2, 3, 5, 8, 16, 32, 64, 128, 256, 512, 1024, 4096]
FRAME_MIN_EVENTS = 20           # a bin's share is drawn for a group of at least so many events

# the charts' colours: ink, secondary ink and grid, and the groups' (dataviz palette, slots 1 and 2, validated on a
# light surface), each with a marker so that colour is not the only cue
INK, INK2, GRID = '#0b0b0b', '#52514e', '#e4e3df'
GROUPS = {'annotated': ('#2a78d6', 'o'), 'novel': ('#eb6834', 's'), 'all': ('#2a78d6', 'o')}


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


def fit(min_jc_array) -> tuple:
    """
    The -m model, as the run and the GUI's preview fit it: general_mixture_model, where a fit that fails is no fit, as
    a degenerate one is, so that the run warns and falls back to -r. 2026-11-13 a failure stopped the run: rMATS files
    holding no junctions raised ZeroDivisionError (reports/2026-11-13-gui-cli-scope.md §11)

    :param min_jc_array: junction read counts (Junction.min_jc) plus 1
    :return: ln(count + 1), the parameters and the cutoff, as general_mixture_model returns them (the last two None if
        there is no fit); and why there is no fit, completing 'Mixture model fit …' (None if there is one)
    """
    ln_sjc = np.log(np.asarray(min_jc_array, dtype=float).ravel())
    if ln_sjc.size == 0:
        return ln_sjc, None, None, 'failed: the rMATS files read hold no junctions'
    try:
        ln_sjc, params, min_count = general_mixture_model(min_jc_array=min_jc_array)
    except Exception as e:
        return ln_sjc, None, None, f'failed ({type(e).__name__}: {e})'
    return ln_sjc, params, min_count, None if min_count is not None else 'was degenerate'


def preview(rmats_dir) -> tuple:
    """
    The -m model of an rMATS folder, fitted as a run with -m fits it: on every rMATS file in the folder, also those -s
    leaves out. It needs no GTF or genome, so the GUI shows it before the run (reports/2026-11-13-gui-cli-scope.md §4).

    :param rmats_dir: the rMATS folder
    :return: every junction's splice type, read count and FDR (RmatsResults.junction_counts), with its rMATS ID,
        whether its alternative region keeps the frame and whether rMATS calls it novel (PR 27; events); then, as fit
        returns them, ln(count + 1), the parameters, the cutoff (None if there is no fit: the run uses -r), and why
        there is no fit
    """
    results = RmatsResults(logging.getLogger('jcast'), rmats_dir, splice_types=(), others_if_present=True)
    if not results.tables:
        raise FileNotFoundError(f'no rMATS file (<splice type>.MATS.JC.txt) in {rmats_dir}')
    results.get_junction_count_array()
    return (events(results), *fit(results.min_jc_array))


def alternative_length(table, splice_type: str) -> pd.Series:
    """
    The length of each event's alternative region, in nt: the skipped exon (SE), the retained intron (RI), or the
    difference between the two forms' exons (MXE, A5SS, A3SS). With a multiple of 3, both forms keep the reading frame
    of the rest of the gene; a random length does in one event of three (reports/2026-11-15-read-count-model-audit.md
    §4.4)

    :param table: an RmatsResults table: rMATS's coordinates, the starts 0-based
    :param splice_type: its splice type
    """
    alt1 = table.alt1_ee - table.alt1_es
    if splice_type == 'SE':
        return alt1
    if splice_type == 'RI':
        return table.down_es - table.anc_ee
    return alt1 - (table.alt2_ee - table.alt2_es)


def novel_events(rmats_dir, splice_type: str):
    """
    The IDs of the events that rMATS calls novel: those in fromGTF.novelJunction.<type>.txt, and with --novelSS in
    fromGTF.novelSpliceSite.<type>.txt. rMATS before 4.1.1 calls some annotated events novel
    (reports/2026-11-15-read-count-model-audit.md §4.1)

    :return: a set of IDs; None if neither file is in the folder
    """
    files = [os.path.join(rmats_dir, f'fromGTF.{kind}.{splice_type}.txt') for kind in ('novelJunction',
                                                                                        'novelSpliceSite')]
    files = [f for f in files if os.path.isfile(f)]
    if not files:
        return None
    return set().union(*(pd.read_csv(f, sep='\t', usecols=['ID'])['ID'] for f in files))


def events(results) -> pd.DataFrame:
    """
    The junctions read (results.counts, after get_junction_count_array), with each event's rMATS ID (junction_id),
    whether its alternative region is a multiple of 3 nt (frame_kept), and whether rMATS calls it novel (novel; NA for
    a splice type whose novel-event files are not in the folder). For the Model page's charts (PR 27)
    """
    ids, frame, novel = [], [], []
    for splice_type, table in results.tables.items():
        ids.append(table.id)
        frame.append(alternative_length(table, splice_type) % 3 == 0)
        called = novel_events(results.dir, splice_type)
        novel.append(pd.Series(pd.NA if called is None else table.id.isin(called), index=table.index,
                               dtype='boolean'))
    if not ids:
        return results.counts.assign(junction_id=pd.Series(dtype=int), frame_kept=pd.Series(dtype=bool),
                                     novel=pd.Series(dtype='boolean'))
    return results.counts.assign(junction_id=pd.concat(ids, ignore_index=True).to_numpy(),
                                 frame_kept=pd.concat(frame, ignore_index=True).to_numpy(),
                                 novel=pd.concat(novel, ignore_index=True).array)


def filter_counts(counts, cutoff: int, qvalue) -> pd.DataFrame:
    """
    How many junctions of each splice type a run with this read count cutoff and -q would translate, and how many its
    filters would skip: fates.filtered, which the run applies to each junction

    :param counts: junctions, as preview returns them
    :param cutoff: the read count cutoff (-r, or the -m model's)
    :param qvalue: (q_lo, q_hi) (-q)
    :return: one row per splice type, in the order of SPLICE_TYPES; columns junctions, skipped_low_count, skipped_fdr,
        kept
    """
    skipped = pd.Series([fates.filtered(n, fdr, cutoff, qvalue) for n, fdr in zip(counts.min_jc, counts.fdr)],
                        index=counts.index, dtype=object)
    by_type = counts.junction_type
    table = pd.DataFrame({'junctions': by_type.value_counts(),
                          'skipped_low_count': (skipped == fates.SKIPPED_LOW_COUNT).groupby(by_type).sum(),
                          'skipped_fdr': (skipped == fates.SKIPPED_P_VALUE).groupby(by_type).sum(),
                          })
    table['kept'] = table.junctions - table.skipped_low_count - table.skipped_fdr
    return table.reindex([t for t in SPLICE_TYPES if t in table.index])


def _groups(counts) -> dict:
    """ The annotated and the novel events, or all of them if rMATS's call is missing for one: name -> mask """
    novel = counts.novel
    if novel.isna().any() or not len(counts):
        return {'all': np.ones(len(counts), dtype=bool)}
    novel = novel.to_numpy(dtype=bool)
    return {'annotated': ~novel, 'novel': novel}


def kept_at_cutoffs(counts, qvalue, splice_types=SPLICE_TYPES, cutoffs=None) -> pd.DataFrame:
    """
    How many junctions a run would translate at each read count cutoff: those of the splice types whose read count is
    at least the cutoff and whose FDR -q keeps, as fates.filtered decides; annotated and novel (PR 27)

    :param counts: junctions, as preview returns them
    :param qvalue: (q_lo, q_hi) (-q)
    :param splice_types: the splice types selected (-s)
    :param cutoffs: the cutoffs; by default about 200, evenly on a log scale from 1 to the read count that 0.1 % of
        the junctions reach (at least 10), past which the counts barely change
    :return: index the cutoffs; a column per group: annotated and novel, or all if rMATS's call is missing
    """
    counts = counts[counts.junction_type.isin(splice_types)]
    q_lo, q_hi = qvalue
    fdr = counts.fdr.to_numpy(dtype=float)
    with np.errstate(invalid='ignore'):
        q_kept = np.where(np.isnan(fdr), (q_lo, q_hi) == fates.ALL_FDRS, (q_lo <= fdr) & (fdr <= q_hi))
    read = counts.min_jc.to_numpy()
    if cutoffs is None:
        top = max(int(np.quantile(read[q_kept], 0.999)) if q_kept.any() else 1, 10)
        cutoffs = np.unique(np.round(np.geomspace(1, top, 200)).astype(int))
    cutoffs = np.asarray(cutoffs)
    kept = {}
    for name, mask in _groups(counts).items():
        n = np.sort(read[mask & q_kept])
        kept[name] = len(n) - np.searchsorted(n, cutoffs, side='left')
    return pd.DataFrame(kept, index=pd.Index(cutoffs, name='cutoff'))


def frame_by_count(counts, splice_types=SPLICE_TYPES) -> pd.DataFrame:
    """
    By read count (FRAME_BINS), the events of the splice types, and the share whose alternative region is a multiple
    of 3 nt: annotated and novel, or all (PR 27; the read-count audit's §4.4 chart, there for the junctions translated)

    :return: index the bins' lower bounds; columns '<group> events' and '<group> frame kept' (a fraction) per group
    """
    counts = counts[counts.junction_type.isin(splice_types)]
    bins = pd.cut(counts.min_jc, FRAME_BINS + [np.inf], right=False, labels=FRAME_BINS)
    table = {}
    for name, mask in _groups(counts).items():
        by_bin = counts.frame_kept[mask].groupby(bins[mask], observed=False)
        table[f'{name} events'] = by_bin.size()
        table[f'{name} frame kept'] = by_bin.mean()
    return pd.DataFrame(table).rename_axis('bin')


def _style(ax) -> None:
    for side in ['top', 'right']:
        ax.spines[side].set_visible(False)
    for side in ['left', 'bottom']:
        ax.spines[side].set_color(INK2)
    ax.tick_params(colors=INK2, labelsize=8.5)
    ax.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.xaxis.set_major_formatter('{x:,.0f}')


def _label(name: str) -> str:
    return 'all events (rMATS\'s novel-event files are not in the folder)' if name == 'all' else f'{name} events'


def _cutoff_line(ax, cutoff) -> tuple:
    line = ax.axvline(max(cutoff, 1), color=INK, linewidth=1.5, linestyle=(0, (4, 3)))
    label = ax.text(0, 1.01, '', transform=ax.get_xaxis_transform(), color=INK, fontsize=8.5)
    return line, label


def move_chart_cutoff(line, label, cutoff, text: str = None) -> None:
    """ Put a chart's cutoff (draw_kept, draw_frame) at another read count; with text, its label's. The label is
    right of the line, or left of it in the right half of the chart """
    x = max(cutoff, 1)
    line.set_xdata([x, x])
    lo, hi = (math.log(max(v, 1e-9)) for v in label.axes.get_xlim())
    right_half = math.log(x) > (lo + hi) / 2
    label.set_x(x / 1.08 if right_half else x * 1.08)
    label.set_horizontalalignment('right' if right_half else 'left')
    label.set_text(f'cutoff = {cutoff} reads' + (f': {text}' if text else ''))


def draw_kept(ax, kept, cutoff) -> tuple:
    """
    The junctions a run would translate at each read count cutoff (kept_at_cutoffs), annotated and novel, and the
    cutoff (PR 27)

    :return: the cutoff's line and label, which move_chart_cutoff moves
    """
    _style(ax)
    for name in kept.columns:
        color, _ = GROUPS[name]
        ax.plot(kept.index, kept[name], color=color, linewidth=2, label=_label(name))
    ax.set_xscale('log')
    ax.xaxis.set_major_formatter('{x:,.0f}')
    ax.yaxis.set_major_formatter('{x:,.0f}')
    ax.set_ylim(bottom=0)
    ax.set_xlabel('Read count cutoff (-r, or the model\'s)', color=INK, fontsize=9)
    ax.set_ylabel('Junctions translated', color=INK, fontsize=9)
    ax.set_title('Junctions a run would translate at each cutoff', color=INK, fontsize=10, loc='left', pad=16)
    ax.legend(frameon=False, fontsize=8.5, labelcolor=INK, loc='upper right')
    line, label = _cutoff_line(ax, cutoff)
    move_chart_cutoff(line, label, cutoff)
    return line, label


def draw_frame(ax, frame, cutoff) -> tuple:
    """
    By read count, the share of events whose alternative region is a multiple of 3 nt (frame_by_count), annotated and
    novel, against one in three, and the cutoff (PR 27). A group's bin of fewer than FRAME_MIN_EVENTS is not drawn

    :return: the cutoff's line and label, which move_chart_cutoff moves
    """
    _style(ax)
    edges = FRAME_BINS + [FRAME_BINS[-1] * 2]
    x = [math.sqrt(lo * max(hi - 1, lo)) for lo, hi in zip(edges, edges[1:])]     # each bin's middle, on a log scale
    groups = [c[:-len(' events')] for c in frame.columns if c.endswith(' events')]
    for name in groups:
        color, marker = GROUPS[name]
        y = np.where(frame[f'{name} events'] >= FRAME_MIN_EVENTS, 100 * frame[f'{name} frame kept'], np.nan)
        ax.plot(x, y, color=color, linewidth=2, marker=marker, markersize=5, label=_label(name))
    ax.axhline(100 / 3, color=INK2, linewidth=1, linestyle=(0, (2, 2)), label='one in three: a random length')
    ax.set_xscale('log')
    ax.xaxis.set_major_formatter('{x:,.0f}')
    ax.set_ylim(0, 100)
    ax.set_xlabel('Read count', color=INK, fontsize=9)
    ax.set_ylabel('% of events', color=INK, fontsize=9)
    ax.set_title('Events whose alternative region is a multiple of 3 nt, by read count', color=INK, fontsize=10,
                 loc='left', pad=16)
    ax.legend(frameon=False, fontsize=8.5, labelcolor=INK, loc='upper right')
    line, label = _cutoff_line(ax, cutoff)
    move_chart_cutoff(line, label, cutoff)
    return line, label


def draw_model(ax, ln_sjc, params, min_count) -> tuple:
    """
    Draw the read count histogram, the two weighted components, their sum and the cutoff on ax, with a second x axis
    in reads. 2026-11-13 plot_general_mixture_model drew on a figure of its own; the GUI draws on its own Axes
    (REVIEW.md §10)

    :param ax: a matplotlib Axes
    :param ln_sjc: Natural log of read counts from rMATS.
    :param params: Output of general_mixture_model; None if there is no fit, for the histogram and the cutoff only
        (the GUI, with -r)
    :param min_count: The read count cutoff, from general_mixture_model or -r.
    :return: the cutoff's line and label, which move_cutoff moves (the GUI)
    """
    ax.hist(ln_sjc, bins=50, histtype='bar', density=True, ec='red', alpha=0.5)
    ax.set_xlabel('Natural logarithm of read count + 1')
    ax.set_ylabel('Frequency')
    ax.set_title('Junction read counts (smaller of summed inclusion and skipping junction reads)')
    ax.set_ylim(0, 1)

    if params is not None:
        shape, rate, mu, sd, w0, w1 = params
        x = np.linspace(0, 10, 2000)
        eval_dist0 = stats.gamma.pdf(x, shape, scale=1 / rate) * w0
        eval_dist1 = stats.norm.pdf(x, mu, sd) * w1
        ax.plot(x, eval_dist0, color="C0")
        ax.plot(x, eval_dist1, color="C3")
        ax.plot(x, eval_dist0 + eval_dist1, color="black")
        ax.fill_between(x, 0, eval_dist0, color="C0", alpha=0.2)
        ax.fill_between(x, 0, eval_dist1, color="C3", alpha=0.2)
    line = ax.axvline(linewidth=2, ls='--', color='black', x=np.log(min_count + 1))
    label = ax.text(0, 0.75, '')
    move_cutoff(line, label, min_count)

    # Second axis for untransformed reads
    tcks = ax.get_xticks()
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(tcks[1:-1])
    ax2.set_xlabel('Read count')
    ax2.set_xticklabels((np.exp(tcks[1:-1]) - 1).round().astype(int))
    return line, label


def move_cutoff(line, label, min_count) -> None:
    """ Put the cutoff's line and label, as draw_model returns them, at another read count """
    x = np.log(min_count + 1)
    line.set_xdata([x, x])
    label.set_x(x + math.log(1.25))
    label.set_text("cutoff = " + str(min_count) + " reads")


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
    # a Figure without pyplot, which keeps no state between runs and needs no GUI backend
    from matplotlib.figure import Figure

    fig = Figure(constrained_layout=True, figsize=(10, 4))
    draw_model(fig.subplots(), ln_sjc, params, min_count)
    fig.savefig(fname=os.path.join(write_dir, filename + '.pdf'))
