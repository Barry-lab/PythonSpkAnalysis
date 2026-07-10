"""Bayesian population decoding of position from spike counts.

Ports the core of ``Decoding/``:

- :func:`spike_counts_per_window` — spike counts per (optionally overlapping,
  optionally Gaussian-smoothed) time window (``decode_getSpkCnt.m``).
- :func:`bayesian_posterior` — naive maximum-likelihood Poisson population
  decoder (``decode_calcPost.m``).
- :func:`bayesian_posterior_with_prior` — full Bayesian decoder with a spatial
  occupancy prior and a sequential movement (closeness) prior
  (``decode_calcBayesPost.m``, Zhang et al. 1998 eq. 36).
- :func:`decode_maximum_likelihood` — turn a posterior into decoded positions
  and measure accuracy against ground truth (``decode_processPost.m``).
- :func:`crunch_decode` — a Recording-level convenience orchestrator (a
  redesign of ``decode_crunchData.m`` around :class:`~barrylab.core.recording.Recording`;
  validate its pieces individually, this end-to-end wrapper is not bit-exact).
"""

import numpy as np
from scipy.ndimage import convolve1d
from scipy.special import gamma

from ._binning import histnd
from .maps import smooth_ratemap, spike_map

# "Small value" added to expected spike counts to avoid literal zeros (a rate
# map bin at exactly 0 Hz makes a single observed spike there impossible,
# i.e. probability 0 -- not useful). decode_calcPost.m and decode_calcBayesPost.m
# share the identical originating comment for this, but as of the current
# upstream MATLAB they've DIVERGED: decode_calcPost.m was fixed in 2020
# (eps^8 -> eps, since eps^8 is so vanishingly tiny it doesn't actually
# prevent the zero-probability problem the comment describes) but that fix
# was never propagated to decode_calcBayesPost.m, which still uses the old,
# ineffective eps^8. Kept as two separate constants matching each function's
# actual current MATLAB source, rather than one shared value -- see FINDINGS.md.
_SMALL_V_CALC_POST = np.finfo(float).eps       # decode_calcPost.m (fixed upstream 2020)
_SMALL_V_BAYES_POST = np.finfo(float).eps ** 8  # decode_calcBayesPost.m (not yet fixed upstream)


def _gaussian_kernel_1d(sigma):
    """Matches MATLAB ``fspecial('gaussian', [1, ceil(5*sigma)], sigma)``."""
    n = int(np.ceil(sigma * 5))
    x = np.arange(n) - (n - 1) / 2.0
    h = np.exp(-(x ** 2) / (2 * sigma ** 2))
    h[h < np.finfo(float).eps * h.max()] = 0.0
    return h / h.sum()


def spike_counts_per_window(spike_pos, start_end_pos, win_width,
                            win_overlap=0.0, win_smooth=0.0):
    """Spike counts per time window, with optional overlap/smoothing.

    Port of ``decode_getSpkCnt.m``. All positions/widths are in position-sample
    units (not seconds), 0-based, matching :attr:`Electrode.pos_index`.

    Parameters
    ----------
    spike_pos : ndarray
        Position-sample index of each spike (repeats allowed).
    start_end_pos : (int, int)
        First and last position sample to consider (inclusive).
    win_width : int
        Window width in position samples.
    win_overlap : float
        Overlap between neighbouring windows as a proportion of ``win_width``
        (``0`` = no overlap, must be ``< 1``).
    win_smooth : float
        If > 0, the sigma (in position samples) of a Gaussian used to smooth
        spike counts before windowing (so counts become non-integer).

    Returns
    -------
    (spike_count, win_lim)
        ``spike_count`` — ``(n_win,)``, possibly fractional if smoothed.
        ``win_lim`` — ``(n_win, 2)`` int, the ``[start, end]`` position sample
        (inclusive) of each window.
    """
    if win_overlap >= 1:
        raise ValueError("win_overlap must be less than 1")
    spike_pos = np.asarray(spike_pos, dtype=float).ravel()

    win_overlap_samp = int(round(win_width * win_overlap))
    n_pos_samp = start_end_pos[1] - start_end_pos[0] + 1
    win_advance = win_width - win_overlap_samp
    n_win = int(np.floor((n_pos_samp - win_width) / win_advance)) + 1

    win_idx = np.arange(n_win)
    win_lim = np.round(
        np.column_stack([win_idx, win_idx]) * win_advance
        + np.array([0, win_width - 1])
        + start_end_pos[0]
    ).astype(np.int64)
    win_pos = win_lim[:, [0]] + np.arange(win_width)  # (n_win, win_width)

    fine_start, fine_end = int(win_lim[0, 0]), int(win_lim[-1, 1])
    n_fine = fine_end - fine_start + 1
    fine_counts = np.zeros(n_fine)
    fine_idx = np.floor(spike_pos - fine_start).astype(np.int64)
    valid = (fine_idx >= 0) & (fine_idx < n_fine)
    np.add.at(fine_counts, fine_idx[valid], 1.0)

    if win_smooth > 0:
        kernel = _gaussian_kernel_1d(win_smooth)
        smoothed = convolve1d(fine_counts, kernel, mode="constant", cval=0.0)
        denom = convolve1d(np.ones(n_fine), kernel, mode="constant", cval=0.0)
        fine_counts = smoothed / denom

    fine_index = win_pos - fine_start
    spike_count = fine_counts[fine_index].sum(axis=1)
    return spike_count, win_lim


def bayesian_posterior(obs_spk, rate_maps, t_bin, normalise=True):
    """Naive maximum-likelihood Poisson population decoder.

    Port of ``decode_calcPost.m``. For each time bin, compares observed spike
    counts across a population of cells against their rate maps and returns
    the (Poisson) posterior probability of the animal being at each location.

    Parameters
    ----------
    obs_spk : ndarray (n_cells, n_time_bins)
        Observed spike count per cell per time bin (need not be integer).
    rate_maps : ndarray (n_y, n_x, n_cells)
        Stacked, smoothed rate maps (Hz); unvisited bins should be NaN.
    t_bin : float
        Time window width in seconds.
    normalise : bool
        If ``True`` (default), each time bin's posterior sums to 1.

    Returns
    -------
    ndarray (n_y, n_x, n_time_bins)
    """
    obs_spk = np.asarray(obs_spk, dtype=float)
    rate_maps = np.asarray(rate_maps, dtype=float)
    n_cells = obs_spk.shape[0]
    if n_cells != rate_maps.shape[2]:
        raise ValueError("number of cells in obs_spk and rate_maps must match")

    expected = (rate_maps + _SMALL_V_CALC_POST) * t_bin  # (ny, nx, n_cells)
    expon = np.exp(-expected)
    factorial = gamma(obs_spk + 1)  # (n_cells, n_time)

    working = (expected[:, :, :, None] ** obs_spk[None, None, :, :]) / factorial[None, None, :, :]
    working = working * expon[:, :, :, None]
    post = np.prod(working, axis=2)  # product over cells -> (ny, nx, n_time)

    if normalise:
        post = post / np.nansum(post, axis=(0, 1), keepdims=True)
    return post


def bayesian_posterior_with_prior(obs_spk, rate_maps, dwell_map, t_bin,
                                  n_step_for_vel=15, speed_scale_fact=5.0):
    """Full Bayesian decoder: spike likelihood x spatial prior x movement prior.

    Port of ``decode_calcBayesPost.m`` (Zhang et al. 1998, eq. 36 + a sequential
    "closeness" prior). Unlike :func:`bayesian_posterior`, this also uses the
    overall occupancy (as a spatial prior) and assumes the animal is unlikely
    to have moved far from its own just-decoded location — applied as a 2-D
    Gaussian centred on the previous decode, with sigma set from the recent
    decoded-position speed. Because each step depends on the previous step's
    decode, this runs a genuine sequential loop (not vectorisable over time).

    Parameters
    ----------
    obs_spk : ndarray (n_cells, n_time_bins)
    rate_maps : ndarray (n_y, n_x, n_cells)
    dwell_map : ndarray (n_y, n_x)
        Smoothed dwell time (s) per bin, NaN where never visited.
    t_bin : float
        Time window width in seconds.
    n_step_for_vel : int
        Number of preceding decoded steps used to estimate speed for the
        movement prior (also how many initial steps use a flat prior).
    speed_scale_fact : float
        Scales the inferred speed into the movement prior's Gaussian sigma.

    Returns
    -------
    ndarray (n_y, n_x, n_time_bins)
    """
    obs_spk = np.asarray(obs_spk, dtype=float)
    rate_maps = np.asarray(rate_maps, dtype=float)
    dwell_map = np.asarray(dwell_map, dtype=float)
    n_cells = obs_spk.shape[0]
    if n_cells != rate_maps.shape[2]:
        raise ValueError("number of cells in obs_spk and rate_maps must match")

    px = dwell_map / np.nansum(dwell_map)

    expected = (rate_maps + _SMALL_V_BAYES_POST) * t_bin
    expon = np.exp(-expected)
    working = (expected[:, :, :, None] ** obs_spk[None, None, :, :]) * expon[:, :, :, None]
    post = np.prod(working, axis=2)
    post = post / np.nansum(post, axis=(0, 1), keepdims=True)
    post = post * px[:, :, None]
    post = post / np.nansum(post, axis=(0, 1), keepdims=True)

    ny, nx = dwell_map.shape
    xx, yy = np.meshgrid(np.arange(nx), np.arange(ny))
    last_few_pos = np.zeros((n_step_for_vel, 2))
    p_new_x_given_x = np.ones((ny, nx)) / (ny * nx)
    i_s = j_s = sigma = None

    n_time = post.shape[2]
    for n in range(n_time):
        post[:, :, n] = post[:, :, n] * p_new_x_given_x

        if n >= n_step_for_vel:
            p_new_x_given_x = np.exp(-(((xx - j_s) ** 2 + (yy - i_s) ** 2) / (2 * sigma ** 2)))
            p_new_x_given_x = p_new_x_given_x / p_new_x_given_x.sum()

        # ML decode of this step's peak. MATLAB's max() ignores NaN (unless
        # everything is NaN); column-major ('F') order matches its tie-break.
        flat_idx = np.nanargmax(post[:, :, n].ravel(order="F"))
        i_s, j_s = np.unravel_index(flat_idx, (ny, nx), order="F")
        last_few_pos = np.vstack([last_few_pos[1:], [i_s, j_s]])
        step_dist = np.sqrt(np.sum(np.diff(last_few_pos, axis=0) ** 2, axis=1))
        speed = step_dist.mean()
        sigma = speed * speed_scale_fact
        if sigma == 0:
            sigma = np.finfo(float).eps  # avoid a degenerate zero-width gaussian

    post = post / np.nansum(post, axis=(0, 1), keepdims=True)
    return post


def decode_maximum_likelihood(post, recording, bin_px, win_lim):
    """Decode ML position from a posterior and measure accuracy.

    Port of ``decode_processPost.m``. Time bins whose posterior is all-NaN are
    left as NaN in the output (use this to skip specific windows).

    .. note::
       The original MATLAB's "if several bins are jointly maximal, pick one at
       random" is not actually implemented (its random-index line operates on
       an already-scalar index from ``max()``, which is a no-op) — ties always
       resolve to the first maximum in MATLAB's column-major order. Reproduced
       here via ``order='F'`` rather than genuine random tie-breaking, per
       goal 1 (see ``FINDINGS.md``).

    Parameters
    ----------
    post : ndarray (n_y, n_x, n_time_bins)
        Posterior from :func:`bayesian_posterior` or
        :func:`bayesian_posterior_with_prior`.
    recording : Recording
        Source of ground-truth position (``recording.xy``).
    bin_px : float
        Spatial bin size in pixels (matching what built ``post``'s grid).
    win_lim : ndarray (n_time_bins, 2)
        ``[start, end]`` position sample (inclusive, 0-based) per time window,
        from :func:`spike_counts_per_window`.

    Returns
    -------
    (err, dec_xy, true_xy, ml_val)
        ``err`` — Euclidean decode error in pixels, ``(n_time_bins,)``.
        ``dec_xy``/``true_xy`` — decoded/true bin-centre position in pixels,
        ``(n_time_bins, 2)``.
        ``ml_val`` — the posterior's peak value per time bin.
    """
    win_lim = np.asarray(win_lim, dtype=np.int64)
    n_win = win_lim.shape[0]
    ny, nx = post.shape[0], post.shape[1]

    win_width = int(win_lim[0, 1] - win_lim[0, 0] + 1)
    win_full_pos_ind = win_lim[:, [0]] + np.arange(win_width)  # (n_win, win_width)

    true_pos_x = recording.xy[win_full_pos_ind, 0].mean(axis=1)
    true_pos_y = recording.xy[win_full_pos_ind, 1].mean(axis=1)
    true_xy_raw = np.column_stack([true_pos_x, true_pos_y])

    window_extent = (ny * bin_px, nx * bin_px)
    _, bins = histnd(true_xy_raw[:, [1, 0]], np.array(window_extent), bin_px)

    xx, yy = np.meshgrid((np.arange(nx) + 0.5) * bin_px, (np.arange(ny) + 0.5) * bin_px)

    true_xy = np.column_stack([xx[bins[:, 0], bins[:, 1]], yy[bins[:, 0], bins[:, 1]]])

    ml_val = np.full(n_win, np.nan)
    dec_xy = np.full((n_win, 2), np.nan)
    for n in range(n_win):
        cur_post = post[:, :, n]
        if np.all(np.isnan(cur_post)):
            continue
        flat_idx = np.nanargmax(cur_post.ravel(order="F"))
        row, col = np.unravel_index(flat_idx, (ny, nx), order="F")
        ml_val[n] = cur_post[row, col]
        dec_xy[n] = [xx[row, col], yy[row, col]]

    err = np.sqrt(np.sum((dec_xy - true_xy) ** 2, axis=1))
    return err, dec_xy, true_xy, ml_val


def crunch_decode(recording, cells, t_bin=0.5, test_fraction=0.1,
                  win_overlap=0.0, win_smooth=5.0, bin_cm=2.0, seed=None):
    """End-to-end decode: build rate maps on a train split, decode a test split.

    A Recording-level redesign of ``decode_crunchData.m`` (not a bit-exact
    port — validate its component functions individually; see
    ``PORTING_PLAN.md``). Splits position samples into contiguous train/test
    blocks, builds rate maps from the train block, decodes the test block, and
    returns per-window decode error.

    Parameters
    ----------
    recording : Recording
    cells : list of (electrode, cluster)
        The population of cells to decode with.
    t_bin : float
        Decoding time window, in seconds.
    test_fraction : float
        Fraction of the session (a contiguous trailing block) held out for
        testing; the rest builds the rate maps.
    win_overlap : float
        Time-window overlap (0 to <1) for the spike-count windows.
    win_smooth : float
        Gaussian sigma (position samples) for spike-count smoothing.
    bin_cm : float
        Spatial bin size in cm for both rate maps and the decode grid.
    seed : int, optional
        Unused (kept for API symmetry with a future randomised split);
        the current split is deterministic (contiguous trailing block).

    Returns
    -------
    dict with ``err`` (cm), ``dec_xy``/``true_xy`` (cm), ``ml_val``, ``n_cells``.
    """
    n_pos = recording.n_pos
    n_test = int(round(n_pos * test_fraction))
    train_end = n_pos - n_test
    train_idx = np.arange(0, train_end)
    test_start, test_end = train_end, n_pos - 1

    bin_px = (recording.pixels_per_metre / 100.0) * bin_cm
    win_width_samp = int(round(recording.pos_sample_rate * t_bin))

    train_mask = np.zeros(n_pos, dtype=bool)
    train_mask[train_idx] = True

    # Occupancy restricted to the train split too (spike_map bins an arbitrary
    # position-index subset, so it doubles as a train-only dwell-count map).
    occ_train = spike_map(recording, train_idx, bin_cm) / recording.pos_sample_rate

    rate_maps = []
    spike_counts = []
    win_lim = None
    for electrode, cluster in cells:
        elec = recording.electrode(electrode)
        spk_idx_all = elec.spike_pos_index(cluster)
        spk_idx_train = spk_idx_all[train_mask[spk_idx_all]]

        spk_map = spike_map(recording, spk_idx_train, bin_cm)
        rm, _, _ = smooth_ratemap(occ_train, spk_map, 5, "boxcar")
        rate_maps.append(rm)

        counts, win_lim = spike_counts_per_window(
            spk_idx_all, (test_start, test_end), win_width_samp, win_overlap, win_smooth,
        )
        spike_counts.append(counts)

    rate_maps = np.stack(rate_maps, axis=2)
    spike_counts = np.stack(spike_counts, axis=0)

    post = bayesian_posterior(spike_counts, rate_maps, t_bin)
    err_px, dec_xy_px, true_xy_px, ml_val = decode_maximum_likelihood(post, recording, bin_px, win_lim)

    ppm = recording.pixels_per_metre
    return {
        "err": err_px / ppm * 100.0,
        "dec_xy": dec_xy_px / ppm * 100.0,
        "true_xy": true_xy_px / ppm * 100.0,
        "ml_val": ml_val,
        "n_cells": len(cells),
    }
