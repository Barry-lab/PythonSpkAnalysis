"""Place x direction (pxd) joint firing model.

Ports ``pxd.m`` and ``wrapper_pxd.m``: a maximum-likelihood factorial model
separating a cell's place (p) and head-direction (d) tuning under a Poisson
firing assumption (Burgess, Cacucci, Lever & O'Keefe, Hippocampus 2005). Unlike
a plain rate map (which conflates place and direction effects when the two are
correlated during behaviour), this fits them jointly: firing rate at
``(direction j, place bin i)`` is modelled as ``p(i) * d(j)``, estimated by
alternating maximum-likelihood updates of ``p`` given ``d`` and vice versa
until the log-likelihood converges.
"""

import warnings

import numpy as np
from scipy.special import gammaln

from ..params import (
    PXD_ACCURACY,
    PXD_DIR_SMOOTH_BINS,
    PXD_MAX_ITER,
    PXD_POS_BIN_CM,
    PXD_POS_SMOOTH_BINS,
    PXD_TOL,
    DIR_BIN_DEG,
)
from ._binning import histnd
from .maps import smooth_ratemap


def _p_estimate(d, spikes, times, tol):
    """p(i1,i2) = sum_j spikes(j,i1,i2) / sum_j(d(j) * times(j,i1,i2))."""
    denom = np.sum(d[:, None, None] * times, axis=0)
    denom = np.where(denom < tol, 1.0 / tol, denom)
    return np.sum(spikes, axis=0) / denom


def _d_estimate(p, spikes, times, tol):
    """d(j) = sum_i1i2 spikes(j,i1,i2) / sum_i1i2(p(i1,i2) * times(j,i1,i2))."""
    denom = np.sum(p[None, :, :] * times, axis=(1, 2))
    denom = np.where(denom < tol, 1.0 / tol, denom)
    return np.sum(spikes, axis=(1, 2)) / denom


def _loglikelihood(p, d, spikes, times, tol):
    """Poisson log-likelihood of spikes under expected = p(i1,i2)*d(j)*times(j,i1,i2)."""
    expected = p[None, :, :] * d[:, None, None] * times
    expected_clamped = np.where(expected < tol, tol, expected)
    return np.sum(spikes * np.log(expected_clamped) - expected - gammaln(spikes + 1))


def pxd(spikes, times, max_iter=PXD_MAX_ITER, accuracy=PXD_ACCURACY, tol=PXD_TOL):
    """Fit the place x direction ML factorial model to pre-binned data.

    Port of ``pxd.m``. Alternates estimating ``d`` given ``p`` and ``p`` given
    ``d`` (Poisson maximum likelihood), stopping when the fractional change in
    log-likelihood drops below ``accuracy`` or ``max_iter`` is reached.

    Parameters
    ----------
    spikes : ndarray (n_dir_bins, n_y_bins, n_x_bins)
        Spike counts per (direction, position) bin, e.g. from :func:`~barrylab.analysis._binning.histnd`.
    times : ndarray, same shape as ``spikes``
        Dwell time in **seconds** per bin.
    max_iter, accuracy, tol : see :mod:`barrylab.params` (``PXD_*``).

    Returns
    -------
    (p, d)
        ``p`` — unsmoothed positional firing map, shape ``(n_y_bins, n_x_bins)``.
        ``d`` — unsmoothed directional firing map, shape ``(n_dir_bins,)``.
        Both are normalised so their predicted total spike count matches the
        actual total. Returns ``(nan, nan)`` (with a warning) if the fit does
        not converge within ``max_iter`` iterations.
    """
    spikes = np.asarray(spikes, dtype=float)
    times = np.asarray(times, dtype=float)
    if spikes.shape != times.shape:
        raise ValueError("spikes and times must have the same shape")

    n_dir, n_y, n_x = spikes.shape
    p = np.ones((n_y, n_x))
    d = np.ones(n_dir)
    fit = 1.0
    converged = False

    for _ in range(max_iter):
        d = _d_estimate(p, spikes, times, tol)
        p = _p_estimate(d, spikes, times, tol)
        prev_fit = fit
        fit = _loglikelihood(p, d, spikes, times, tol)
        if abs(prev_fit - fit) < -accuracy * fit:
            converged = True
            break

    if not converged:
        warnings.warn("pxd did not converge, returning nan")
        return float("nan"), float("nan")

    total_spikes = spikes.sum()
    time_per_pos_bin = times.sum(axis=0)
    predicted = np.sum(p * time_per_pos_bin)
    p = p * (total_spikes / predicted)

    time_per_dir_bin = times.sum(axis=(1, 2))
    predicted_dir = np.sum(d * time_per_dir_bin)
    d = d * (total_spikes / predicted_dir)

    return p, d


def place_by_direction(recording, electrode, cluster,
                       bin_cm=PXD_POS_BIN_CM, dir_bin_deg=DIR_BIN_DEG,
                       smooth=True, pos_smooth_bins=PXD_POS_SMOOTH_BINS,
                       dir_smooth_bins=PXD_DIR_SMOOTH_BINS,
                       max_iter=PXD_MAX_ITER, accuracy=PXD_ACCURACY, tol=PXD_TOL):
    """Bin one cluster's spikes by (direction, position) and fit the pxd model.

    Port of ``wrapper_pxd.m``: bins raw position/direction/spike data into the
    3-D histogram :func:`pxd` needs, fits the model, and (by default) smooths
    the two marginals the same way the original did — a boxcar over the
    positional map (edge-corrected so unvisited bins don't bleed in) and a
    plain circular boxcar over the directional map (no edge correction; the
    original assumes every direction bin has some dwell time).

    .. note::
       Faithfully reproduces a bug in ``wrapper_pxd.m``: the spatial bin range
       is computed as ``[max(y), max(y)]`` instead of ``[max(y), max(x)]`` (a
       copy-paste error — the min-value row is correct). For a non-square
       tracked window this mis-scales the X extent. See ``FINDINGS.md``.

    Parameters
    ----------
    recording : Recording
    electrode, cluster : int
    bin_cm : float
        Spatial bin size in cm (coarser than a plain rate map — see
        :mod:`barrylab.params` ``PXD_POS_BIN_CM`` docstring).
    dir_bin_deg : float
        Directional bin width in degrees.
    smooth : bool
        If ``True`` (default) return the smoothed maps; if ``False``, the raw
        ML fit.
    pos_smooth_bins, dir_smooth_bins : int
        Boxcar widths used when ``smooth=True``.
    max_iter, accuracy, tol : see :func:`pxd`.

    Returns
    -------
    (p, d)
        ``p`` is ``(n_y_bins, n_x_bins)`` Hz, NaN where unvisited. ``d`` is
        ``(n_dir_bins,)`` Hz.
    """
    elec = recording.electrode(electrode)
    spk_idx = elec.spike_pos_index(cluster)

    bin_px = (recording.pixels_per_metre / 100.0) * bin_cm
    max_y = float(np.max(recording.xy[:, 1]))
    min_y = float(np.min(recording.xy[:, 1]))
    min_x = float(np.min(recording.xy[:, 0]))
    # NB: reproduces wrapper_pxd.m's bug -- max_y used for BOTH the y and x
    # max range (should be max_x for the third column); see docstring above.
    data_ranges = np.array([[360.0, max_y, max_y], [0.0, min_y, min_x]])
    bin_size = np.array([dir_bin_deg, bin_px, bin_px])

    spk_dir_pos = np.column_stack([
        recording.direction[spk_idx], recording.xy[spk_idx, 1], recording.xy[spk_idx, 0],
    ])
    all_dir_pos = np.column_stack([recording.direction, recording.xy[:, 1], recording.xy[:, 0]])

    spike_bin, _ = histnd(spk_dir_pos, data_ranges, bin_size)
    occ_bin, _ = histnd(all_dir_pos, data_ranges, bin_size)
    occ_bin = occ_bin / recording.pos_sample_rate  # seconds

    p, d = pxd(spike_bin, occ_bin, max_iter=max_iter, accuracy=accuracy, tol=tol)
    if not smooth or (np.isscalar(p) and np.isnan(p)):
        return p, d

    visited = (occ_bin.sum(axis=0) != 0).astype(float)  # any dwell time at this (y,x)?
    smooth_p, _, _ = smooth_ratemap(visited, p * visited, pos_smooth_bins, "boxcar", "norm")

    always_visited = np.ones_like(d)
    smooth_d, _, _ = smooth_ratemap(always_visited, d, dir_smooth_bins, "boxcar", "circular")

    return smooth_p, smooth_d
