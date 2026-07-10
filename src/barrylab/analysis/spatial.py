"""Spatial-coding metrics.

Ports ``skaggsInfo.m`` (Skaggs et al. spatial information) and
``skaggsAdaptiveBinning.m`` (adaptive smoothing, the recommended precursor to a
spatial-information estimate).
"""

import numpy as np
from scipy.ndimage import binary_fill_holes, convolve, label

from ..params import ADAPTIVE_SMOOTH_ALPHA


def spatial_information(rates, dwell_times):
    """Skaggs spatial information.

        I = sum_x p(x) * r(x) * log2( r(x) / mean_r )

    Parameters
    ----------
    rates : ndarray
        Firing rate per bin (Hz). May contain NaN for unvisited bins.
    dwell_times : ndarray
        Occupancy per bin (seconds), same shape as ``rates``.

    Returns
    -------
    (bits_per_spike, bits_per_second)
        ``(nan, nan)`` if the cell never fires.

    Notes
    -----
    Faithful port of ``skaggsInfo.m``. For a well-behaved estimate feed it an
    adaptively-binned rate map (see :func:`adaptive_binning`); a boxcar-smoothed
    map can give spurious values.
    """
    rates = np.asarray(rates, dtype=float)
    dwell_times = np.asarray(dwell_times, dtype=float)
    if rates.shape != dwell_times.shape:
        raise ValueError("rates and dwell_times must have the same shape")
    if np.nanmax(rates) <= 0:
        return np.nan, np.nan

    rates = np.where(np.isnan(rates), 0.0, rates)
    dwell_times = np.where(np.isnan(dwell_times), 0.0, dwell_times)

    duration = dwell_times.sum()
    mean_rate = np.sum(rates * dwell_times) / duration

    p_x = dwell_times / duration
    p_r = rates / mean_rate
    weighted = p_x * rates
    active = weighted > 0
    bits_per_second = np.sum(weighted[active] * np.log2(p_r[active]))
    bits_per_spike = bits_per_second / mean_rate
    return bits_per_spike, bits_per_second


def adaptive_binning(dwell_times, spikes, alpha=ADAPTIVE_SMOOTH_ALPHA):
    """Skaggs adaptive smoothing of a binned map. Port of ``skaggsAdaptiveBinning.m``.

    Each bin is smoothed with the smallest disc whose radius ``r`` satisfies
    ``alpha / (sqrt(n_spikes) * dwell) <= r``. Returns ``(rate, smooth_spk,
    smooth_pos)`` with unvisited bins as NaN.
    """
    dwell_times = np.asarray(dwell_times, dtype=float)
    spikes = np.asarray(spikes, dtype=float)

    if spikes.sum() == 0:
        empty = np.where(dwell_times == 0, np.nan, dwell_times)
        return (np.where(dwell_times == 0, np.nan, spikes), empty, empty)

    smooth_pos = np.zeros_like(dwell_times)
    smooth_spk = np.zeros_like(dwell_times)
    visited = (dwell_times > 0).astype(float)

    done = dwell_times == 0  # unvisited bins need no smoothing
    is_1d = dwell_times.ndim == 1 or 1 in dwell_times.shape
    boundary = "wrap" if is_1d else "constant"
    r_bump = 1.0 if is_1d else 0.5

    r = 1.0
    while not done.all():
        if r > max(dwell_times.shape) / 2.0:
            unfinished = ~done
            smooth_spk[unfinished] = np.nan
            smooth_pos[unfinished] = np.nan
            break

        kernel = _disc_kernel(r, is_1d)
        f_spk = convolve(spikes, kernel, mode=boundary, cval=0.0)
        f_pos = convolve(dwell_times, kernel, mode=boundary, cval=0.0)
        f_vis = convolve(visited, kernel, mode=boundary, cval=0.0)

        with np.errstate(divide="ignore", invalid="ignore"):
            passed = alpha / (np.sqrt(f_spk) * f_pos) <= r
        passed &= ~done
        smooth_spk[passed] = f_spk[passed] / f_vis[passed]
        smooth_pos[passed] = f_pos[passed] / f_vis[passed]
        done |= passed
        r += r_bump

    with np.errstate(divide="ignore", invalid="ignore"):
        rate = smooth_spk / smooth_pos
    unvisited = dwell_times == 0
    rate[unvisited] = np.nan
    smooth_pos[unvisited] = np.nan
    smooth_spk[unvisited] = np.nan
    return rate, smooth_spk, smooth_pos


def _disc_kernel(r, is_1d):
    """Flat disc (2-D) or boxcar (1-D) kernel, matching ``fspecial('disk', r)`` usage."""
    if is_1d:
        return np.ones(int(1 + r * 2))
    n = int(np.ceil(r))
    yy, xx = np.mgrid[-n:n + 1, -n:n + 1]
    disc = (np.sqrt(xx ** 2 + yy ** 2) <= r).astype(float)
    return disc


def true_ppm_bounding_box(occupancy_bins, bin_size_px, env_size_m):
    """Infer the true pixels-per-metre from a square environment's occupied area.

    Port of ``get_true_PPM_boundingBox.m``. The `.set` file's recorded PPM is
    sometimes wrong; if the environment is square and its side length is known,
    the true PPM can be inferred from how many bins the animal actually visited.
    Finds the largest occupied (4-connected) region, fills any interior holes,
    and derives PPM from ``sqrt(filled_area) * bin_size_px / env_size_m``.

    Only reliable if the environment is well sampled — under-visited corners
    make the visited area (and so the inferred PPM) an overestimate.

    Parameters
    ----------
    occupancy_bins : ndarray
        Binned dwell time/counts (e.g. from :func:`~barrylab.analysis.maps.occupancy_map`);
        a bin is "visited" if its value is > 0.
    bin_size_px : float
        Bin size in camera pixels, as used to build ``occupancy_bins``.
    env_size_m : float
        Side length of the (square) environment, in metres.

    Returns
    -------
    (ppm, i_range, j_range)
        ``ppm`` is the inferred pixels-per-metre. ``i_range``/``j_range`` are
        0-based Python ``range`` objects over the occupied region's rows/cols
        (for cropping). **Note:** these faithfully reproduce an off-by-one in
        the original MATLAB that excludes the region's last row and column —
        see ``FINDINGS.md``. Returns ``(nan, range(0,0), range(0,0))`` if no
        bin was ever visited.
    """
    visited = np.asarray(occupancy_bins) > 0
    labelled, n = label(visited)  # scipy default = 4-connectivity in 2D
    if n == 0:
        return float("nan"), range(0, 0), range(0, 0)

    filled_areas = np.empty(n)
    bboxes = []
    for i in range(1, n + 1):
        mask = labelled == i
        filled_areas[i - 1] = binary_fill_holes(mask).sum()
        rows, cols = np.where(mask)
        bboxes.append((rows.min(), rows.max(), cols.min(), cols.max()))

    best = int(np.argmax(filled_areas))
    filled_area = filled_areas[best]
    r0, r1, c0, c1 = bboxes[best]

    ppm = float(np.sqrt(filled_area) * bin_size_px / env_size_m)
    # MATLAB: iRange = ijCorner(1):ijCorner(1)+ijWidth(1)-1, which (after
    # translating the 1-based BoundingBox arithmetic) excludes the region's
    # final row/col. Reproduced here rather than "fixed" — see FINDINGS.md.
    return ppm, range(r0, r1), range(c0, c1)
