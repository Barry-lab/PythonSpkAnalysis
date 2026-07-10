"""Spatial binning and rate maps.

Ports ``histnd.m`` (the binning primitive), ``bin_data.m`` (position/spike
binning) and ``make_smooth_ratemap.m`` (smoothing with edge normalisation).
``rate_map`` bins occupancy and spikes, smooths each ignoring unvisited bins,
then divides. Smoothing is configurable: a boxcar (default, 5 bins — the lab
standard) or a Gaussian (by sigma), or none.

All functions consume a :class:`~barrelab.core.recording.Recording` or plain
arrays; none of them touch a file path.
"""

import numpy as np
from scipy.ndimage import convolve, gaussian_filter

from ..params import RATEMAP_BIN_CM, RATEMAP_GAUSS_SIGMA, RATEMAP_SMOOTH_BINS

_EPS = np.finfo(float).eps


def _histnd_2d(data_yx, ranges_yx, bin_size):
    """2-D histogram matching ``histnd.m``.

    Parameters
    ----------
    data_yx : ndarray (n, 2)
        Points as (row=y, col=x) — the same column order the MATLAB passes.
    ranges_yx : sequence (2,)
        Maximum permitted (y, x) values; the map spans ``[0, range]``.
    bin_size : float
        Square bin size in the same units as ``data_yx``.

    Returns
    -------
    ndarray (n_y_bins, n_x_bins)
        Bin counts. Bins are ``ceil(value / bin_size)`` (1-based internally);
        out-of-range points are clamped into the edge bins and rows containing
        NaN are dropped, exactly as in the original.
    """
    data_yx = np.asarray(data_yx, dtype=float)
    quotient = data_yx / bin_size
    # A value of exactly 0 would fall in bin 0; nudge it into bin 1.
    quotient = np.where(quotient == 0, quotient + _EPS, quotient)
    bins = np.ceil(quotient)  # NaNs preserved

    max_bin = np.ceil(np.asarray(ranges_yx, dtype=float) / bin_size).astype(int)

    valid = ~np.isnan(bins).any(axis=1)
    bins = bins[valid].astype(np.int64)
    bins[bins < 1] = 1
    bins[:, 0] = np.minimum(bins[:, 0], max_bin[0])
    bins[:, 1] = np.minimum(bins[:, 1], max_bin[1])

    out = np.zeros((max_bin[0], max_bin[1]))
    np.add.at(out, (bins[:, 0] - 1, bins[:, 1] - 1), 1.0)
    return out


def _bin_size_px(recording, bin_cm):
    return (recording.pixels_per_metre / 100.0) * bin_cm


def occupancy_map(recording, bin_cm=RATEMAP_BIN_CM):
    """Dwell time (seconds) per spatial bin."""
    bin_px = _bin_size_px(recording, bin_cm)
    counts = _histnd_2d(recording.xy[:, [1, 0]], recording.window_extent, bin_px)
    return counts / recording.pos_sample_rate


def spike_map(recording, pos_index, bin_cm=RATEMAP_BIN_CM):
    """Spike counts per spatial bin for spikes at the given position indices."""
    bin_px = _bin_size_px(recording, bin_cm)
    spike_xy = recording.xy[pos_index][:, [1, 0]]
    return _histnd_2d(spike_xy, recording.window_extent, bin_px)


def smooth_ratemap(bin_pos, bin_spk, kern_size=RATEMAP_SMOOTH_BINS,
                   kern_type="boxcar", smooth_type="norm"):
    """Smooth binned dwell + spikes into a rate map. Port of ``make_smooth_ratemap.m``.

    Occupancy and spike maps are smoothed independently with edge normalisation
    (each bin divided by the smoothed visited mask) so unvisited bins do not
    bleed in, then divided to give rate. Unvisited bins are set to NaN.

    Parameters
    ----------
    bin_pos, bin_spk : ndarray
        Binned dwell time (s) and spike counts, same shape.
    kern_size : float
        Smoothing amount, interpreted by ``kern_type``: for ``"boxcar"`` the box
        width in bins (5 is the lab standard); for ``"gaussian"`` the standard
        deviation in bins.
    kern_type : {"boxcar", "gaussian"}
        Smoothing kernel.
    smooth_type : {"norm", "circular"}
        ``"norm"`` (default) for ordinary maps; ``"circular"`` wraps the edges,
        for 1-D circular variables such as head direction.

    Returns
    -------
    (rate, smooth_pos, smooth_spk)
    """
    bin_pos = np.asarray(bin_pos, dtype=float)
    bin_spk = np.asarray(bin_spk, dtype=float)
    mode = "wrap" if smooth_type == "circular" else "constant"

    if kern_type == "boxcar":
        # imfilter's default boundary is zero-padding; a symmetric ones kernel
        # convolved with mode="constant" reproduces it (verified vs MATLAB).
        kernel = np.ones([int(round(kern_size))] * bin_pos.ndim)
        smooth = lambda a: convolve(a, kernel, mode=mode, cval=0.0)  # noqa: E731
    elif kern_type == "gaussian":
        # truncated at 4 sigma (scipy default) — numerically identical to
        # MATLAB's full-map fspecial('gaussian', ...) for any usable sigma.
        smooth = lambda a: gaussian_filter(a, sigma=kern_size, mode=mode, cval=0.0)  # noqa: E731
    else:
        raise ValueError("kern_type must be 'boxcar' or 'gaussian', got %r" % (kern_type,))

    visited = (bin_pos != 0).astype(float)
    denom = smooth(visited)
    with np.errstate(divide="ignore", invalid="ignore"):
        smooth_pos = smooth(bin_pos) / denom
        smooth_spk = smooth(bin_spk) / denom
        rate = smooth_spk / smooth_pos

    unvisited = visited == 0
    rate[unvisited] = np.nan
    smooth_pos[unvisited] = np.nan
    smooth_spk[unvisited] = np.nan
    return rate, smooth_pos, smooth_spk


def rate_map(recording, electrode, cluster, bin_cm=RATEMAP_BIN_CM,
             smoothing="boxcar", smooth_size=None):
    """Firing-rate map (Hz) for one cluster.

    With the defaults this reproduces ``easy_make_ratemap.m`` (a 5-bin boxcar).

    Parameters
    ----------
    recording : Recording
    electrode : int
        Tetrode/group id.
    cluster : int
        Cluster id within that electrode.
    bin_cm : float
        Spatial bin size in cm.
    smoothing : {"boxcar", "gaussian", None}
        Smoothing method. ``None`` returns the raw (unsmoothed) rate map.
    smooth_size : float, optional
        Smoothing amount: boxcar width in bins, or gaussian sigma in bins.
        Defaults to 5 bins for boxcar and 2 bins (sigma) for gaussian.

    Returns
    -------
    ndarray (n_y_bins, n_x_bins) of firing rate in Hz, NaN where unvisited.
    """
    elec = recording.electrode(electrode)
    occ = occupancy_map(recording, bin_cm)
    spk = spike_map(recording, elec.spike_pos_index(cluster), bin_cm)

    if smoothing is None:
        with np.errstate(divide="ignore", invalid="ignore"):
            rate = spk / occ
        rate[occ == 0] = np.nan
        return rate

    if smooth_size is None:
        smooth_size = RATEMAP_SMOOTH_BINS if smoothing == "boxcar" else RATEMAP_GAUSS_SIGMA
    rate, _, _ = smooth_ratemap(occ, spk, smooth_size, smoothing)
    return rate
