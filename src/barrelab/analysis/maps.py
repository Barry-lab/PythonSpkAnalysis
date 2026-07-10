"""Spatial binning and rate maps.

Ports ``histnd.m`` (the binning primitive), ``bin_data.m`` (position/spike
binning) and ``make_smooth_ratemap.m`` (boxcar smoothing with edge
normalisation). ``rate_map`` reproduces the historic ``easy_make_ratemap.m``
pipeline: bin occupancy and spikes, smooth each with a 5-bin boxcar ignoring
unvisited bins, then divide.

All functions consume a :class:`~barrelab.core.recording.Recording` or plain
arrays; none of them touch a file path.
"""

import numpy as np
from scipy.ndimage import convolve

from ..params import RATEMAP_BIN_CM, RATEMAP_SMOOTH_BINS

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


def smooth_ratemap(bin_pos, bin_spk, kern_size=RATEMAP_SMOOTH_BINS, kern_type="boxcar"):
    """Smooth binned dwell + spikes into a rate map. Port of ``make_smooth_ratemap.m``.

    Occupancy and spike maps are smoothed independently with edge normalisation
    (each bin divided by the smoothed visited-mask) so that unvisited bins do
    not bleed in, then divided to give rate. Unvisited bins are set to NaN.

    Returns ``(rate, smooth_pos, smooth_spk)``.
    """
    bin_pos = np.asarray(bin_pos, dtype=float)
    bin_spk = np.asarray(bin_spk, dtype=float)

    if kern_type == "boxcar":
        kernel = np.ones((kern_size, kern_size))
    elif kern_type == "gaussian":
        kernel = _gaussian_kernel(bin_pos.shape, kern_size)
    else:
        raise ValueError("unknown kern_type %r" % (kern_type,))

    visited = (bin_pos != 0).astype(float)
    # imfilter's default boundary is zero-padding; convolve with a symmetric
    # kernel reproduces it.
    denom = convolve(visited, kernel, mode="constant", cval=0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        smooth_pos = convolve(bin_pos, kernel, mode="constant", cval=0.0) / denom
        smooth_spk = convolve(bin_spk, kernel, mode="constant", cval=0.0) / denom
        rate = smooth_spk / smooth_pos

    unvisited = visited == 0
    rate[unvisited] = np.nan
    smooth_pos[unvisited] = np.nan
    smooth_spk[unvisited] = np.nan
    return rate, smooth_pos, smooth_spk


def _gaussian_kernel(shape, sigma):
    ky = np.arange(shape[0]) - (shape[0] - 1) / 2.0
    kx = np.arange(shape[1]) - (shape[1] - 1) / 2.0
    gy = np.exp(-(ky ** 2) / (2 * sigma ** 2))
    gx = np.exp(-(kx ** 2) / (2 * sigma ** 2))
    kernel = np.outer(gy, gx)
    return kernel / kernel.sum()


def rate_map(recording, electrode, cluster, bin_cm=RATEMAP_BIN_CM,
             smooth_bins=RATEMAP_SMOOTH_BINS):
    """Smoothed 2-D firing-rate map (Hz) for one cluster.

    Reproduces ``easy_make_ratemap.m``.

    Parameters
    ----------
    recording : Recording
    electrode : int
        Tetrode/group id.
    cluster : int
        Cluster id within that electrode.
    bin_cm : float
        Spatial bin size in cm.
    smooth_bins : int
        Boxcar smoothing width in bins.

    Returns
    -------
    ndarray (n_y_bins, n_x_bins) of firing rate in Hz, NaN where unvisited.
    """
    elec = recording.electrode(electrode)
    occ = occupancy_map(recording, bin_cm)
    spk = spike_map(recording, elec.spike_pos_index(cluster), bin_cm)
    rate, _, _ = smooth_ratemap(occ, spk, smooth_bins, "boxcar")
    return rate
