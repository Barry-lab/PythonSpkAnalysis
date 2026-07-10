"""Head-direction analysis: tuning curves and directional metrics.

Ports the lab's ``Universal_Matlab/Analysis_Dir`` canon:

- :func:`head_direction_map` — firing rate as a function of head direction
  (the directional analogue of a rate map: bin by direction, smooth circularly,
  divide), built with ``bin_data('direction')`` + ``make_smooth_ratemap`` logic.
- :func:`circular_mean` — weighted circular mean (``cb_circMean.m``).
- :func:`resultant_vector` — mean resultant length and preferred direction, the
  standard Rayleigh measure of tuning strength.
- :func:`directionality` — KL divergence of the tuning curve from uniform, the
  lab's directionality index (``dir_kldiv.m`` / ``kldiv.m``).
- :func:`directional_fft` — FFT of the tuning curve to detect multimodal (e.g.
  hexadirectional) tuning (``circ_FTT.m``).
"""

import numpy as np
from scipy.ndimage import gaussian_filter1d

from ..params import DIR_BIN_DEG, DIR_SMOOTH_BINS
from .maps import smooth_ratemap

_EPS = np.finfo(float).eps


def _histnd_1d(values, range_max, bin_size):
    """1-D histogram matching ``histnd.m`` (bins are ``ceil(value / bin_size)``)."""
    v = np.asarray(values, dtype=float)
    q = v / bin_size
    q = np.where(q == 0, q + _EPS, q)
    bins = np.ceil(q)
    max_bin = int(np.ceil(range_max / bin_size))
    valid = ~np.isnan(bins)
    bins = bins[valid].astype(np.int64)
    bins[bins < 1] = 1
    bins = np.minimum(bins, max_bin)
    out = np.zeros(max_bin)
    np.add.at(out, bins - 1, 1.0)
    return out


def directional_occupancy(recording, bin_deg=DIR_BIN_DEG):
    """Dwell time (seconds) per head-direction bin."""
    counts = _histnd_1d(recording.direction, 360.0, bin_deg)
    return counts / recording.pos_sample_rate


def head_direction_map(recording, electrode, cluster, bin_deg=DIR_BIN_DEG,
                       smoothing="boxcar", smooth_size=DIR_SMOOTH_BINS):
    """Firing rate (Hz) as a function of head direction, for one cluster.

    The directional analogue of :func:`~barrelab.analysis.maps.rate_map`:
    occupancy and spikes are binned over direction and smoothed **circularly**
    (the ends wrap) before dividing.

    Parameters
    ----------
    recording : Recording
    electrode, cluster : int
    bin_deg : float
        Directional bin width in degrees (default 6 -> 60 bins).
    smoothing : {"boxcar", "gaussian", None}
    smooth_size : float
        Boxcar width or gaussian sigma in bins.

    Returns
    -------
    ndarray
        Tuning curve of length ``360 / bin_deg``; bin ``i`` covers directions
        ``[i*bin_deg, (i+1)*bin_deg)`` degrees.
    """
    elec = recording.electrode(electrode)
    occ = directional_occupancy(recording, bin_deg)
    spk = _histnd_1d(recording.direction[elec.spike_pos_index(cluster)], 360.0, bin_deg)
    if smoothing is None:
        with np.errstate(divide="ignore", invalid="ignore"):
            rate = spk / occ
        rate[occ == 0] = np.nan
        return rate
    rate, _, _ = smooth_ratemap(occ, spk, smooth_size, smoothing, smooth_type="circular")
    return rate


def circular_mean(angles, weights=None):
    """Weighted circular mean of angles in radians. Port of ``cb_circMean.m``.

    Returns the mean angle in ``[0, 2*pi)``. NaNs in either input are treated as
    mutual and ignored.
    """
    angles = np.asarray(angles, dtype=float)
    w = np.ones_like(angles) if weights is None else np.asarray(weights, dtype=float)
    w = np.where(np.isnan(angles), np.nan, w)
    angles = np.where(np.isnan(w), np.nan, angles)
    sum_w = np.nansum(w)
    mean = np.arctan2(np.nansum(np.sin(angles) * w) / sum_w,
                      np.nansum(np.cos(angles) * w) / sum_w)
    return np.mod(mean, 2 * np.pi)


def resultant_vector(tuning_curve):
    """Mean resultant length and preferred direction of a tuning curve.

    The standard Rayleigh measure of directional tuning: a value near 0 means
    uniform (untuned) firing, near 1 means sharply tuned.

    Returns
    -------
    (r, preferred_deg)
        ``r`` in [0, 1]; ``preferred_deg`` in [0, 360) (``nan`` if the cell is
        silent).
    """
    rates = np.nan_to_num(np.asarray(tuning_curve, dtype=float))
    n = len(rates)
    theta = (np.arange(n) + 0.5) * (2 * np.pi / n)  # bin-centre angles
    total = rates.sum()
    if total == 0:
        return 0.0, np.nan
    x = np.sum(rates * np.cos(theta)) / total
    y = np.sum(rates * np.sin(theta)) / total
    return float(np.hypot(x, y)), float(np.mod(np.degrees(np.arctan2(y, x)), 360.0))


def _kldiv(p1, p2, mode=None):
    """Kullback-Leibler / Jensen-Shannon divergence (bits). Port of ``kldiv.m``."""
    p1 = np.asarray(p1, dtype=float)
    p2 = np.asarray(p2, dtype=float)
    if mode == "js":
        q = (p1 + p2) / 2
        return 0.5 * (np.sum(p1 * (np.log2(p1) - np.log2(q)))
                      + np.sum(p2 * (np.log2(p2) - np.log2(q))))
    if mode == "sym":
        kl1 = np.sum(p1 * (np.log2(p1) - np.log2(p2)))
        kl2 = np.sum(p2 * (np.log2(p2) - np.log2(p1)))
        return (kl1 + kl2) / 2
    return float(np.sum(p1 * (np.log2(p1) - np.log2(p2))))


def directionality(tuning_curve):
    """KL divergence of a tuning curve from a uniform circle. Port of ``dir_kldiv.m``.

    A directionality index: how far the (smoothed) directional firing departs
    from uniform. Tom Wills found it correlates very highly with Skaggs
    bits/spike; >0.15-0.25 is a common "directional" cut-off. Should be run on a
    smoothed tuning curve.
    """
    rm = np.asarray(tuning_curve, dtype=float).copy()
    # KL is undefined at 0; nudge zeros (and any unvisited NaN bins) up slightly.
    rm[~np.isfinite(rm)] = 0.0
    rm[rm == 0] = 1e-5
    p = rm / rm.sum()
    n = len(rm)
    uniform = np.ones(n) / n
    return _kldiv(p, uniform)


def directional_fft(tuning_curve, smooth_hz=None):
    """FFT of a directional tuning curve. Port of ``circ_FTT.m``.

    Detects multimodal directional tuning (e.g. the 6-fold "hexadirectional"
    signal). The tuning curve should be **unsmoothed**.

    Returns
    -------
    (power, freqs, rel_freq, rel_power)
        ``power`` at each frequency component (mean removed, so DC ~ 0);
        ``freqs`` the matching cycles-per-circle; ``rel_freq = [2,3,4,5,7]`` and
        ``rel_power`` their power relative to the 6-cycle component.
    """
    dr = np.asarray(tuning_curve, dtype=float).ravel()
    samp_freq = len(dr)
    dr = dr - dr.mean()

    fft_len = 360  # pad so a clean 6-cycle component is resolvable
    fftd = np.fft.fft(dr, fft_len)
    power = (fftd * np.conj(fftd)).real / fft_len
    power = power[: fft_len // 2 + 1]
    freqs = (samp_freq / 2) * np.linspace(0, 1, fft_len // 2 + 1)

    if smooth_hz:
        kern_sig = round(smooth_hz * ((fft_len // 2 + 1) / (samp_freq / 2)))
        if kern_sig > 0:
            power = gaussian_filter1d(power, kern_sig, mode="reflect")

    rel_freq = np.array([2, 3, 4, 5, 7])
    rel_power = (power / power[6])[rel_freq]  # 6-cycle component is power[6]
    return power, freqs, rel_freq, rel_power


def head_direction_stats(recording, electrode, cluster, bin_deg=DIR_BIN_DEG):
    """Convenience: tuning curve + summary head-direction metrics for one cluster.

    Returns
    -------
    dict with ``tuning_curve``, ``resultant_length``, ``preferred_direction``
    (deg) and ``directionality`` (KL divergence).
    """
    tc = head_direction_map(recording, electrode, cluster, bin_deg=bin_deg)
    r, preferred = resultant_vector(tc)
    return {
        "tuning_curve": tc,
        "resultant_length": r,
        "preferred_direction": preferred,
        "directionality": directionality(tc),
    }
