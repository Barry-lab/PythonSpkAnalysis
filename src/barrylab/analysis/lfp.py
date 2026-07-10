"""LFP power spectrum analysis.

Ports ``eeg_power_spectra.m``: an FFT-based power spectrum of the local field
potential, restrictable to specific position samples, with theta-band peak
detection and signal-to-noise. Requires :class:`~barrylab.core.recording.Recording`
loaded with ``load_dacq(..., load_lfp=True)``.
"""

import numpy as np
from scipy.ndimage import convolve1d
from scipy.signal import detrend


def _gaussian_kernel_1d(width_hz, sigma_hz, freq_res_hz):
    """MATLAB ``fspecial('gaussian', [1, kernelLength], sigma_bins)`` where
    ``kernelLength = count(freqs <= width_hz)`` and ``sigma_bins =
    kernelLength * (sigma_hz/width_hz)``."""
    kernel_length = max(1, int(round(width_hz / freq_res_hz)) + 1)
    sigma_bins = kernel_length * (sigma_hz / width_hz)
    x = np.arange(kernel_length) - (kernel_length - 1) / 2.0
    h = np.exp(-(x ** 2) / (2 * sigma_bins ** 2))
    h[h < np.finfo(float).eps * h.max()] = 0.0
    return h / h.sum()


def power_spectrum(eeg, pos2use=None, smooth_kern=(2.0, 0.5), theta_range=(6.0, 12.0),
                   upper_freq=25.0, s2n_width=1.0, eeg_freq=250.0):
    """FFT power spectrum of an LFP trace, with theta peak and signal-to-noise.

    Port of ``eeg_power_spectra.m``.

    Parameters
    ----------
    eeg : ndarray
        LFP samples in volts (e.g. ``recording.lfp[ch]["samples"]``).
    pos2use : ndarray, optional
        0-based position-sample indices to restrict the analysis to (e.g. only
        high-speed running periods). ``None`` uses the whole trace.
    smooth_kern : (float, float)
        ``(width_hz, sigma_hz)`` for Gaussian-smoothing the power spectrum.
    theta_range : (float, float)
        Frequency band (Hz) to search for the peak.
    upper_freq : float
        Truncate the returned spectrum above this frequency (Hz).
    s2n_width : float
        Width (Hz) of the band around the theta peak used for signal-to-noise
        (half on each side).
    eeg_freq : float
        Sampling rate of ``eeg``, in Hz (250 for `.eeg`, ~4800 for `.egf`).

    Returns
    -------
    dict with ``power``, ``freqs`` (both truncated at ``upper_freq``),
    ``max_power``, ``max_freq`` (peak within ``theta_range``), ``s2n``.
    """
    eeg = np.asarray(eeg, dtype=float).ravel()
    eeg = detrend(eeg)

    if pos2use is not None:
        pos2use = np.asarray(pos2use, dtype=np.int64).ravel()
        eeg_v_pos_rate = int(round(eeg_freq / 50.0))
        offsets = np.arange(eeg_v_pos_rate)
        eeg_idx = np.sort((eeg_v_pos_rate * pos2use[:, None] + offsets[None, :]).ravel())
        eeg = eeg[eeg_idx]

    fft_length = 1 << (len(eeg) - 1).bit_length()  # next power of 2 >= len(eeg)
    fft_result = np.fft.fft(eeg, fft_length)
    freqs = np.arange(len(fft_result)) * (eeg_freq / len(fft_result))
    freqs = freqs[freqs <= upper_freq]

    power = (fft_result * np.conj(fft_result)).real / len(fft_result)
    power = power[: len(freqs)]

    freq_res = eeg_freq / fft_length
    kernel = _gaussian_kernel_1d(smooth_kern[0], smooth_kern[1], freq_res)
    power = convolve1d(power, kernel, mode="reflect")

    theta_mask = (freqs > theta_range[0]) & (freqs < theta_range[1])
    max_power = float(np.max(power[theta_mask]))
    max_freq = float(np.mean(freqs[theta_mask][power[theta_mask] == max_power]))

    peak_band = (freqs > max_freq - s2n_width / 2.0) & (freqs < max_freq + s2n_width / 2.0)
    min_freq_band = freqs >= 2.0  # lowerFrq, hardcoded in the original
    s2n = float(np.nanmean(power[peak_band]) / np.nanmean(power[~peak_band & min_freq_band]))

    return {"power": power, "freqs": freqs, "max_power": max_power,
            "max_freq": max_freq, "s2n": s2n}
