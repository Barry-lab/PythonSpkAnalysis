"""LFP power spectrum, instantaneous phase/frequency, and ripple detection.

Ports ``eeg_power_spectra.m`` (FFT-based power spectrum), ``eeg_instant_freq_power.m``
(bandpass + Hilbert instantaneous phase/amplitude/frequency), and
``detect_ripples.m`` (Csicsvari-style high-frequency ripple/SWR detection).
Requires :class:`~barrylab.core.recording.Recording` loaded with
``load_dacq(..., load_lfp=True)``.
"""

import warnings

import numpy as np
from scipy.ndimage import convolve1d
from scipy.signal import detrend, filtfilt, firwin, hilbert


def _gaussian_kernel_1d(width_hz, sigma_hz, freq_res_hz):
    """MATLAB ``fspecial('gaussian', [1, kernelLength], sigma_bins)`` where
    ``kernelLength = count(freqs <= width_hz)`` and ``sigma_bins =
    kernelLength * (sigma_hz/width_hz)``."""
    kernel_length = max(1, int(round(width_hz / freq_res_hz)) + 1)
    sigma_bins = kernel_length * (sigma_hz / width_hz)
    x = np.arange(kernel_length) - (kernel_length - 1) / 2.0
    h = np.exp(-(x**2) / (2 * sigma_bins**2))
    h[h < np.finfo(float).eps * h.max()] = 0.0
    return h / h.sum()


def power_spectrum(
    eeg,
    pos2use=None,
    smooth_kern=(2.0, 0.5),
    theta_range=(6.0, 12.0),
    upper_freq=25.0,
    s2n_width=1.0,
    eeg_freq=250.0,
):
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

    return {
        "power": power,
        "freqs": freqs,
        "max_power": max_power,
        "max_freq": max_freq,
        "s2n": s2n,
    }


def instant_freq_power(eeg, filter_range, eeg_freq=250.0, pos_sample_rate=50.0):
    """Bandpass-filtered instantaneous phase/amplitude/frequency of an LFP trace.

    Port of ``eeg_instant_freq_power.m``: FIR-bandpass ``eeg`` in a chosen band
    (e.g. theta, ``(6, 12)`` Hz) using a Blackman-windowed FIR filter design
    matching MATLAB's ``fir1(eegFreq, [lo hi]/nyquist, blackman(eegFreq+1))``
    (a 1-second-duration filter, i.e. its order always equals the sample rate,
    however large -- this is what the original does, not a port simplification),
    zero-phase filter with ``filtfilt``, then take the analytic signal (Hilbert
    transform) to recover instantaneous phase, amplitude, and frequency.

    Parameters
    ----------
    eeg : ndarray
        Raw LFP samples, in volts for meaningful amplitude/power (a 2-column
        ``[samples, timestamp]`` array is also accepted -- only the first
        column is used, matching the MATLAB source's defensive check).
    filter_range : (float, float)
        Bandpass edges in Hz, e.g. ``(6, 12)`` for theta or ``(150, 250)`` for
        ripples.
    eeg_freq : float
        Sampling rate of ``eeg`` in Hz. Defaults to 250 (standard `.eeg`);
        pass the true rate (e.g. ``recording.lfp[ch]["sample_rate"]``, ~4800
        for `.egf`) otherwise. Must be (very close to) an integer -- both
        the FIR filter order and the Blackman window length are literally
        ``round(eeg_freq)``, matching MATLAB's ``fir1(eegSampFreq, ...)``.
    pos_sample_rate : float or None
        Rate (Hz) to downsample ``freq``/``amp`` to, matching position
        samples (50 Hz for DACQ `.pos`; pass ``recording.pos_sample_rate``).
        Pass ``None`` to skip the downsampling step entirely (e.g. when only
        the full-rate outputs are needed, as in :func:`detect_ripples`).

    Returns
    -------
    dict with:
        ``filt_eeg`` -- bandpass-filtered eeg (full ``eeg_freq`` rate).
        ``phase`` -- unwrapped instantaneous phase (radians), monotonically
            increasing (not wrapped back to 0 at 2*pi); 0 phase is defined at
            an EEG peak (a property of the Hilbert-transform analytic signal
            convention, not extra logic). NaN wherever amplitude is below a
            fixed ``1e-7`` floor (the original's fix for spurious phase
            estimates at near-zero amplitude, which otherwise cluster
            suspiciously around ~4.2 rad rather than being uniformly random).
        ``amp`` -- instantaneous amplitude (this is the analytic signal's
            modulus -- true *amplitude*, not squared power, matching the
            MATLAB source's own comment). NaN under the same amplitude floor.
        ``freq`` -- instantaneous frequency (Hz), full ``eeg_freq`` rate. NaN
            at the first/last sample (undefined -- see note below) and
            wherever the Hilbert transform produces a spurious negative
            frequency (regressing phase; the MATLAB source treats these, and
            their matching phase/amp samples, as unreliable).
        ``freq_pos_rate`` / ``amp_pos_rate`` -- ``freq``/``amp`` downsampled
            (``nanmean`` over non-overlapping blocks of size
            ``eeg_freq/pos_sample_rate``) to ``pos_sample_rate``. ``None`` if
            ``pos_sample_rate`` is ``None``.

    Notes
    -----
    **Instantaneous-frequency centring, ported without the MATLAB indirection.**
    MATLAB computes ``diff(phase)`` (length N-1, each element the frequency
    "between" two adjacent samples) then re-centres it back onto all N
    original sample points via a two-argument ``interp2`` call querying
    half-integer positions ``0.5, 1.5, ..., N+0.5`` against an
    ``[freqEEG, ones(...)]`` matrix -- an awkward way of writing "linearly
    interpolate the diff array onto the midpoints between its own samples,
    and leave the two out-of-range endpoints as NaN". Traced through by hand
    and reproduced directly (equivalent result, no 2-D interpolation
    indirection needed): ``freq[0] = freq[-1] = nan``, and
    ``freq[i] = mean(d[i-1], d[i])`` for interior ``i``, where ``d`` is the
    raw phase-diff array.

    **``eeg_freq`` must yield an exact number of downsampling blocks.**
    MATLAB's ``reshape(freqEEG, downSampConst, [])`` requires the array
    length to be an *exact* multiple of ``downSampConst`` and errors
    otherwise. Real recordings occasionally have a handful of trailing
    samples beyond the last full block (rounding during acquisition), so
    rather than replicate that crash, trailing leftover samples here are
    simply dropped before reshaping -- a deliberate design decision (there
    is no well-defined "correct" MATLAB behaviour to preserve for this case,
    since MATLAB itself would just error); see ``FINDINGS.md``.
    """
    eeg = np.asarray(eeg, dtype=float)
    if eeg.ndim > 1:
        # Defensive column-orientation + drop-timestamp-column check,
        # mirroring the MATLAB source exactly (eeg_convert_to_volts already
        # drops any timestamp column upstream, but the original re-checks).
        if eeg.shape[1] > eeg.shape[0]:
            eeg = eeg.T
        eeg = eeg[:, 0]
    eeg = eeg.ravel()
    n_samples = eeg.shape[0]

    # Detrend (remove best-fit line / DC shift), same convention as power_spectrum.
    eeg = detrend(eeg)

    # --- Bandpass filter design + zero-phase filtering ----------------------
    # fir1(eegSampFreq, [lo hi]/nyquist, blackman(eegSampFreq+1)) -> a
    # Blackman-windowed FIR bandpass of order eegSampFreq (eegSampFreq+1
    # taps). scipy.signal.firwin's default frequency-normalisation (fs=) and
    # scaling (scale=True, matching fir1's default) reproduce this design.
    numtaps = int(round(eeg_freq)) + 1
    fir_coeffs = firwin(numtaps, filter_range, window="blackman", pass_zero=False, fs=eeg_freq)
    filt_eeg = filtfilt(fir_coeffs, [1.0], eeg)
    analytic = hilbert(filt_eeg)

    phase = np.unwrap(np.angle(analytic))
    amp = np.abs(analytic)

    # Low-amplitude points have spurious phase/frequency -> nan them out.
    amp_filt = 1e-7
    low_amp = amp < amp_filt
    perc_low_amp = 100.0 * np.count_nonzero(low_amp) / len(amp)
    if perc_low_amp > 10.0:
        warnings.warn(
            f">10% ({perc_low_amp:.2g}%) of EEG points have v. low amplitude - set to nan.",
            stacklevel=2,
        )
    phase = np.where(low_amp, np.nan, phase)
    amp = np.where(low_amp, np.nan, amp)

    # Instantaneous frequency (Hz) -- see docstring note on the interp2
    # indirection this replaces.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)  # nan arithmetic
        d_phase = np.diff(phase) * (eeg_freq / (2.0 * np.pi))
    freq = np.full(n_samples, np.nan)
    if n_samples > 2:
        freq[1:-1] = (d_phase[:-1] + d_phase[1:]) / 2.0

    # Spurious (Hilbert-artefact) negative frequencies -> nan, and null the
    # matching phase/amp too. NaN comparisons are False in both MATLAB and
    # numpy, so the pre-existing nan positions are left untouched by this.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        neg = freq < 0
    phase = np.where(neg, np.nan, phase)
    amp = np.where(neg, np.nan, amp)
    freq = np.where(neg, np.nan, freq)

    if pos_sample_rate is None:
        freq_pos_rate = None
        amp_pos_rate = None
    else:
        down_samp_const = eeg_freq / pos_sample_rate
        if not np.isclose(down_samp_const, round(down_samp_const)):
            raise ValueError(
                f"eeg_freq/pos_sample_rate ({down_samp_const:.6f}) is not an integer -- cannot "
                "downsample to whole blocks."
            )
        down_samp_const = int(round(down_samp_const))
        n_blocks = n_samples // down_samp_const
        n_used = n_blocks * down_samp_const  # drop any trailing partial block
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)  # all-nan blocks
            freq_pos_rate = np.nanmean(freq[:n_used].reshape(n_blocks, down_samp_const), axis=1)
            amp_pos_rate = np.nanmean(amp[:n_used].reshape(n_blocks, down_samp_const), axis=1)

    return {
        "filt_eeg": filt_eeg,
        "phase": phase,
        "amp": amp,
        "freq": freq,
        "freq_pos_rate": freq_pos_rate,
        "amp_pos_rate": amp_pos_rate,
    }


def detect_ripples(eeg, eeg_freq, rip_freq=(150.0, 250.0), rip_thresh=5.0):
    """Detect ripple (sharp-wave ripple / SWR) events in a high-frequency LFP trace.

    Port of ``detect_ripples.m`` (Csicsvari-style ripple detection): finds
    contiguous regions where smoothed ripple-band power exceeds
    ``mean + rip_thresh*std``, expands each out to the surrounding region
    above the looser ``mean + 0.5*std``, and keeps the resulting events that
    are both long enough and actually contain some genuinely high-power
    activity.

    Requires a high sample rate (e.g. the `.egf`, ~4.8 kHz channel) -- per
    the MATLAB source's own warning, detecting a 150-250 Hz band needs eeg
    sampled at >=500 Hz to avoid Nyquist problems.

    Parameters
    ----------
    eeg : ndarray
        Unfiltered LFP (any scale -- detection is relative to its own mean
        and std, so it doesn't need volts conversion).
    eeg_freq : float
        Sample rate of ``eeg`` in Hz (see :func:`instant_freq_power` for the
        constraint that this must be very close to an integer).
    rip_freq : (float, float)
        Ripple frequency band, Hz. Default ``(150, 250)``.
    rip_thresh : float
        Detection threshold, in SD above the mean ripple-band power.
        **NB: the MATLAB source's own header comment claims a default of 4
        ("Csisvari uses 4"), but the code's actual default (``dfRipThresh``)
        is 5 -- a genuine doc/code mismatch in the original, confirmed by
        reading the source. This port's default (5) matches the code that
        actually runs, not the stale comment -- see ``FINDINGS.md``.**

    Returns
    -------
    dict with:
        ``rip_ind`` -- ``(n_ripples, 2)`` 0-based ``[start, end]`` sample
            index (into ``eeg``) of each ripple, both inclusive.
        ``power`` -- ``(len(eeg),)`` smoothed ripple-band power (squared
            amplitude, 0.05s boxcar-smoothed) for the whole trace.
        ``pk_power`` -- ``(n_ripples,)`` peak power within each ripple.
        ``mean_freq`` -- ``(n_ripples,)`` mean instantaneous frequency (Hz)
            of the filtered eeg over each ripple's extent.
        ``duration`` -- ``(n_ripples,)`` duration of each ripple, in seconds.

    Notes
    -----
    **A ripple still in progress at the very end of ``eeg`` is silently
    dropped.** If ripple-band power crosses back below the low threshold for
    the last time and then rises again but never falls before the trace
    ends, the MATLAB source's ``ripStr = ripStr(1:length(ripEnd))`` keeps
    only the *first* ``length(ripEnd)`` start indices, silently discarding
    this final, unterminated event rather than reporting it against the end
    of the recording. Confirmed by reading the source (not assumed) and
    faithfully reproduced here, since it is well-defined (if lossy)
    behaviour, not a crash -- see ``FINDINGS.md``.

    **An empty candidate-ripple list is handled explicitly, unlike the
    MATLAB.** ``ripStr(1)>ripEnd(1)`` indexes element 1 of what can be an
    *empty* array (e.g. no sample anywhere in the trace crosses the low
    threshold -- a perfectly plausible outcome, not just a theoretical edge
    case, for a short or genuinely ripple-free recording), which errors in
    MATLAB. "No ripples found" is a well-defined result worth returning
    rather than a crash worth reproducing, so this port returns empty
    arrays for that case instead.

    **The MATLAB's ``[valid, pkPower, meanFreq]=deal(zeros(length(ripStr)))``
    is a confirmed, genuine sloppy-coding bug (``zeros(n)`` with a scalar `n`
    makes an n-by-n matrix, not the evidently-intended n-by-1 column) -- but
    traced through by hand and confirmed numerically harmless.** The
    following loop only ever linear-indexes positions ``1..n``, which for an
    n-by-n matrix stored column-major fall entirely within column 1 -- i.e.
    exactly the same memory positions a proper n-by-1 vector would use -- so
    every value the loop writes, and everything later read back out via
    ``valid(n)``/logical-indexing, is bit-identical to what a correctly
    n-by-1-shaped implementation would give. The only real cost is a
    wasted, potentially large ``n^2`` allocation. This port simply uses
    correctly-shaped 1-D arrays throughout (there is no output-changing
    "bug" here to faithfully preserve, just a memory-wasteful implementation
    detail not worth reproducing) -- see ``FINDINGS.md`` for the full
    index-tracing argument.
    """
    rip_low_thresh = 0.5
    # 60ms. NB the MATLAB source's own inline comment claims "40ms = 0.004"
    # here, which is wrong on two independent counts (60ms != 40ms, and
    # 40ms = 0.04s, not 0.004s) -- a stale/miscopied comment, not a code bug;
    # the actual *value* used (0.06s = 60ms) is what's ported.
    rip_min_dur = 60.0 / 1000.0

    result = instant_freq_power(eeg, rip_freq, eeg_freq=eeg_freq, pos_sample_rate=None)
    amp = result["amp"]
    freq = result["freq"]

    power = amp**2
    power = np.where(np.isnan(power), 0.0, power)

    # Smooth with a ~half-ripple-width (0.05s) boxcar. filter2(kern, X) with a
    # symmetric all-ones (here normalised) column kernel in MATLAB's default
    # 'same'/zero-padded mode is exactly convolve1d(..., mode="constant",
    # cval=0.0) -- already established and validated in line_trajectory.py
    # (see FINDINGS.md), reused here rather than re-derived.
    kern_len = int(round(eeg_freq / 20.0))
    kern = np.ones(kern_len) / kern_len
    power = convolve1d(power, kern, mode="constant", cval=0.0)

    mean_power = np.nanmean(power)
    std_power = np.nanstd(power, ddof=1)  # MATLAB (nan)std default is N-1
    rips_full = power > (mean_power + rip_thresh * std_power)
    rips_low = power > (mean_power + rip_low_thresh * std_power)

    # find(diff(x)==1)/==-1 rising/falling edge indices. The index returned
    # is the sample *just before* the transition in both MATLAB (1-based) and
    # here (0-based) -- diff shifts alignment identically either way, so no
    # extra +/-1 correction is needed translating this specific idiom.
    edges = np.diff(rips_low.astype(np.int8))
    rip_str = np.flatnonzero(edges == 1)
    rip_end = np.flatnonzero(edges == -1)

    if rip_str.size == 0 or rip_end.size == 0:
        # See docstring: MATLAB crashes here on an empty array; "no
        # candidate ripples anywhere" is a well-defined result instead.
        return {
            "rip_ind": np.empty((0, 2), dtype=np.int64),
            "power": power,
            "pk_power": np.empty(0),
            "mean_freq": np.empty(0),
            "duration": np.empty(0),
        }

    if rip_str[0] > rip_end[0]:
        # Trace starts already above the low threshold (an orphan "end"
        # with no matching "start" inside the trace) -- drop it.
        rip_end = rip_end[1:]
    # Drop any trailing, unterminated "start" -- see docstring note.
    rip_str = rip_str[: len(rip_end)]

    duration = (rip_end - rip_str) / eeg_freq
    long_enough = duration > rip_min_dur
    rip_str = rip_str[long_enough]
    rip_end = rip_end[long_enough]

    n_candidates = len(rip_str)
    valid = np.zeros(n_candidates, dtype=bool)
    pk_power = np.zeros(n_candidates)
    mean_freq = np.zeros(n_candidates)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)  # all-nan freq slices
        for n in range(n_candidates):
            s, e = rip_str[n], rip_end[n] + 1  # +1: MATLAB's a:b range is inclusive
            valid[n] = np.any(rips_full[s:e])
            pk_power[n] = np.max(power[s:e])
            mean_freq[n] = np.nanmean(freq[s:e])

    rip_str = rip_str[valid]
    rip_end = rip_end[valid]
    pk_power = pk_power[valid]
    mean_freq = mean_freq[valid]
    duration = (rip_end - rip_str) / eeg_freq
    rip_ind = np.column_stack([rip_str, rip_end]).astype(np.int64)

    return {
        "rip_ind": rip_ind,
        "power": power,
        "pk_power": pk_power,
        "mean_freq": mean_freq,
        "duration": duration,
    }
