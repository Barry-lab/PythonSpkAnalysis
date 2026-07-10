"""Tests for LFP power spectrum, instantaneous freq/power, and ripple detection."""

import numpy as np
import pytest

from barrylab.analysis.lfp import detect_ripples, instant_freq_power, power_spectrum


def test_power_spectrum_recovers_known_frequency():
    fs = 250.0
    t = np.arange(0, 60, 1.0 / fs)
    rng = np.random.RandomState(0)
    signal = 1e-4 * np.sin(2 * np.pi * 8.0 * t) + 1e-6 * rng.randn(len(t))
    result = power_spectrum(signal, theta_range=(6.0, 12.0), eeg_freq=fs)
    assert abs(result["max_freq"] - 8.0) < 0.5
    assert result["s2n"] > 3.0  # a clean sinusoid should stand out sharply


def test_power_spectrum_shapes_and_truncation():
    fs = 250.0
    rng = np.random.RandomState(1)
    signal = rng.randn(int(fs * 10))
    result = power_spectrum(signal, upper_freq=25.0, eeg_freq=fs)
    assert np.all(result["freqs"] <= 25.0)
    assert result["power"].shape == result["freqs"].shape


def test_power_spectrum_pos2use_restricts_to_subset():
    fs = 250.0
    rng = np.random.RandomState(2)
    n_pos = 200
    signal = 1e-4 * np.sin(2 * np.pi * 8.0 * np.arange(0, n_pos * fs / 50) / fs)
    signal = signal + 1e-6 * rng.randn(len(signal))
    pos2use = np.arange(50, 150)
    result_subset = power_spectrum(signal, pos2use=pos2use, eeg_freq=fs)
    result_full = power_spectrum(signal, eeg_freq=fs)
    # both should find the same embedded frequency; subset just uses fewer samples
    assert abs(result_subset["max_freq"] - result_full["max_freq"]) < 1.0


def test_power_spectrum_on_real_lfp(example_trial):
    import barrylab as bl

    rec = bl.load_dacq(example_trial, load_lfp=[1])
    ch = rec.lfp[1]
    result = bl.power_spectrum(ch["samples"], eeg_freq=ch["sample_rate"])
    # a foraging rat's hippocampal LFP should show a clear, dominant theta peak
    assert 6.0 <= result["max_freq"] <= 11.0
    assert result["s2n"] > 2.0


# --- instant_freq_power ------------------------------------------------------


def test_instant_freq_power_recovers_known_frequency():
    fs = 250.0
    t = np.arange(0, 60, 1.0 / fs)
    rng = np.random.RandomState(0)
    signal = 1e-4 * np.sin(2 * np.pi * 8.0 * t) + 1e-6 * rng.randn(len(t))
    result = instant_freq_power(signal, (6.0, 12.0), eeg_freq=fs, pos_sample_rate=50.0)
    # away from the very first/last samples the recovered instantaneous
    # frequency should sit tightly around the known 8 Hz carrier
    assert abs(np.nanmean(result["freq"]) - 8.0) < 0.1
    assert abs(np.nanmean(result["freq_pos_rate"]) - 8.0) < 0.1
    assert np.isnan(result["freq"][0]) and np.isnan(result["freq"][-1])
    assert result["amp"].shape == result["freq"].shape == signal.shape
    assert result["filt_eeg"].shape == signal.shape


def test_instant_freq_power_pos_rate_none_skips_downsampling():
    fs = 250.0
    rng = np.random.RandomState(3)
    signal = rng.randn(int(fs * 5))
    result = instant_freq_power(signal, (6.0, 12.0), eeg_freq=fs, pos_sample_rate=None)
    assert result["freq_pos_rate"] is None
    assert result["amp_pos_rate"] is None


def test_instant_freq_power_downsample_requires_integer_ratio():
    fs = 250.0
    rng = np.random.RandomState(4)
    signal = rng.randn(int(fs * 5))
    with pytest.raises(ValueError):
        instant_freq_power(signal, (6.0, 12.0), eeg_freq=fs, pos_sample_rate=37.0)


def test_instant_freq_power_drops_trailing_partial_block_without_crashing():
    # MATLAB's reshape() would error unless len(eeg) is an exact multiple of
    # eeg_freq/pos_sample_rate; this port instead drops the trailing partial
    # block (a deliberate design decision -- see FINDINGS.md).
    fs = 250.0
    rng = np.random.RandomState(5)
    n = int(fs * 5) + 3  # 3 samples beyond an exact multiple of 5 (250/50)
    signal = rng.randn(n)
    result = instant_freq_power(signal, (6.0, 12.0), eeg_freq=fs, pos_sample_rate=50.0)
    assert result["freq_pos_rate"].shape[0] == n // 5


def test_instant_freq_power_on_real_lfp(example_trial):
    import barrylab as bl

    rec = bl.load_dacq(example_trial, load_lfp=[1])
    ch = rec.lfp[1]
    fs = ch["sample_rate"]
    # first 60s only -- keeps the (deliberately 1-second-order) FIR filter
    # tractable while still being a meaningful real-data check.
    chunk = ch["samples"][: int(fs * 60)]
    result = instant_freq_power(chunk, (6.0, 12.0), eeg_freq=fs, pos_sample_rate=None)
    # a foraging rat's hippocampal LFP should show textbook theta-band
    # instantaneous frequency -- the same real-data sanity check already
    # used for power_spectrum, replicated here (see FINDINGS.md).
    assert 6.0 <= np.nanmean(result["freq"]) <= 11.0
    assert np.mean(np.isnan(result["freq"])) < 0.05


# --- detect_ripples -----------------------------------------------------------


def _inject_bursts(
    fs,
    duration,
    burst_times,
    burst_dur=0.15,
    burst_freq=200.0,
    burst_amp=5e-5,
    noise_amp=1e-6,
    seed=1,
):
    n = int(fs * duration)
    t = np.arange(n) / fs
    rng = np.random.RandomState(seed)
    eeg = noise_amp * rng.randn(n)
    for bt in burst_times:
        idx = (t >= bt) & (t < bt + burst_dur)
        eeg[idx] += burst_amp * np.sin(2 * np.pi * burst_freq * (t[idx] - bt))
    return eeg


def test_detect_ripples_recovers_injected_bursts():
    fs = 1200.0
    burst_times = [3.0, 8.5, 15.2]
    eeg = _inject_bursts(fs, 20.0, burst_times)
    result = detect_ripples(eeg, fs, rip_freq=(150.0, 250.0), rip_thresh=5.0)

    assert result["rip_ind"].shape[0] == len(burst_times)
    detected_starts = np.sort(result["rip_ind"][:, 0] / fs)
    for expected, detected in zip(sorted(burst_times), detected_starts):
        # the low-threshold expansion pads the detected window a bit either
        # side of the injected burst, but the start should still land close by
        assert abs(expected - detected) < 0.1
    # duration must exceed the 60ms minimum and the recovered ripple-band
    # frequency should sit close to the injected 200 Hz oscillation
    assert np.all(result["duration"] > 0.06)
    assert np.all(np.abs(result["mean_freq"] - 200.0) < 15.0)
    assert result["power"].shape == eeg.shape


def test_detect_ripples_no_activity_returns_empty():
    fs = 1200.0
    rng = np.random.RandomState(6)
    eeg = 1e-6 * rng.randn(int(fs * 10))  # pure background noise, no bursts
    result = detect_ripples(eeg, fs)
    # exercises the empty-candidate-list path that crashes in the MATLAB
    # (ripStr(1)>ripEnd(1) indexing an empty array) -- see FINDINGS.md
    assert result["rip_ind"].shape == (0, 2)
    assert result["pk_power"].shape == (0,)
    assert result["duration"].shape == (0,)


def test_detect_ripples_default_threshold_matches_documented_code_default():
    # The MATLAB source's own docstring claims a default of 4, but its actual
    # code default (dfRipThresh) is 5 -- confirmed by reading the source (see
    # FINDINGS.md). This pins the port's default to that real behaviour.
    fs = 1200.0
    eeg = _inject_bursts(fs, 20.0, [5.0])
    default_result = detect_ripples(eeg, fs)
    explicit_5 = detect_ripples(eeg, fs, rip_thresh=5.0)
    explicit_4 = detect_ripples(eeg, fs, rip_thresh=4.0)
    np.testing.assert_array_equal(default_result["rip_ind"], explicit_5["rip_ind"])
    # a looser (4 SD) threshold must never detect *fewer* ripples than a
    # stricter (5 SD) one
    assert explicit_4["rip_ind"].shape[0] >= explicit_5["rip_ind"].shape[0]


def test_detect_ripples_on_real_lfp(example_trial):
    import barrylab as bl

    rec = bl.load_dacq(example_trial, load_lfp=[1])
    ch = rec.lfp[1]
    fs = ch["sample_rate"]
    chunk = ch["samples"][: int(fs * 120)]
    result = detect_ripples(chunk, fs)
    # plausibility, not a precise oracle: any detected events should have a
    # sensible duration and a mean frequency that actually sits inside the
    # requested ripple band, not clustering at its edges (which would
    # suggest a filter-design or indexing bug rather than real events)
    assert np.all(result["duration"] > 0.06)
    if result["mean_freq"].size:
        assert np.all(result["mean_freq"] >= 150.0) and np.all(result["mean_freq"] <= 250.0)
