"""Tests for LFP power spectrum analysis."""

import numpy as np

from barrylab.analysis.lfp import power_spectrum


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
