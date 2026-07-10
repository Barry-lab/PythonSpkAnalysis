"""Tests for Skaggs spatial information (hand-computed reference values)."""

import numpy as np

from barrelab.analysis.spatial import spatial_information


def test_spatial_information_hand_example():
    # Two equally-occupied bins; one silent, one at 10 Hz.
    #   mean_rate = 5;  p_x = [0.5, 0.5];  p_r = [0, 2]
    #   bits/s   = 0.5 * 10 * log2(2) = 5
    #   bits/spk = 5 / 5 = 1
    bits_per_spike, bits_per_sec = spatial_information([0.0, 10.0], [0.5, 0.5])
    assert np.isclose(bits_per_spike, 1.0)
    assert np.isclose(bits_per_sec, 5.0)


def test_spatial_information_silent_cell_is_nan():
    bps, bpsec = spatial_information([0.0, 0.0], [0.5, 0.5])
    assert np.isnan(bps) and np.isnan(bpsec)


def test_spatial_information_uniform_firing_is_zero():
    # A cell that fires identically everywhere carries no spatial information.
    rates = np.full((4, 4), 3.0)
    dwell = np.full((4, 4), 0.1)
    bps, _ = spatial_information(rates, dwell)
    assert np.isclose(bps, 0.0)


def test_spatial_information_ignores_nan_bins():
    a = spatial_information([np.nan, 0.0, 10.0], [np.nan, 0.5, 0.5])
    b = spatial_information([0.0, 10.0], [0.5, 0.5])
    assert np.allclose(a, b)
