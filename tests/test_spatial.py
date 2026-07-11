"""Tests for Skaggs spatial information (hand-computed reference values)."""

import numpy as np

from barrylab.analysis.spatial import adaptive_binning, spatial_information


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


# --------------------------------------------------------------------------
# adaptive_binning (skaggsAdaptiveBinning.m) -- previously untested, which is
# how the empty-spike-map return bug (2nd output was the dwell map, not the
# spike map) slipped through.
# --------------------------------------------------------------------------
def test_adaptive_binning_returns_rate_spk_pos_order():
    rng = np.random.RandomState(0)
    dwell = rng.uniform(0.2, 1.0, size=(12, 12))
    spikes = rng.poisson(2.0, size=(12, 12)).astype(float)
    rate, smooth_spk, smooth_pos = adaptive_binning(dwell, spikes)
    assert rate.shape == smooth_spk.shape == smooth_pos.shape == dwell.shape
    # rate should equal smooth_spk / smooth_pos where finite
    finite = np.isfinite(rate) & np.isfinite(smooth_spk) & np.isfinite(smooth_pos)
    assert np.allclose(rate[finite], smooth_spk[finite] / smooth_pos[finite])
    # smoothed dwell should be positive and on the order of the raw dwell,
    # NOT accidentally the spike counts
    assert np.all(smooth_pos[finite] > 0)


def test_adaptive_binning_empty_spike_map():
    # Regression: with no spikes, the 2nd return (smooth_spk) must be the
    # all-zero spike map, not the dwell map (the fixed bug).
    dwell = np.ones((6, 6)) * 0.5
    dwell[0, 0] = 0.0  # one unvisited bin
    spikes = np.zeros((6, 6))
    rate, smooth_spk, smooth_pos = adaptive_binning(dwell, spikes)
    assert np.isnan(rate[0, 0]) and np.isnan(smooth_spk[0, 0]) and np.isnan(smooth_pos[0, 0])
    # visited bins: rate and smooth_spk are 0 (no spikes); smooth_pos is dwell
    visited = dwell > 0
    assert np.allclose(rate[visited], 0.0)
    assert np.allclose(smooth_spk[visited], 0.0)  # NOT the dwell value (0.5)
    assert np.allclose(smooth_pos[visited], 0.5)


def test_adaptive_binning_unvisited_bins_are_nan():
    dwell = np.ones((8, 8))
    dwell[2:4, 2:4] = 0.0  # an unvisited patch
    spikes = np.ones((8, 8))
    rate, _, _ = adaptive_binning(dwell, spikes)
    assert np.all(np.isnan(rate[2:4, 2:4]))
    assert np.all(np.isfinite(rate[dwell > 0]))
