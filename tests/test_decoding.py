"""Tests for Bayesian population decoding."""

import numpy as np
from scipy.stats import poisson

from barrylab.analysis.decoding import (
    bayesian_posterior,
    bayesian_posterior_with_prior,
    decode_maximum_likelihood,
    spike_counts_per_window,
)


def test_bayesian_posterior_matches_poisson_pmf():
    # 1 cell, 2 spatial bins, rates 2Hz and 5Hz, t_bin=1s, 3 observed spikes
    rm = np.array([[[2.0], [5.0]]])  # (ny=1, nx=2, n_cells=1)
    obs = np.array([[3.0]])
    post = bayesian_posterior(obs, rm, t_bin=1.0, normalise=False)
    expected = [poisson.pmf(3, 2.0), poisson.pmf(3, 5.0)]
    assert np.allclose(post.ravel(), expected, rtol=1e-6)


def test_bayesian_posterior_normalises_to_one():
    rng = np.random.RandomState(0)
    rm = rng.uniform(0.5, 5.0, size=(5, 5, 3))
    obs = rng.poisson(2, size=(3, 4)).astype(float)
    post = bayesian_posterior(obs, rm, t_bin=0.5, normalise=True)
    assert np.allclose(post.sum(axis=(0, 1)), 1.0)


def test_bayesian_posterior_rejects_cell_mismatch():
    import pytest
    with pytest.raises(ValueError):
        bayesian_posterior(np.ones((2, 5)), np.ones((3, 3, 3)), t_bin=0.5)


def test_spike_counts_per_window_no_overlap():
    spikes = np.array([0, 0, 5, 5, 9, 15])
    counts, win_lim = spike_counts_per_window(spikes, (0, 15), win_width=5)
    assert win_lim.tolist() == [[0, 4], [5, 9], [10, 14]]
    assert counts.tolist() == [2.0, 3.0, 0.0]


def test_spike_counts_per_window_overlap_sums_correctly():
    spikes = np.array([2.0])  # single spike at possamp 2
    counts, win_lim = spike_counts_per_window(spikes, (0, 9), win_width=4, win_overlap=0.5)
    # the spike at possamp 2 should be counted in every window whose [start,end] contains it
    contains_spike = (win_lim[:, 0] <= 2) & (2 <= win_lim[:, 1])
    assert np.array_equal(counts > 0, contains_spike)
    assert counts[contains_spike].tolist() == [1.0] * contains_spike.sum()


def test_spike_counts_per_window_rejects_full_overlap():
    import pytest
    with pytest.raises(ValueError):
        spike_counts_per_window(np.array([0.0]), (0, 10), win_width=4, win_overlap=1.0)


def test_spike_counts_per_window_smoothing_preserves_total_mass():
    spikes = np.array([5.0, 5.0, 5.0])
    counts_raw, _ = spike_counts_per_window(spikes, (0, 19), win_width=20, win_smooth=0.0)
    counts_smooth, _ = spike_counts_per_window(spikes, (0, 19), win_width=20, win_smooth=1.5)
    # a single window spanning the whole range: smoothing shouldn't create or destroy spikes
    assert np.isclose(counts_raw[0], 3.0)
    assert np.isclose(counts_smooth[0], 3.0, atol=0.05)


def test_decode_recovers_position_from_clean_place_fields():
    """The definitive correctness check: clean (non-periodic) place fields with
    known ground truth should decode far better than chance. Real grid cells
    are used elsewhere in this suite (see test_decoding integration below) but
    are a poor standalone correctness check -- their periodicity makes
    absolute-position decoding genuinely ambiguous even with perfect code."""
    rng = np.random.RandomState(0)
    ny, nx, n_cells = 20, 20, 12
    yy, xx = np.mgrid[0:ny, 0:nx]
    centres = rng.uniform(2, ny - 3, size=(n_cells, 2))
    sigma = 2.5
    rate_maps = np.zeros((ny, nx, n_cells))
    for c in range(n_cells):
        cy, cx = centres[c]
        rate_maps[:, :, c] = 15.0 * np.exp(-((yy - cy) ** 2 + (xx - cx) ** 2) / (2 * sigma ** 2)) + 0.2

    t_bin = 0.5
    n_windows = 150
    true_bins = rng.randint(3, ny - 3, size=(n_windows, 2))
    obs_spk = np.zeros((n_cells, n_windows))
    for w in range(n_windows):
        r, c = true_bins[w]
        obs_spk[:, w] = rng.poisson(rate_maps[r, c, :] * t_bin)

    post = bayesian_posterior(obs_spk, rate_maps, t_bin)
    dec_bins = np.zeros((n_windows, 2), dtype=int)
    for w in range(n_windows):
        flat = np.nanargmax(post[:, :, w].ravel(order="F"))
        dec_bins[w] = np.unravel_index(flat, (ny, nx), order="F")

    err_bins = np.sqrt(np.sum((dec_bins - true_bins) ** 2, axis=1))
    chance_level = np.sqrt(ny ** 2 + nx ** 2) / 3  # rough chance-level scale
    assert err_bins.mean() < chance_level / 2


def test_bayesian_posterior_with_prior_normalises_and_is_finite():
    rng = np.random.RandomState(1)
    ny, nx, n_cells, n_time = 10, 10, 3, 25
    rate_maps = rng.uniform(0.5, 5.0, size=(ny, nx, n_cells))
    dwell = rng.uniform(0.5, 2.0, size=(ny, nx))
    obs = rng.poisson(2, size=(n_cells, n_time)).astype(float)
    post = bayesian_posterior_with_prior(obs, rate_maps, dwell, t_bin=0.5)
    assert post.shape == (ny, nx, n_time)
    sums = np.nansum(post, axis=(0, 1))
    assert np.allclose(sums, 1.0, atol=1e-6)


def test_decode_maximum_likelihood_zero_error_when_post_matches_true_bin():
    """If the posterior is a clean delta at the animal's true bin, decode
    error should be ~0."""
    class _FakeRecording:
        pass

    ny, nx, bin_px = 5, 5, 10.0
    n_pos = 200
    rng = np.random.RandomState(2)
    rec = _FakeRecording()
    # place the animal at a fixed known pixel location for every sample
    rec.xy = np.tile([2.5 * bin_px, 3.5 * bin_px], (n_pos, 1))

    win_lim = np.array([[0, 19], [20, 39]])
    post = np.zeros((ny, nx, 2))
    post[3, 2, :] = 1.0  # row=3 (y), col=2 (x) -> matches the fixed true location

    err, dec_xy, true_xy, ml_val = decode_maximum_likelihood(post, rec, bin_px, win_lim)
    assert np.allclose(err, 0.0, atol=1e-9)
    assert np.allclose(ml_val, 1.0)
