"""Tests for the general N-D histnd binning and the place x direction (pxd) model."""

import numpy as np

from barrylab.analysis._binning import histnd
from barrylab.analysis.pxd import pxd


def test_histnd_2d_matches_hand_computation():
    data = np.array([[0.5, 0.5], [0.5, 1.5], [2.5, 2.5]])
    counts, bins = histnd(data, np.array([3.0, 3.0]), 1.0)
    assert counts.shape == (3, 3)
    assert counts.sum() == 3
    assert counts[0, 0] == 1
    assert counts[0, 1] == 1
    assert counts[2, 2] == 1
    assert bins.shape == (3, 2)


def test_histnd_3d_basic():
    # 2 direction bins x 2x2 position bins; one point per bin
    data = np.array([
        [0.5, 0.5, 0.5],
        [0.5, 0.5, 1.5],
        [1.5, 1.5, 0.5],
        [1.5, 1.5, 1.5],
    ])
    counts, _ = histnd(data, np.array([2.0, 2.0, 2.0]), 1.0)
    assert counts.shape == (2, 2, 2)
    assert counts.sum() == 4
    assert counts[0, 0, 0] == 1 and counts[1, 1, 1] == 1


def test_histnd_drops_nan_rows():
    data = np.array([[1.0, 1.0, 1.0], [np.nan, 1.0, 1.0]])
    counts, _ = histnd(data, np.array([3.0, 3.0, 3.0]), 1.0)
    assert counts.sum() == 1


def test_histnd_clamps_out_of_range():
    data = np.array([[99.0, 99.0], [-5.0, -5.0]])
    counts, _ = histnd(data, np.array([3.0, 3.0]), 1.0)
    assert counts[-1, -1] == 1
    assert counts[0, 0] == 1


def test_pxd_recovers_separable_place_and_direction_tuning():
    # Build a synthetic cell whose true rate factorises as p(pos) * d(dir):
    # a single spatial hotspot combined with a directional preference.
    rng = np.random.RandomState(0)
    n_dir, n_y, n_x = 12, 8, 8

    true_p = np.ones((n_y, n_x)) * 0.5
    true_p[3:5, 3:5] = 4.0  # a hotspot
    true_d = np.ones(n_dir)
    true_d[0] = 3.0  # preferred direction

    times = rng.uniform(0.5, 1.5, size=(n_dir, n_y, n_x))  # uneven but full coverage
    expected = true_p[None, :, :] * true_d[:, None, None] * times
    spikes = rng.poisson(expected)

    p, d = pxd(spikes, times)
    assert p.shape == (n_y, n_x)
    assert d.shape == (n_dir,)
    assert np.all(np.isfinite(p)) and np.all(np.isfinite(d))

    # recovered maps should correlate strongly with the ground truth
    assert np.corrcoef(p.ravel(), true_p.ravel())[0, 1] > 0.8
    assert np.corrcoef(d, true_d)[0, 1] > 0.7

    # normalisation: pxd.m rescales p and d SEPARATELY so each marginal's
    # projected total matches total spikes (not the full 3-D joint product,
    # which the algorithm doesn't constrain once p and d are independently
    # rescaled at the end).
    total_spikes = spikes.sum()
    assert np.isclose(np.sum(p * times.sum(axis=0)), total_spikes)
    assert np.isclose(np.sum(d * times.sum(axis=(1, 2))), total_spikes)


def test_pxd_rejects_mismatched_shapes():
    import pytest
    with pytest.raises(ValueError):
        pxd(np.ones((2, 3, 3)), np.ones((2, 4, 4)))


def test_place_by_direction_on_example(example_trial):
    import barrylab as bl

    rec = bl.load_dacq(example_trial)
    p, d = bl.place_by_direction(rec, 2, 3)
    assert p.ndim == 2 and d.ndim == 1
    assert np.any(np.isfinite(p)) and np.all(np.isfinite(d))

    # sanity: the pxd place map should correlate with the plain rate map at a
    # matching bin size, and the direction map with the plain tuning curve --
    # both already independently validated.
    rm = bl.rate_map(rec, 2, 3, bin_cm=4.0)
    if rm.shape == p.shape:
        both = ~np.isnan(p) & ~np.isnan(rm)
        assert np.corrcoef(p[both], rm[both])[0, 1] > 0.7

    hd = bl.head_direction_map(rec, 2, 3, bin_deg=6.0, smooth_size=3)
    both_d = np.isfinite(d) & np.isfinite(hd)
    assert np.corrcoef(d[both_d], hd[both_d])[0, 1] > 0.6


def test_place_by_direction_uses_correct_x_extent_on_non_square_window():
    # wrapper_pxd.m used to have a copy-paste bug where the spatial bin range
    # used max(y) for BOTH the y and x extents -- invisible on the shipped
    # example trial (whose tracked window happens to be near-square) so this
    # deliberately non-square synthetic case is needed to actually exercise
    # the fix. X spans 0-300, Y spans 0-100 -- a bug reusing max(y)=100 for
    # the x extent would badly truncate the x-axis bin count/coverage.
    import barrylab as bl
    from barrylab.core.recording import Electrode, Recording

    rng = np.random.RandomState(0)
    n_pos = 4000
    x = rng.uniform(0, 300, n_pos)
    y = rng.uniform(0, 100, n_pos)
    xy = np.column_stack([x, y])
    direction = rng.uniform(0, 360, n_pos)

    spk_idx = rng.choice(n_pos, size=300, replace=False)
    elec = Electrode(identifier=1, spike_times=spk_idx / 50.0, pos_index=spk_idx,
                     clusters=np.ones(300, dtype=int))
    rec = Recording(name="synthetic", xy=xy, direction=direction,
                    pos_sample_rate=50.0, pixels_per_metre=100.0,
                    window=(0.0, 300.0, 0.0, 100.0), electrodes={1: elec})

    p, d = bl.place_by_direction(rec, 1, 1, bin_cm=10.0)
    ny, nx = p.shape
    # x extent (300 px = 3 m = 30cm bins of 10cm) should give ~3x the number
    # of x bins as the y extent (100 px = 1 m) -- NOT equal, which is what
    # the max(y)-reused-for-x bug would have produced.
    assert nx > ny * 2
