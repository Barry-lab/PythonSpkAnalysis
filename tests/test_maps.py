"""Synthetic (data-free) invariant tests for binning and smoothing."""

import numpy as np
import pytest

from barrylab.analysis.maps import _histnd_2d, smooth_ratemap


def test_histnd_basic_counts():
    # points at (y, x); bin size 1 -> bin index == ceil(value)
    data = np.array([[0.5, 0.5], [0.5, 1.5], [2.5, 2.5]])
    out = _histnd_2d(data, [3, 3], 1.0)
    assert out.shape == (3, 3)
    assert out.sum() == 3
    assert out[0, 0] == 1  # (0.5,0.5) -> bin (1,1)
    assert out[0, 1] == 1  # (0.5,1.5) -> bin (1,2)
    assert out[2, 2] == 1  # (2.5,2.5) -> bin (3,3)


def test_histnd_clamps_out_of_range():
    data = np.array([[99.0, 99.0], [-5.0, -5.0]])
    out = _histnd_2d(data, [3, 3], 1.0)
    assert out[-1, -1] == 1  # high point clamps to last bin
    assert out[0, 0] == 1  # low point clamps to first bin


def test_histnd_drops_nan_rows():
    data = np.array([[1.0, 1.0], [np.nan, 1.0]])
    out = _histnd_2d(data, [3, 3], 1.0)
    assert out.sum() == 1


def test_smooth_ratemap_uniform_gives_constant_rate():
    # every bin visited for 0.02 s with exactly 1 spike -> 50 Hz everywhere
    dwell = np.full((5, 5), 0.02)
    spikes = np.ones((5, 5))
    rate, _, _ = smooth_ratemap(dwell, spikes, kern_size=3)
    assert np.allclose(rate, 50.0)


def test_smooth_ratemap_marks_unvisited_nan():
    dwell = np.ones((4, 4))
    dwell[0, 0] = 0  # unvisited
    spikes = np.ones((4, 4))
    rate, _, _ = smooth_ratemap(dwell, spikes, kern_size=3)
    assert np.isnan(rate[0, 0])
    assert np.isfinite(rate[2, 2])


def test_smooth_ratemap_gaussian_uniform_gives_constant_rate():
    # edge normalisation makes a uniform field exact for the gaussian too
    dwell = np.full((6, 6), 0.02)
    spikes = np.ones((6, 6))
    rate, _, _ = smooth_ratemap(dwell, spikes, kern_size=2, kern_type="gaussian")
    assert np.allclose(rate, 50.0)


def test_smooth_ratemap_gaussian_spreads_a_spike():
    dwell = np.full((15, 15), 0.02)
    spikes = np.zeros((15, 15))
    spikes[7, 7] = 10  # a single hot bin
    tight, _, _ = smooth_ratemap(dwell, spikes, kern_size=1, kern_type="gaussian")
    wide, _, _ = smooth_ratemap(dwell, spikes, kern_size=3, kern_type="gaussian")
    assert np.all(np.isfinite(tight)) and np.all(np.isfinite(wide))
    # the peak stays over the hot bin, and more smoothing lowers/spreads it
    assert np.unravel_index(np.nanargmax(tight), tight.shape) == (7, 7)
    assert wide[7, 7] < tight[7, 7]


def test_smooth_ratemap_rejects_bad_kernel():
    with pytest.raises(ValueError):
        smooth_ratemap(np.ones((3, 3)), np.ones((3, 3)), kern_type="triangle")
