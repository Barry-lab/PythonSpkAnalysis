"""Tests for head-direction analysis."""

import numpy as np

from barrelab.analysis.direction import (
    circular_mean,
    directional_fft,
    directionality,
    head_direction_map,
    head_direction_stats,
    resultant_vector,
)


def test_circular_mean_simple():
    # mean of 30 and 60 degrees is 45
    assert np.isclose(circular_mean([np.radians(30), np.radians(60)]), np.radians(45))


def test_circular_mean_wraps_across_zero():
    # mean of 350 and 10 degrees is 0 (== 2*pi on the circle), not 180
    m = circular_mean([np.radians(350), np.radians(10)])
    circular_dist_to_zero = min(m, 2 * np.pi - m)
    assert circular_dist_to_zero < 1e-9


def test_circular_mean_weighted():
    # heavy weight on 0 pulls the mean toward 0
    m = circular_mean([0.0, np.pi / 2], weights=[9.0, 1.0])
    assert 0 < m < np.radians(10)


def test_resultant_vector_uniform_is_zero():
    r, _ = resultant_vector(np.ones(60))
    assert r < 1e-9


def test_resultant_vector_delta_is_one():
    tc = np.zeros(60)
    tc[10] = 1.0
    r, preferred = resultant_vector(tc)
    assert np.isclose(r, 1.0)
    assert np.isclose(preferred, (10 + 0.5) * 6.0)  # bin centre, 6 deg bins


def test_resultant_vector_silent_cell():
    r, preferred = resultant_vector(np.zeros(60))
    assert r == 0.0 and np.isnan(preferred)


def test_directionality_uniform_is_zero():
    assert directionality(np.ones(60)) < 1e-6


def test_directionality_peaked_is_high():
    tc = np.zeros(60)
    tc[10] = 1.0
    assert directionality(tc) > 1.0


def test_directional_fft_detects_six_cycle():
    theta = np.arange(360) * (2 * np.pi / 360)
    dir_rate = 5.0 + np.cos(6 * theta)  # a pure 6-cycle (hexadirectional) signal
    power, freqs, rel_freq, rel_power = directional_fft(dir_rate)
    assert len(power) == 181 and len(freqs) == 181
    assert np.isclose(power[0], 0.0, atol=1e-6)   # DC removed
    assert np.argmax(power) == 6                  # 6-cycle component dominates
    assert np.max(rel_power) < 0.01               # 2,3,4,5,7 cycles negligible


def test_head_direction_map_on_example(example_trial):
    import barrelab as bl

    rec = bl.load_dacq(example_trial)
    tc = head_direction_map(rec, 2, 3)
    assert tc.shape == (60,)
    assert np.all(np.isfinite(tc)) and np.all(tc >= 0)

    stats = head_direction_stats(rec, 2, 3)
    assert 0.0 <= stats["resultant_length"] <= 1.0
    assert 0.0 <= stats["preferred_direction"] < 360.0
    assert stats["directionality"] >= 0.0
