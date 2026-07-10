"""Tests for head-direction analysis."""

import numpy as np

from barrylab.analysis.direction import (
    circular_mean,
    directional_fft,
    directionality,
    head_direction_map,
    head_direction_stats,
    led2v1_dir_offset,
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


def test_circular_mean_wrap_convention_same_angle_different_range():
    # a small negative-going mean: [0, 2pi) wraps it near 2pi; (-pi, pi] keeps it negative
    angles = [np.radians(-5), np.radians(5), np.radians(-3)]
    m_pos = circular_mean(angles, wrap="0_2pi")
    m_signed = circular_mean(angles, wrap="-pi_pi")
    assert m_pos > np.radians(300)          # wrapped up near 2*pi
    assert -np.radians(10) < m_signed < 0   # stayed small and negative
    # same physical angle on the circle either way
    assert np.isclose(np.mod(m_pos, 2 * np.pi), np.mod(m_signed, 2 * np.pi))


def test_circular_mean_rejects_bad_wrap():
    import pytest
    with pytest.raises(ValueError):
        circular_mean([0.0, 1.0], wrap="bogus")


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
    import barrylab as bl

    rec = bl.load_dacq(example_trial)
    tc = head_direction_map(rec, 2, 3)
    assert tc.shape == (60,)
    assert np.all(np.isfinite(tc)) and np.all(tc >= 0)

    stats = head_direction_stats(rec, 2, 3)
    assert 0.0 <= stats["resultant_length"] <= 1.0
    assert 0.0 <= stats["preferred_direction"] < 360.0
    assert stats["directionality"] >= 0.0


def test_led2v1_dir_offset_recovers_constant_bias():
    rng = np.random.RandomState(0)
    led1 = rng.uniform(0, 360, 2000)
    true_offset = 12.0
    led2 = np.mod(led1 + true_offset + rng.normal(0, 3.0, 2000), 360)
    mean_offset, peak_offset = led2v1_dir_offset(led2, led1)
    assert abs(mean_offset - true_offset) < 1.0
    assert abs(peak_offset - true_offset) < 2.0  # 1-degree-binned, so a bit coarser


def test_led2v1_dir_offset_recovers_negative_bias_without_wrapping():
    # a negative bias must come back negative (near -8), NOT wrap to ~352 --
    # this is exactly the circular_mean range-convention issue.
    rng = np.random.RandomState(1)
    led1 = rng.uniform(0, 360, 2000)
    led2 = np.mod(led1 - 8.0 + rng.normal(0, 2.0, 2000), 360)
    mean_offset, _ = led2v1_dir_offset(led2, led1)
    assert -10.0 < mean_offset < -6.0


def test_led2v1_dir_offset_zero_bias():
    rng = np.random.RandomState(2)
    led1 = rng.uniform(0, 360, 2000)
    led2 = np.mod(led1 + rng.normal(0, 2.0, 2000), 360)
    mean_offset, peak_offset = led2v1_dir_offset(led2, led1)
    assert abs(mean_offset) < 1.0
    assert abs(peak_offset) < 2.0


def test_led2v1_dir_offset_ignores_nans():
    led1 = np.array([10.0, 20.0, np.nan, 40.0] * 50)
    led2 = np.array([15.0, 25.0, 99.0, 45.0] * 50)  # constant +5 offset where valid
    mean_offset, _ = led2v1_dir_offset(led2, led1)
    assert abs(mean_offset - 5.0) < 0.5
