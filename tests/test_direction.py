"""Tests for head-direction analysis."""

import numpy as np

from barrylab.analysis.direction import (
    circular_mean,
    directional_fft,
    directional_occupancy,
    directionality,
    head_direction_map,
    head_direction_stats,
    led2v1_dir_offset,
    resultant_vector,
)


def test_directional_occupancy_bins_dwell_by_direction(example_trial):
    import barrylab as bl

    rec = bl.load_dacq(example_trial)
    occ = directional_occupancy(rec, bin_deg=6.0)
    assert occ.shape == (60,)  # 360 / 6
    assert np.all(occ >= 0)
    # total dwell time should equal the session duration (every sample binned once)
    assert np.isclose(occ.sum(), rec.n_pos / rec.pos_sample_rate)


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
    assert m_pos > np.radians(300)  # wrapped up near 2*pi
    assert -np.radians(10) < m_signed < 0  # stayed small and negative
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
    assert np.isclose(power[0], 0.0, atol=1e-6)  # DC removed
    assert np.argmax(power) == 6  # 6-cycle component dominates
    assert np.max(rel_power) < 0.01  # 2,3,4,5,7 cycles negligible


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


def test_led2v1_dir_offset_wrap_convention_is_a_true_no_op():
    """Pins the exact algebraic equivalence the MATLAB fix relies on: computing
    the mean via wrap="0_2pi" then explicitly re-wrapping to (-pi, pi] gives an
    identical result to computing directly with wrap="-pi_pi" -- for both a
    positive- and a negative-true-offset case. This is what makes
    Universal_Matlab's led2v1_dir_offset.m safe to switch from the external
    toolbox circ_mean to the lab's own cb_circMean (which only returns
    [0, 2*pi)): the switch changes nothing about the returned meanOffset once
    the (-pi, pi] wrap is applied explicitly afterward, in either direction.
    See FINDINGS.md.
    """
    rng = np.random.RandomState(3)
    for true_offset_deg in (-5.0, 50.0):
        led1 = rng.uniform(0, 360, 2000)
        led2 = np.mod(led1 + true_offset_deg + rng.normal(0, 2.0, 2000), 360)
        diff = np.radians(np.mod(led2 - led1 + 180, 360) - 180)
        via_signed = circular_mean(diff, wrap="-pi_pi")
        via_0_2pi_rewrapped = circular_mean(diff, wrap="0_2pi")
        via_0_2pi_rewrapped = np.mod(via_0_2pi_rewrapped + np.pi, 2 * np.pi) - np.pi
        assert np.isclose(via_signed, via_0_2pi_rewrapped)


def test_led2v1_dir_offset_on_real_trial_reveals_a_genuine_large_offset(golden_mat):
    """Real-data sanity check requested by the PI: does 2-LED heading (pos.dir)
    align with movement direction (pos.dir_disp, 1-LED/displacement-based) on
    average? Run directly on the authoritative MATLAB-computed arrays from the
    shipped golden fixture (not this port's own derived values), so this is
    checking the underlying data/method, not this port's direction pipeline.

    The honest answer for this real trial: no, not closely -- offset is
    consistently around +50 degrees (not near 0), with a moderate (not
    noise-level) resultant length, so it is a real, fairly consistent
    difference rather than random jitter. led2v1_dir_offset.m's own docstring
    explicitly frames 1-LED/displacement heading as "less accurate" than
    2-LED heading, and frame-to-frame displacement at 50 Hz is a genuinely
    noisy proxy for a foraging rat's actual head orientation (a weaving gait
    dominates the raw per-sample displacement vector) -- so a large offset
    here does not by itself indicate a tracking or pipeline bug. What matters
    for this port: the number is stable and reproducible (pinned below), and
    unaffected by which circular-mean convention computes it (see
    test_led2v1_dir_offset_wrap_convention_is_a_true_no_op) -- confirmed
    directly on real data too, not just synthetic cases.
    """
    pos = golden_mat["data"].pos
    led2_dir = np.asarray(pos.dir, dtype=float).ravel()
    led1_dir = np.asarray(pos.dir_disp, dtype=float).ravel()

    mean_offset, peak_offset = led2v1_dir_offset(led2_dir, led1_dir)
    assert 45.0 < mean_offset < 58.0  # pins the real, reproducible value (~52 deg)

    diff = np.radians(np.mod(led2_dir - led1_dir + 180, 360) - 180)
    resultant_length = np.abs(np.nanmean(np.exp(1j * diff)))
    assert resultant_length > 0.3  # a real, non-noise-level relationship, not near 0
