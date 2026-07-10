"""Tests for circular-linear regression (theta phase precession)."""

import numpy as np

from barrylab.analysis.phase_precession import _kempter_corr, circular_linear_regression

# Fisher (1993) dataset B.20 (p.252), example 6.9 (p.160) -- embedded in
# circ_regress_cb.m's own docstring as a worked example.
_FISHER_X = np.array(
    [
        107,
        46,
        33,
        67,
        122,
        69,
        43,
        30,
        12,
        25,
        37,
        69,
        5,
        83,
        68,
        38,
        21,
        1,
        71,
        60,
        71,
        71,
        57,
        53,
        38,
        70,
        7,
        48,
        7,
        21,
        27,
    ],
    dtype=float,
)
_FISHER_T_DEG = np.array(
    [
        67,
        66,
        74,
        61,
        58,
        60,
        100,
        89,
        171,
        166,
        98,
        60,
        197,
        98,
        86,
        123,
        165,
        133,
        101,
        105,
        71,
        84,
        75,
        98,
        83,
        71,
        74,
        91,
        38,
        200,
        56,
    ],
    dtype=float,
)
_FISHER_T = _FISHER_T_DEG * np.pi / 180.0


def test_kempter_corr_matches_fisher_book_value():
    """The book's R=-0.316 comes from correlating 2*atan(x) directly against
    t (Fisher's own circular transform of the linear variable) -- NOT from
    circ_regress_cb's own optimally-fitted-slope transform, which is a
    different quantity for this (non-phase-precession) dataset. This test
    validates the underlying circular-circular correlation formula/circular
    mean machinery directly against the textbook oracle."""
    theta = 2 * np.arctan(_FISHER_X)
    r = _kempter_corr(theta, _FISHER_T)
    assert np.isclose(r, -0.316, atol=1e-3)


def test_circular_linear_regression_recovers_known_precession():
    rng = np.random.RandomState(0)
    n = 200
    x = rng.uniform(0, 1, n)
    true_slope = -2 * np.pi  # one full cycle of precession across the field
    t = np.mod(true_slope * x + 0.5 + rng.normal(0, 0.15, n), 2 * np.pi)

    result = circular_linear_regression(x, t, k=200, pp_dir=0, random_state=1)
    assert np.isclose(result["slope_opt"], true_slope, rtol=0.05)
    assert result["corr"] < -0.9  # strong negative circular-circular correlation
    assert result["psim"] < 0.01  # highly significant


def test_circular_linear_regression_no_precession_is_less_significant_than_real():
    # A single random null draw can legitimately land p<0.05 ~5% of the time
    # by chance -- that's what a calibrated permutation test means, so don't
    # assert an absolute threshold on one draw. Instead check the comparative,
    # reliable invariant: genuine precession is far more significant than none.
    rng = np.random.RandomState(2)
    n = 100
    x = rng.uniform(0, 1, n)
    t_random = rng.uniform(0, 2 * np.pi, n)  # phase independent of position
    t_precessing = np.mod(-4.0 * x + rng.normal(0, 0.1, n), 2 * np.pi)

    result_random = circular_linear_regression(x, t_random, k=300, pp_dir=0, random_state=3)
    result_real = circular_linear_regression(x, t_precessing, k=300, pp_dir=0, random_state=3)

    assert 0.0 <= result_random["psim"] <= 1.0
    assert result_real["psim"] < result_random["psim"]
    assert result_real["psim"] < 0.01


def test_circular_linear_regression_too_few_points_returns_nan():
    result = circular_linear_regression(np.array([1.0, 2.0]), np.array([0.1, 0.2]), k=10)
    assert all(np.isnan(v) for v in result.values())


def test_circular_linear_regression_drops_nan_pairs():
    rng = np.random.RandomState(4)
    n = 100
    x = rng.uniform(0, 1, n)
    t = np.mod(-3.0 * x + rng.normal(0, 0.1, n), 2 * np.pi)
    x_with_nan = x.copy()
    x_with_nan[:5] = np.nan
    result_clean = circular_linear_regression(x[5:], t[5:], k=100, random_state=5)
    result_nan = circular_linear_regression(x_with_nan, t, k=100, random_state=5)
    assert np.isclose(result_clean["slope_opt"], result_nan["slope_opt"])


def test_circular_linear_regression_rejects_out_of_range_phase():
    import pytest

    with pytest.raises(ValueError):
        circular_linear_regression(np.array([1.0, 2.0, 3.0]), np.array([1.0, 2.0, 10.0]), k=10)


def test_circular_linear_regression_pp_dir_constrains_sign():
    rng = np.random.RandomState(6)
    n = 100
    x = rng.uniform(0, 1, n)
    t = rng.uniform(0, 2 * np.pi, n)
    result_pos = circular_linear_regression(x, t, k=50, pp_dir=1, random_state=7)
    result_neg = circular_linear_regression(x, t, k=50, pp_dir=-1, random_state=7)
    assert result_pos["slope_opt"] >= 0
    assert result_neg["slope_opt"] <= 0
