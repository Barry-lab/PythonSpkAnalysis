"""Tests for elliptical grid-distortion fitting and correction."""

import numpy as np

from barrylab.analysis.grid_ellipse import _fit_ellipse, fit_grid_ellipse, regularise_elliptic_sac


def _matlab_rotation(phi):
    """The rotation-matrix convention gridEllipse_fit.m's own subfunction
    uses internally (R = [cos sin; -sin cos], the transpose of a standard
    CCW rotation matrix) -- needed to construct a synthetic ellipse whose
    recovered `phi` can be checked bit-exact against a known value."""
    return np.array([[np.cos(phi), np.sin(phi)], [-np.sin(phi), np.cos(phi)]])


def test_fit_ellipse_recovers_known_ellipse_exactly():
    true_a, true_b, true_phi = 10.0, 6.0, np.radians(30.0)
    theta = np.linspace(0, 2 * np.pi, 30, endpoint=False)
    pts = _matlab_rotation(true_phi) @ np.vstack([true_a * np.cos(theta), true_b * np.sin(theta)])
    x, y = pts[0] + 50, pts[1] + 50

    fit = _fit_ellipse(x, y)
    assert np.isclose(fit["a"], true_a, atol=1e-6)
    assert np.isclose(fit["b"], true_b, atol=1e-6)
    assert np.isclose(fit["phi"], true_phi, atol=1e-6)


def test_fit_ellipse_recovers_a_circle():
    theta = np.linspace(0, 2 * np.pi, 20, endpoint=False)
    x, y = 5 * np.cos(theta), 5 * np.sin(theta)
    fit = _fit_ellipse(x, y)
    assert np.isclose(fit["a"], fit["b"], atol=1e-6)


def test_fit_grid_ellipse_on_example_grid_cell(example_trial):
    import barrylab as bl

    rec = bl.load_dacq(example_trial)
    for cluster in (1, 3):
        rm = bl.rate_map(rec, 2, cluster)
        sac = bl.spatial_autocorrelogram(rm)
        result = fit_grid_ellipse(sac)
        assert 0.0 <= result["eccentricity"] < 1.0
        assert result["ab_scale"][0] >= result["ab_scale"][1]  # a (major) >= b (minor)
        assert np.all(np.isfinite(result["xy_scale"]))


def test_fit_grid_ellipse_nan_on_blob_without_six_peaks():
    yy, xx = np.mgrid[0:41, 0:41]
    blob = np.exp(-((xx - 20) ** 2 + (yy - 20) ** 2) / (2 * 5.0 ** 2))
    result = fit_grid_ellipse(blob)
    assert np.isnan(result["eccentricity"])


def test_regularise_elliptic_sac_circular_input_is_unchanged():
    rng = np.random.RandomState(0)
    sac = rng.rand(21, 21)
    out = regularise_elliptic_sac(sac, ab_scale=(5.0, 5.0), orient_deg=17.0)
    assert np.array_equal(out, sac)


def test_regularise_elliptic_sac_runs_on_example_grid_cell(example_trial):
    import barrylab as bl

    rec = bl.load_dacq(example_trial)
    rm = bl.rate_map(rec, 2, 3)
    sac = bl.spatial_autocorrelogram(rm)
    ellipse = bl.fit_grid_ellipse(sac)
    corrected = bl.regularise_elliptic_sac(sac, ellipse["ab_scale"], np.degrees(ellipse["orientation"]))

    assert corrected.shape == sac.shape
    assert np.mean(np.isnan(corrected)) < 0.2  # mostly real values, not collapsed to NaN
    # gridness on the corrected map should still be a sane, finite number
    g, _ = bl.gridness(corrected)
    assert np.isfinite(g)
