"""Tests for elliptical grid-distortion fitting and correction."""

import numpy as np

from barrylab.analysis.grid import fast_gridness, spatial_autocorrelogram
from barrylab.analysis.grid_ellipse import (
    _fit_ellipse,
    brandon_regularise,
    fit_grid_ellipse,
    regularise_elliptic_sac,
)


def _hex_grid_map(n=81, spacing=10.0, sigma=2.0, stretch_x=1.0):
    yy, xx = np.mgrid[0:n, 0:n]
    a1 = np.array([spacing * stretch_x, 0.0])
    a2 = np.array([spacing * 0.5 * stretch_x, spacing * np.sqrt(3) / 2])
    centre = np.array([n / 2.0, n / 2.0])
    m = np.zeros((n, n))
    for i in range(-5, 6):
        for j in range(-5, 6):
            c = i * a1 + j * a2 + centre
            m += np.exp(-((xx - c[0]) ** 2 + (yy - c[1]) ** 2) / (2 * sigma ** 2))
    return m


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


# --------------------------------------------------------------------------
# brandon_regularise (brandonRegEcc.m) -- the MATLAB source's own docstring
# calls this method "clunky" and admits "at least one error in the
# equations"; validated structurally rather than against a hand-derived
# analytical ground truth, matching how the (better) least-squares ellipse
# correction was itself validated (FINDINGS.md).
# --------------------------------------------------------------------------
def test_brandon_regularise_too_few_peaks_returns_input_unchanged():
    # A strictly monotonic radial ramp gives exactly one detected peak (see
    # the identical construction/rationale in test_grid.py and
    # test_grid_phase.py -- a flat/tied background trivially satisfies
    # equality-based local-max detection almost everywhere, so it does NOT
    # reliably produce "too few peaks").
    n = 21
    yy, xx = np.mgrid[0:n, 0:n]
    c = n // 2
    ramp = -np.hypot(xx - c, yy - c) + 0.001 * (xx + 1.3 * yy)
    r1, r2 = brandon_regularise(ramp)
    assert np.array_equal(r1, ramp)
    assert np.array_equal(r2, ramp)


def test_brandon_regularise_shape_and_no_crash_on_real_grid_cell(example_trial):
    import barrylab as bl

    rec = bl.load_dacq(example_trial)
    rm = bl.rate_map(rec, 2, 3)
    sac = bl.spatial_autocorrelogram(rm)
    r1, r2 = brandon_regularise(sac)
    assert r1.shape == sac.shape
    assert r2.shape == sac.shape
    assert np.mean(np.isnan(r1)) < 0.2
    assert np.mean(np.isnan(r2)) < 0.2


def test_brandon_regularise_regular_hex_grid_gives_near_identity_method2():
    # For a perfectly regular (non-elliptical) hexagonal lattice, the six
    # central peaks are all equidistant from the centre, so closeR == distR
    # (== R, say). Method 2's minor-axis formula then reduces exactly:
    # sin^2(theta) / (1/R^2 - cos^2(theta)/R^2) = sin^2(theta) / (sin^2(theta)/R^2) = R^2,
    # i.e. m2MinR == m2MajR == R algebraically, regardless of theta -- so
    # regularise_elliptic_sac's own a==b shortcut should fire and return the
    # SAC completely unchanged. This is a real, checkable mathematical
    # property of the method (not just "looks about right").
    sac = spatial_autocorrelogram(_hex_grid_map(stretch_x=1.0))
    _, r2 = brandon_regularise(sac)
    both_nan = np.isnan(r2) & np.isnan(sac)
    matches = np.isclose(r2, sac, atol=1e-9) | both_nan
    assert matches.all()


def test_brandon_regularise_elongated_hex_grid_improves_gridness_for_a_method():
    # An elongated (elliptical) grid should have poor gridness before
    # correction; at least one of the two Brandon methods should improve it
    # -- mirroring how the least-squares ellipse correction was validated
    # (it didn't reliably improve gridness on both real test cells either;
    # this method is explicitly flagged by its own MATLAB docstring as less
    # reliable still).
    sac = spatial_autocorrelogram(_hex_grid_map(stretch_x=1.6))
    g_before, _, _ = fast_gridness(sac)
    r1, r2 = brandon_regularise(sac)
    g1, _, _ = fast_gridness(r1)
    g2, _, _ = fast_gridness(r2)
    assert max(g1, g2) > g_before
