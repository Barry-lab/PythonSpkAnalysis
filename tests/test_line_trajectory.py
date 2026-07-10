"""Tests for best-fit-line replay trajectory scoring (``lineTraj_*.m``)."""

import numpy as np

from barrylab.analysis.line_trajectory import best_line, line_residual, reget_line_band


def _extract_pmat(event):
    """Mirrors the extraction recipe in ``pMat_data.mat READ ME.txt``."""
    nan_ind = np.flatnonzero(~np.isnan(event[0, :]))
    with np.errstate(invalid="ignore"):
        zero_ind = np.flatnonzero(event.sum(axis=0) > 0)
    bins2use = np.intersect1d(nan_ind, zero_ind)
    if bins2use.size == 0:
        return None
    pmat = event[:, bins2use.min() : bins2use.max() + 1].copy()
    pmat[np.isnan(pmat)] = 0.0
    return pmat


def test_best_line_roundtrips_via_reget_and_residual(pmat_data):
    """This is exactly what ``lineTraj_testScript.m`` checks: rebuilding the
    best line's band from its own reported gradient/intercept and rescoring
    it against the pMat must reproduce the original residual. No MATLAB
    runtime is available on this machine to generate independent reference
    numbers from this fixture, so this end-to-end self-consistency check
    (which a bug in any of the three ported functions would very likely
    break) is the strongest available oracle -- see FINDINGS.md.
    """
    n_events = pmat_data.shape[2]
    n_checked = 0
    for i in range(n_events):
        pmat = _extract_pmat(pmat_data[:, :, i])
        if pmat is None:
            continue
        best_res, best_spd, best_y_interc, _ = best_line(pmat, 0.005, 2, 14)
        if np.isnan(best_res):
            continue
        n_checked += 1
        ln_bnd = reget_line_band(pmat.shape, best_spd[0], best_y_interc, 14)
        res_test = line_residual(pmat, ln_bnd)
        assert abs(best_res - res_test) < 1e-9, "event %d: %r vs %r" % (i, best_res, res_test)

    assert n_checked > 50  # sanity: most events should yield a fit


def test_best_line_residual_in_unit_range(pmat_data):
    for i in range(20):
        pmat = _extract_pmat(pmat_data[:, :, i])
        if pmat is None:
            continue
        best_res, _, _, _ = best_line(pmat, 0.005, 2, 14)
        if np.isnan(best_res):
            continue
        assert 0.0 <= best_res <= 1.0


def test_best_line_zero_gradient_line():
    # A pMat with all its mass on a single flat row should be found as
    # ~gradient-0, high-residual.
    pmat = np.zeros((50, 30))
    pmat[25, :] = 1.0
    best_res, best_spd, best_y_interc, best_band = best_line(pmat, 0.005, 2, 14)
    assert best_res > 0.9
    assert best_band[25, :].sum() > 0


def test_line_residual_full_coverage_gives_one():
    pmat = np.random.RandomState(0).rand(10, 10)
    band = np.ones_like(pmat)
    assert np.isclose(line_residual(pmat, band), 1.0)


def test_line_residual_no_coverage_gives_zero():
    pmat = np.random.RandomState(0).rand(10, 10)
    band = np.zeros_like(pmat)
    assert line_residual(pmat, band) == 0.0


def test_reget_line_band_is_binary_mask():
    band = reget_line_band((40, 20), spd=1.5, y_interc=10, y_range=14)
    assert band.shape == (40, 20)
    assert set(np.unique(band)).issubset({0.0, 1.0})


def test_best_line_too_thin_to_fit_any_line_returns_nan():
    # A single-row pMat: every tested (nonzero) gradient needs several rows
    # of vertical extent to place its per-column points, so no line can
    # collect enough of its points on the pMat's one real row to clear the
    # minimum-length cutoff. best_res ends up exactly 0 (not via 0/0), which
    # correctly triggers the explicit "no fit" -> nan substitution.
    pmat = np.random.RandomState(2).rand(1, 50)
    best_res, best_spd, best_y_interc, best_band = best_line(pmat, 0.005, 2, 14)
    assert np.isnan(best_res)
    assert np.all(np.isnan(best_spd))
    assert np.isnan(best_y_interc)


def test_best_line_all_zero_pmat_leaves_spd_uninitialised_nan_quirk():
    # Documents a genuine MATLAB quirk (lineTraj_decode.m), faithfully
    # replicated rather than fixed -- see FINDINGS.md. The code's own
    # comment claims "if there is no data... bestRes will be 0", and an
    # explicit `if bstRes == 0` branch exists to null out bstSpd/bstYInterc
    # for that case. But for a literally all-zero pMat, bstRes is computed
    # as 0/sum(pMat)=0/0, which is NaN under IEEE 754 -- not 0 -- so the
    # safety branch's `== 0` check silently fails to trigger. bstRes still
    # ends up NaN (a caller checking isnan(bstRes) is unaffected), but
    # bstSpd/bstYInterc are left as whatever the (meaningless, tied) first
    # candidate gradient happened to be, rather than also being nan.
    pmat = np.zeros((50, 30))
    with np.errstate(invalid="ignore"):
        best_res, best_spd, best_y_interc, best_band = best_line(pmat, 0.005, 2, 14)
    assert np.isnan(best_res)
    assert not np.any(np.isnan(best_spd))
    assert not np.isnan(best_y_interc)
