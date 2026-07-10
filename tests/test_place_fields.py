"""Tests for place-field stats, spatial coherence, and spatial correlation."""

import numpy as np

from barrylab.analysis._stats import nan_corr
from barrylab.analysis.place_fields import (
    field_stats,
    spatial_coherence,
    spatial_correlation,
)
from barrylab.analysis.spatial import true_ppm_bounding_box


def test_nan_corr_basic_pearson():
    x = np.array([1.0, 2.0, 3.0, 4.0])
    y = np.array([2.0, 4.0, 6.0, 8.0])
    r, p = nan_corr(x, y)
    assert np.isclose(r, 1.0)
    assert p < 0.01


def test_nan_corr_drops_nan_and_inf_pairs():
    x = np.array([1.0, 2.0, np.nan, 4.0, np.inf])
    y = np.array([2.0, 4.0, 6.0, 8.0, 10.0])
    r, p = nan_corr(x, y)
    assert np.isclose(r, 1.0)  # only the 3 finite pairs remain, still perfectly linear


def test_nan_corr_all_invalid_returns_nan():
    x = np.array([np.nan, np.nan])
    y = np.array([1.0, 2.0])
    r, p = nan_corr(x, y)
    assert np.isnan(r) and np.isnan(p)


def test_nan_corr_spearman_and_kendall_run():
    x = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    y = np.array([1.0, 4.0, 9.0, 16.0, 25.0])  # monotonic, not linear
    r_pear, _ = nan_corr(x, y, "pearson")
    r_spear, _ = nan_corr(x, y, "spearman")
    r_kend, _ = nan_corr(x, y, "kendall")
    assert np.isclose(r_spear, 1.0)  # perfectly monotonic
    assert np.isclose(r_kend, 1.0)
    assert r_pear < 1.0  # not perfectly linear


def test_field_stats_silent_cell_is_off():
    rm = np.full((10, 10), np.nan)
    rm[3:6, 3:6] = 0.0
    stats = field_stats(rm, min_field_bins=1, min_rate=0.1)
    assert stats["cell_on"] is False
    assert stats["field_area"] == 0


def test_field_stats_detects_a_clean_field():
    rm = np.full((20, 20), 0.5)  # low background everywhere visited
    rm[8:12, 8:12] = 5.0  # a clear field, mean well above 2x background
    stats = field_stats(rm, min_field_bins=10, min_rate=1.0)
    assert stats["cell_on"] is True
    assert stats["field_area"] == 16
    assert np.isclose(stats["mean_rate"], 5.0)
    assert stats["max_field_area"] == 16


def test_field_stats_rejects_field_below_min_size():
    rm = np.full((20, 20), 0.5)
    rm[8:9, 8:9] = 5.0  # a single-bin "field"
    stats = field_stats(rm, min_field_bins=5, min_rate=1.0)
    assert stats["cell_on"] is False


def test_field_stats_rejects_field_below_min_rate():
    rm = np.full((20, 20), 0.5)
    rm[8:12, 8:12] = 1.05  # above 2x mean but below the min_rate criterion
    stats = field_stats(rm, min_field_bins=10, min_rate=50.0)
    assert stats["cell_on"] is False


def test_spatial_coherence_smooth_field_is_higher_than_noise():
    yy, xx = np.mgrid[0:30, 0:30]
    smooth = np.exp(-((xx - 15) ** 2 + (yy - 15) ** 2) / (2 * 5.0 ** 2))
    rng = np.random.RandomState(0)
    noisy = rng.rand(30, 30)
    assert spatial_coherence(smooth) > spatial_coherence(noisy)


def test_spatial_coherence_handles_unvisited_bins():
    # a non-constant map so the correlation is well-defined (a perfectly
    # uniform map has zero variance -> Pearson r is genuinely undefined/NaN,
    # matching MATLAB's corr() in that degenerate case)
    yy, xx = np.mgrid[0:10, 0:10]
    rm = (xx + 2.0 * yy).astype(float)
    rm[0, 0] = np.nan
    sc = spatial_coherence(rm)
    assert np.isfinite(sc)


def test_spatial_correlation_identical_maps_is_one():
    rng = np.random.RandomState(1)
    rm = rng.rand(15, 15)
    rm[rm < 0.1] = np.nan
    r, n = spatial_correlation(rm, rm.copy())
    assert np.isclose(r, 1.0)
    assert n > 0


def test_spatial_correlation_uncorrelated_maps_near_zero():
    rng = np.random.RandomState(2)
    a = rng.rand(20, 20)
    b = rng.rand(20, 20)
    r, _ = spatial_correlation(a, b)
    assert abs(r) < 0.3


def test_spatial_correlation_excludes_mutual_zeros_and_nans():
    a = np.array([[1.0, 0.0, np.nan], [2.0, 0.0, 3.0]])
    b = np.array([[2.0, 0.0, 5.0], [4.0, 1.0, np.nan]])
    r, n = spatial_correlation(a, b)

    # Independently recompute the expected mask/result (mutual-zero OR either-NaN
    # excluded) to check against, rather than hand-deriving the answer by algebra.
    mask = ((a == 0) & (b == 0)) | np.isnan(a) | np.isnan(b)
    valid = ~mask
    expected_r = np.corrcoef(a[valid], b[valid])[0, 1]
    assert n == int(valid.sum())
    assert np.isclose(r, expected_r)


def test_spatial_correlation_resizes_mismatched_shapes():
    rng = np.random.RandomState(3)
    small = rng.rand(10, 10) + 1.0  # avoid zeros
    # upsample "by hand" as the function would, then correlate against itself
    r, n = spatial_correlation(small, small[:, :8])  # different width
    assert n > 0
    assert -1.0 <= r <= 1.0


def test_spatial_correlation_stack_maps_false_returns_per_cell():
    rng = np.random.RandomState(4)
    a = rng.rand(10, 10, 3) + 1.0
    r, n = spatial_correlation(a, a.copy(), stack_maps=False)
    assert r.shape == (3,)
    assert np.allclose(r, 1.0)


def test_true_ppm_bounding_box_square_region():
    occ = np.zeros((50, 50))
    occ[10:30, 10:30] = 1.0  # a 20x20 occupied square
    ppm, i_range, j_range = true_ppm_bounding_box(occ, bin_size_px=5.0, env_size_m=1.0)
    # filled_area = 400 -> sqrt = 20 -> ppm = 20 * 5 / 1 = 100
    assert np.isclose(ppm, 100.0)
    assert list(i_range) == list(range(10, 29))  # off-by-one preserved (see FINDINGS.md)
    assert list(j_range) == list(range(10, 29))


def test_true_ppm_bounding_box_no_visited_bins():
    occ = np.zeros((10, 10))
    ppm, i_range, j_range = true_ppm_bounding_box(occ, bin_size_px=5.0, env_size_m=1.0)
    assert np.isnan(ppm)
    assert list(i_range) == [] and list(j_range) == []


def test_true_ppm_bounding_box_picks_largest_region():
    occ = np.zeros((30, 30))
    occ[1:3, 1:3] = 1.0     # small region, 4 bins
    occ[10:20, 10:20] = 1.0  # large region, 100 bins
    ppm, i_range, _ = true_ppm_bounding_box(occ, bin_size_px=1.0, env_size_m=1.0)
    assert 10 in i_range  # picked the large region, not the small one
    assert 1 not in i_range
