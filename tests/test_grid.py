"""Tests for the spatial autocorrelogram and gridness."""

import numpy as np

from barrylab.analysis.grid import (
    gridness,
    grid_properties,
    spatial_autocorrelogram,
)


def _hex_grid_map(n=61, spacing=12.0, sigma=2.5):
    """A synthetic rate map of Gaussian fields on a hexagonal lattice."""
    yy, xx = np.mgrid[0:n, 0:n]
    a1 = np.array([spacing, 0.0])
    a2 = np.array([spacing * 0.5, spacing * np.sqrt(3) / 2])
    centre = np.array([n / 2, n / 2])
    m = np.zeros((n, n))
    for i in range(-5, 6):
        for j in range(-5, 6):
            c = i * a1 + j * a2 + centre
            m += np.exp(-((xx - c[0]) ** 2 + (yy - c[1]) ** 2) / (2 * sigma ** 2))
    return m


def test_sac_central_is_one_and_symmetric():
    rng = np.random.RandomState(0)
    smooth = rng.rand(20, 25)
    sac = spatial_autocorrelogram(smooth)
    assert sac.shape == (39, 49)  # (2n-1)
    cy, cx = sac.shape[0] // 2, sac.shape[1] // 2
    assert np.isclose(sac[cy, cx], 1.0)
    # the central (zero-lag) bin is the global maximum
    assert np.isclose(np.nanmax(sac), sac[cy, cx])
    # an autocorrelogram is point-symmetric about its centre (corners are NaN
    # where the overlap is too small — compare with equal_nan)
    assert np.allclose(sac, sac[::-1, ::-1], atol=1e-9, equal_nan=True)


def test_sac_handles_nan_bins():
    # a structured (non-constant) map with an unvisited bin
    yy, xx = np.mgrid[0:15, 0:15]
    m = (xx + 2 * yy).astype(float)
    m[0, 0] = np.nan
    sac = spatial_autocorrelogram(m)
    cy, cx = sac.shape[0] // 2, sac.shape[1] // 2
    assert np.isclose(sac[cy, cx], 1.0)  # NaN ignored, central still defined


def test_gridness_of_synthetic_hex_grid_is_high():
    sac = spatial_autocorrelogram(_hex_grid_map())
    standard, expanding = gridness(sac)
    assert standard > 0.5, "a clean hexagonal grid should score high gridness"
    assert expanding > 0.5


def test_grid_properties_of_hex_grid():
    sac = spatial_autocorrelogram(_hex_grid_map(spacing=12.0))
    props = grid_properties(sac)
    # six peaks at ~one lattice spacing from the centre
    assert 9 < props["scale"] < 15
    assert props["gridness"] > 0.5
    # three axes, ~60 deg apart
    assert len(props["orientation"]) == 3


def test_gridness_of_blob_is_not_high():
    # a single field has no hexagonal structure
    yy, xx = np.mgrid[0:41, 0:41]
    blob = np.exp(-((xx - 20) ** 2 + (yy - 20) ** 2) / (2 * 5.0 ** 2))
    standard, _ = gridness(spatial_autocorrelogram(blob))
    assert not (standard > 0.5)  # nan or low, but never grid-like


def test_gridness_discriminates_example_grid_cells(example_trial):
    """The example README lists t2c1 and t2c3 as grid cells, t2c4/t2c5 as not."""
    import barrylab as bl

    rec = bl.load_dacq(example_trial)
    grid_cells = {cl: bl.grid_score(rec, 2, cl)[0] for cl in (1, 3)}
    non_grid = {cl: bl.grid_score(rec, 2, cl)[0] for cl in (4, 5)}
    assert min(grid_cells.values()) > 0.4, grid_cells
    assert max(non_grid.values()) < min(grid_cells.values()), (grid_cells, non_grid)
