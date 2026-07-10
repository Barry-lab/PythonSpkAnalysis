"""Tests for the viz/ layer (image_rm.m / tintColorMap.m). Requires the
``plot`` extra (matplotlib); skips cleanly if unavailable.

Per the task brief, this is a plotting function: "runs and produces a
sane-looking artifact" is the whole validation bar, not numerical
ground-truth -- no MATLAB comparison is meaningful for a figure.
"""

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib", reason="matplotlib not installed (barrylab[plot])")
matplotlib.use("Agg")  # headless backend for CI/test environments

from barrylab.viz.colormaps import TINT_CMAP  # noqa: E402
from barrylab.viz.plots import image_rm  # noqa: E402


def test_tint_cmap_has_five_levels_matching_matlab_literal():
    assert TINT_CMAP.N == 5
    expected = [
        (0.00, 0.00, 0.75),
        (0.00, 0.75, 0.90),
        (0.10, 1.00, 0.40),
        (1.00, 1.00, 0.00),
        (0.95, 0.00, 0.00),
    ]
    for i, rgb in enumerate(expected):
        np.testing.assert_allclose(TINT_CMAP(i)[:3], rgb)


def test_image_rm_runs_on_synthetic_ratemap_with_nans():
    rm = np.array([[1.0, 2.0, np.nan], [3.0, np.nan, 5.0], [0.0, 1.0, 2.0]])
    im = image_rm(rm)
    assert im.get_array().shape == rm.shape
    assert im.cmap.N == 5
    alpha = np.asarray(im.get_alpha())
    np.testing.assert_array_equal(alpha, (~np.isnan(rm)).astype(float))
    assert im.axes.get_title() == "5.0Hz"


def test_image_rm_accepts_custom_colormap_and_axes():
    import matplotlib.pyplot as plt

    rm = np.array([[0.0, 1.0], [2.0, 3.0]])
    fig, ax = plt.subplots()
    cmap = plt.get_cmap("viridis", 10)
    im = image_rm(rm, cmap=cmap, ax=ax)
    assert im.axes is ax
    assert im.cmap.N == 10


def test_image_rm_all_unvisited_gives_nan_title():
    rm = np.full((3, 3), np.nan)
    im = image_rm(rm)
    assert im.axes.get_title() == "nanHz"


def test_image_rm_runs_on_real_ratemap(example_trial):
    import barrylab as bl

    rec = bl.load_dacq(example_trial)
    rm = bl.rate_map(rec, 2, 3)
    im = image_rm(rm)
    assert im.get_array().shape == rm.shape
    assert np.isfinite(np.nanmax(rm))
