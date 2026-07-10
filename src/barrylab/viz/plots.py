"""Ratemap plotting. Port of ``image_rm.m`` (+ ``tintColorMap.m``, see
:mod:`barrylab.viz.colormaps`).

Requires the ``plot`` extra (``pip install "barrylab[plot]"`` -- matplotlib).
As with the rest of this ``viz`` package, this module is *not* imported by
the top-level ``barrylab`` package.
"""

import matplotlib.pyplot as plt
import numpy as np

from .colormaps import TINT_CMAP


def image_rm(rm, cmap=None, ax=None):
    """Draw a rate map with the lab's standard colour scheme.

    Port of ``image_rm.m``: draws ``rm`` with ``imagesc``-equivalent scaling
    (matplotlib's ``imshow`` normalises continuous data onto a colormap's
    discrete levels the same way MATLAB's ``imagesc`` does), marks unvisited
    bins (NaN in ``rm``) fully transparent -- reproducing MATLAB's
    ``set(h,'AlphaData',~unvisted)`` -- and titles the axes with the peak
    firing rate to one decimal place (``title([maxRate 'Hz'])``).

    Parameters
    ----------
    rm : ndarray (2d)
        Firing rate per spatial bin; NaN for unvisited bins.
    cmap : matplotlib colormap, optional
        Colour map to draw with. Default: the 5-level tint map
        (:data:`barrylab.viz.colormaps.TINT_CMAP`), matching the MATLAB
        default (``color2use='tintColorMap'`` when not otherwise specified).
        Pass e.g. ``matplotlib.pyplot.get_cmap("viridis", 10)`` for a
        different/finer map (MATLAB's ``image_rm(rm, 'parula(10)')``).
    ax : matplotlib.axes.Axes, optional
        Axes to draw into. A new figure/axes is created if omitted.

    Returns
    -------
    matplotlib.image.AxesImage
    """
    rm = np.asarray(rm, dtype=float)
    if cmap is None:
        cmap = TINT_CMAP
    if ax is None:
        _, ax = plt.subplots()

    unvisited = np.isnan(rm)
    im = ax.imshow(rm, cmap=cmap)
    # Explicit per-pixel alpha (matches MATLAB's AlphaData trick exactly);
    # belt-and-braces alongside matplotlib's own default fully-transparent
    # cmap.set_bad() behaviour for NaN.
    im.set_alpha((~unvisited).astype(float))
    ax.set_aspect("equal")
    ax.set_axis_off()

    peak_rate = np.nan if unvisited.all() else np.nanmax(rm)
    ax.set_title(f"{peak_rate:.1f}Hz")
    return im
