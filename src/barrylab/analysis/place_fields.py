"""Place-field and spatial-stability metrics.

Ports the remaining ``Universal_Matlab/Analysis_Spatial_Field`` functions:

- :func:`field_stats` — basic place-field metrics: field size, in-field rate,
  and whether the cell meets the lab's criteria for an "active" place cell
  (``pcFieldStats.m``).
- :func:`spatial_coherence` — correlation of each bin's rate with its
  neighbours' mean rate, a smoothness/signal-quality measure (``spatialCoherance.m``).
- :func:`spatial_correlation` — Pearson correlation between two rate maps
  (optionally resized to match), the standard remapping/stability measure
  (``spatialCorrelation.m``).
"""

import numpy as np
from scipy.ndimage import convolve, label, map_coordinates


def field_stats(rate_map_smoothed, min_field_bins, min_rate):
    """Basic place-field metrics from a smoothed rate map. Port of ``pcFieldStats.m``.

    A field is a contiguous (4-connected) region above twice the mean rate,
    covering at least ``min_field_bins`` bins.

    Parameters
    ----------
    rate_map_smoothed : ndarray
        Smoothed rate map; unvisited bins should be NaN.
    min_field_bins : int
        Minimum field size, in bins, to count as a valid field.
    min_rate : float
        Minimum mean in-field rate for the cell to be considered active.

    Returns
    -------
    dict with keys:
        ``cell_on`` (bool) — meets both the field-size and rate criteria,
        ``field_area`` (int) — combined size of valid fields, in bins,
        ``mean_rate`` (float) — mean rate within valid fields,
        ``max_field_area`` (int) — size of the largest valid field.
    """
    rm = np.asarray(rate_map_smoothed, dtype=float)
    off = {"cell_on": False, "field_area": 0, "mean_rate": 0.0, "max_field_area": 0}

    if np.nansum(rm) == 0:
        return off

    field_thresh = 2.0 * np.nanmean(rm)
    with np.errstate(invalid="ignore"):  # NaN > thresh is correctly False
        above = rm > field_thresh
    labelled, n = label(above)  # scipy default = 4-connectivity in 2D
    if n == 0:
        return off

    areas = np.bincount(labelled.ravel())[1:]  # areas[i-1] = size of label i
    valid = np.where(areas >= min_field_bins)[0] + 1  # 1-based label ids
    if valid.size == 0:
        return off

    in_field = np.isin(labelled, valid)
    mean_rate = float(np.mean(rm[in_field]))
    if mean_rate <= min_rate:
        return off

    return {
        "cell_on": True,
        "field_area": int(areas[valid - 1].sum()),
        "mean_rate": mean_rate,
        "max_field_area": int(areas[valid - 1].max()),
    }


def spatial_coherence(rate_map_unsmoothed):
    """Correlation of each bin's rate with the mean of its 8 neighbours.

    Port of ``spatialCoherance.m`` — a measure of spatial smoothness/signal
    quality (Muller & Kubie 1989). Should be run on an **unsmoothed** rate map.
    Unvisited bins (NaN) are excluded from the correlation but do contribute to
    (dilute) their neighbours' means where they are themselves a neighbour of a
    visited bin.

    Returns
    -------
    float
        The Pearson coefficient (not Fisher z-transformed — matches Tint).
        ``nan`` if fewer than 2 valid bins remain.
    """
    rm = np.asarray(rate_map_unsmoothed, dtype=float)
    visited = ~np.isnan(rm)
    rm_filled = np.where(visited, rm, 0.0)

    neighbour_kernel = np.ones((3, 3))
    neighbour_kernel[1, 1] = 0.0
    neighbour_sum = convolve(rm_filled, neighbour_kernel, mode="constant", cval=0.0)
    n_visited_neigh = convolve(visited.astype(float), neighbour_kernel, mode="constant", cval=0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        mean_neigh = neighbour_sum / n_visited_neigh

    valid = visited & np.isfinite(mean_neigh)
    if valid.sum() < 2:
        return float("nan")
    return float(np.corrcoef(rm[valid], mean_neigh[valid])[0, 1])


def _resize_bilinear(arr, out_shape):
    """Resize a 2-D array by bilinear interpolation over its pixel-index grid.

    Equivalent to MATLAB ``interp2(arr, linspace(1,w,ow), linspace(1,h,oh)',
    '*linear')``: NaNs propagate into any interpolated value that touches them
    (matches the source comment "when interpolating nans will spread - this is
    fine"). A no-op if ``out_shape`` already matches ``arr.shape`` (the common,
    exactly-validated case).
    """
    if arr.shape == tuple(out_shape):
        return arr.astype(float, copy=True)
    src_h, src_w = arr.shape
    out_h, out_w = out_shape
    yi = np.linspace(0, src_h - 1, out_h)
    xi = np.linspace(0, src_w - 1, out_w)
    yy, xx = np.meshgrid(yi, xi, indexing="ij")
    return map_coordinates(arr, [yy, xx], order=1, mode="nearest")


def spatial_correlation(rate_map1, rate_map2, stack_maps=True):
    """Pearson correlation between two rate maps. Port of ``spatialCorrelation.m``.

    Maps that differ in shape are bilinearly resized to their largest common
    shape first (a no-op when they already match). Bins that are mutually zero
    or where either map is unvisited (NaN) are excluded — standard practice for
    place-cell remapping/stability analysis.

    Parameters
    ----------
    rate_map1, rate_map2 : ndarray
        2-D (single map) or 3-D (stacked, last axis = cell) rate maps; the
        stacked axis must match in length between the two.
    stack_maps : bool
        If ``True`` (default), all cells are pooled into one correlation. If
        ``False``, each cell is correlated separately.

    Returns
    -------
    (spatial_corr, n_valid_bins)
        ``spatial_corr`` is a scalar if ``stack_maps`` else an array (one value
        per cell, ``nan`` where fewer than 2 valid bins exist for that cell).
        ``n_valid_bins`` is the total count of bins used.
    """
    a = np.atleast_3d(np.asarray(rate_map1, dtype=float))
    b = np.atleast_3d(np.asarray(rate_map2, dtype=float))
    if a.shape[2] != b.shape[2]:
        raise ValueError("stacked rate maps must match in size along the 3rd (cell) dimension")
    stack_size = a.shape[2]
    max_shape = (max(a.shape[0], b.shape[0]), max(a.shape[1], b.shape[1]))

    interp_a = np.stack(
        [_resize_bilinear(a[:, :, k], max_shape) for k in range(stack_size)], axis=2
    )
    interp_b = np.stack(
        [_resize_bilinear(b[:, :, k], max_shape) for k in range(stack_size)], axis=2
    )

    mask = ((interp_a == 0) & (interp_b == 0)) | np.isnan(interp_a) | np.isnan(interp_b)
    n_valid_bins = int(np.sum(~mask))

    if stack_maps:
        valid = ~mask
        spatial_corr = float(np.corrcoef(interp_a[valid], interp_b[valid])[0, 1])
    else:
        spatial_corr = np.full(stack_size, np.nan)
        for k in range(stack_size):
            valid_k = ~mask[:, :, k]
            if valid_k.sum() >= 2:
                spatial_corr[k] = np.corrcoef(
                    interp_a[:, :, k][valid_k], interp_b[:, :, k][valid_k]
                )[0, 1]

    return spatial_corr, n_valid_bins
