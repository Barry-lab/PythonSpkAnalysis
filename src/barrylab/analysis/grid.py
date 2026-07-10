"""Grid-cell analysis: spatial autocorrelogram (SAC) and gridness.

Ports the lab's grid-cell canon from ``Universal_Matlab/Analysis_Grid`` (and
``xPearson.m`` from ``Analysis_General``):

- :func:`spatial_autocorrelogram` — the 2-D Pearson autocorrelation of a rate
  map (``xPearson.m``), NaN-aware, normalised so the central peak is 1.
- :func:`gridness` — standard (Sargolini 2006) and expanding (Brandon 2011 /
  Killian 2012) gridness from a SAC (``multiGridness.m``).
- :func:`grid_properties` — grid scale, orientation, peak coordinates and mean
  peak r (``sacProps.m``).
- :func:`grid_score` — convenience: rate map -> SAC -> standard gridness.

The MATLAB image-processing calls (``imregionalmax``, ``bwlabel``,
``regionprops``, ``imrotate``) are reproduced with ``scipy.ndimage`` so no
scikit-image dependency is needed. Values are validated structurally against a
known grid cell; to lock them bit-exact, regenerate fixtures with
``validation/generate_fixtures.m`` (see FINDINGS.md).
"""

import numpy as np
from scipy.ndimage import center_of_mass, label, maximum_filter, rotate
from scipy.signal import correlate2d

from .maps import rate_map

_EIGHT = np.ones((3, 3))  # 8-connectivity structuring element
_ROT_ANGLES = (60, 120, 30, 90, 150)


# --------------------------------------------------------------------------
# Spatial autocorrelogram (xPearson.m)
# --------------------------------------------------------------------------
def _filter2(h, x, mode="full"):
    """MATLAB ``filter2(h, x, mode)`` == 2-D correlation of ``h`` over ``x``."""
    return correlate2d(x, h, mode=mode)


def _pearson_all_lags(seq1, seq2, mode="full"):
    """Pearson correlation at every 2-D lag, ignoring NaNs. Port of ``xPearson.m``.

    Returns ``(coef, n_overlap)`` where ``n_overlap`` is the number of finite
    bins contributing to each lag.
    """
    seq1 = np.asarray(seq1, dtype=float)
    seq2 = np.asarray(seq2, dtype=float)
    ones1 = (~np.isnan(seq1)).astype(float)
    ones2 = (~np.isnan(seq2)).astype(float)
    s1 = np.where(np.isnan(seq1), 0.0, seq1)
    s2 = np.where(np.isnan(seq2), 0.0, seq2)

    s1xs2 = _filter2(s1, s2, mode)
    sum1 = _filter2(s1, ones2, mode)
    sum2 = _filter2(ones1, s2, mode)
    sum1sq = _filter2(s1 ** 2, ones2, mode)
    sum2sq = _filter2(ones1, s2 ** 2, mode)
    n = _filter2(ones1, ones2, mode)
    n_sq = n ** 2

    with np.errstate(divide="ignore", invalid="ignore"):
        std1 = np.sqrt(sum1sq / n - sum1 ** 2 / n_sq)
        std2 = np.sqrt(sum2sq / n - sum2 ** 2 / n_sq)
        covar = s1xs2 / n - (sum1 * sum2) / n_sq
        coef = covar / (std1 * std2)
    return np.real(coef), n


def spatial_autocorrelogram(rate_map_2d, min_overlap=None):
    """Spatial autocorrelogram (SAC) of a rate map.

    The SAC is the 2-D Pearson autocorrelation at all spatial lags; its central
    bin (zero lag) is 1 by construction, and a grid cell shows a hexagonal ring
    of six positive peaks around the centre. Output is ``(2*ny-1, 2*nx-1)``.

    Parameters
    ----------
    rate_map_2d : ndarray
        A (smoothed) firing-rate map; NaN bins (unvisited) are ignored.
    min_overlap : int, optional
        If given, SAC bins computed from fewer than this many overlapping map
        bins are set to NaN (Hafting et al. exclude very-low-overlap bins).

    Returns
    -------
    ndarray
        The spatial autocorrelogram.
    """
    coef, n = _pearson_all_lags(rate_map_2d, rate_map_2d, "full")
    if min_overlap is not None:
        coef = np.where(n < min_overlap, np.nan, coef)
    return coef


# --------------------------------------------------------------------------
# Peak finding on the SAC (shared by sacProps.m and multiGridness.m)
# --------------------------------------------------------------------------
def _sac_peaks(sac):
    """Find the central peak and its closest neighbours in a SAC.

    Returns a dict of peak geometry, or ``None`` if no positive peak is found.
    Coordinates follow MATLAB's 1-based, y-down convention for the geometry;
    ``peak_rc`` gives 0-based (row, col) for array indexing.
    """
    sac = np.real(np.asarray(sac, dtype=float))
    filled = np.where(np.isnan(sac), -1.0, sac)
    tmp = np.where(filled <= 0, -1.0, filled)

    # imregionalmax, restricted to positive peaks (the SAC's fields are sharp).
    is_max = (tmp == maximum_filter(tmp, size=3, mode="nearest")) & (tmp > 0)
    lbl, n = label(is_max, structure=_EIGHT)
    if n == 0:
        return None

    com = np.array(center_of_mass(is_max, lbl, range(1, n + 1)))  # (row, col), 0-based
    x1 = com[:, 1] + 1.0  # to MATLAB 1-based
    y1 = com[:, 0] + 1.0
    cy, cx = np.ceil(np.array(sac.shape) / 2.0)  # central bin, 1-based
    xy_central = np.column_stack([x1 - cx, y1 - cy])
    dist = np.sqrt((xy_central ** 2).sum(axis=1))

    order = np.argsort(dist)
    closest = order[:7] if len(order) >= 7 else order
    return {
        "filled": filled,
        "central": (cy, cx),
        "label_ids": (closest + 1),                       # 1-based label ids
        "dist": dist[closest],
        "xy_central": xy_central[closest],                # includes centre as row 0
        "peak_rc": np.rint(com[closest]).astype(int),     # 0-based (row, col)
    }


def _peak_fields(sac_filled, peaks):
    """Label each bin with the id of the peak whose half-height field covers it.

    Port of ``sf_findPeakExtent``: for each peak, take the connected region above
    half the peak height that contains the peak. Overlapping fields resolve to
    the higher id (MATLAB ``max`` over the field stack).
    """
    field = np.zeros(sac_filled.shape)
    for k, pid in enumerate(peaks["label_ids"]):
        row, col = peaks["peak_rc"][k]
        half = sac_filled[row, col] * 0.5
        above, _ = label(sac_filled > half, structure=_EIGHT)
        comp = above[row, col]
        if comp == 0:
            continue
        mask = above == comp
        field[mask] = np.maximum(field[mask], pid)
    return field


def _rotational_gridness(sac_filled, mask):
    """Gridness over a mask: min(r@60,120) - max(r@30,90,150) of rotated SACs."""
    sel_a = sac_filled[mask]
    corrs = []
    for angle in _ROT_ANGLES:
        rot = rotate(sac_filled, angle, reshape=False, order=1, mode="constant", cval=0.0)
        sel_b = rot[mask]
        valid = np.isfinite(sel_a) & np.isfinite(sel_b)
        if valid.sum() < 2:
            corrs.append(np.nan)
        else:
            corrs.append(np.corrcoef(sel_a[valid], sel_b[valid])[0, 1])
    return min(corrs[0], corrs[1]) - max(corrs[2], corrs[3], corrs[4])


# --------------------------------------------------------------------------
# Gridness (multiGridness.m)
# --------------------------------------------------------------------------
def gridness(sac):
    """Standard and expanding gridness of a SAC. Port of ``multiGridness.m``.

    Parameters
    ----------
    sac : ndarray
        A spatial autocorrelogram (from :func:`spatial_autocorrelogram`).

    Returns
    -------
    (standard, expanding)
        ``standard`` — Sargolini (2006) gridness over a disc out to the outer
        edge of the six central fields. ``expanding`` — the maximum gridness
        over a family of expanding annuli (Brandon 2011 / Killian 2012). The two
        are not directly comparable in absolute value. ``(nan, nan)`` if fewer
        than two peaks are found.
    """
    peaks = _sac_peaks(sac)
    if peaks is None or len(peaks["label_ids"]) < 2:
        return float("nan"), float("nan")

    filled = peaks["filled"]
    field = _peak_fields(filled, peaks)
    cy, cx = peaks["central"]
    yy, xx = np.mgrid[1:filled.shape[0] + 1, 1:filled.shape[1] + 1]
    centre_radi = np.sqrt((xx - cx) ** 2 + (yy - cy) ** 2)
    max_allow = int(np.min(np.floor(np.array(filled.shape) / 2.0)))
    central_id = peaks["label_ids"][0]

    # Standard: disc out to the furthest edge of the six fields, minus the centre.
    outer = 0.0
    for pid in peaks["label_ids"]:
        m = field == pid
        if m.any():
            outer = max(outer, centre_radi[m].max())
    outer = min(outer, max_allow)
    standard_mask = (centre_radi <= outer) & (field != central_id)
    standard = _rotational_gridness(filled, standard_mask)

    # Expanding: max gridness over annuli from inner+2 out to the SAC edge.
    inner = int(round(np.mean(np.sqrt((peaks["xy_central"] ** 2).sum(axis=0))) / 2.0))
    best = -np.inf
    for outer_r in range(inner + 2, max_allow + 1):
        annulus = (centre_radi <= outer_r) & (centre_radi >= inner)
        g = _rotational_gridness(filled, annulus)
        if np.isfinite(g):
            best = max(best, g)
    expanding = best if np.isfinite(best) else float("nan")
    return standard, expanding


def grid_properties(sac):
    """Grid scale, orientation and peak geometry from a SAC. Port of ``sacProps.m``.

    Returns
    -------
    dict with keys:
        ``scale`` — grid scale in SAC bins (median distance of the 6 central
            peaks from the centre),
        ``gridness`` — standard gridness,
        ``orientation`` — sorted axis orientations (deg CCW from +x); the first
            is the grid's overall orientation,
        ``peak_coord`` — (n, 2) peak positions as (x, y), origin top-left,
        ``mean_r`` — mean SAC value at the six central peaks.
    """
    peaks = _sac_peaks(sac)
    if peaks is None or len(peaks["label_ids"]) < 2:
        return {"scale": np.nan, "gridness": np.nan, "orientation": np.array([]),
                "peak_coord": np.empty((0, 2)), "mean_r": np.nan}

    filled = peaks["filled"]
    scale = float(np.median(peaks["dist"][1:]))  # exclude central peak

    # Orientation: polar angle of peaks above the x-axis (y negated for image frame).
    xy = peaks["xy_central"]
    theta = np.arctan2(-xy[:, 1], xy[:, 0])
    keep = (theta >= 0)
    keep[0] = False  # drop the central point
    orientation = np.degrees(np.sort(theta[keep]))

    # Mean r at the peaks.
    mean_r = float(np.mean([filled[r, c] for r, c in peaks["peak_rc"]]))

    standard, _ = gridness(sac)
    # (x, y) with origin top-left = 0-based col, row + 1 back to MATLAB 1-based.
    peak_coord = np.column_stack([peaks["peak_rc"][:, 1] + 1, peaks["peak_rc"][:, 0] + 1])
    return {"scale": scale, "gridness": standard, "orientation": orientation,
            "peak_coord": peak_coord, "mean_r": mean_r}


def grid_score(recording, electrode, cluster, bin_cm=2.0, smoothing="boxcar",
               smooth_size=None, min_overlap=None):
    """Convenience: rate map -> SAC -> standard gridness for one cluster.

    Returns
    -------
    (standard_gridness, scale_cm)
        Standard gridness and grid scale converted to cm (bins * ``bin_cm``).
    """
    rm = rate_map(recording, electrode, cluster, bin_cm=bin_cm,
                  smoothing=smoothing, smooth_size=smooth_size)
    sac = spatial_autocorrelogram(rm, min_overlap=min_overlap)
    props = grid_properties(sac)
    return props["gridness"], props["scale"] * bin_cm
