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
- :func:`fast_gridness` — a leaner, DISTINCT standard/expanding-gridness
  implementation used inside the FieldShuffle significance loop
  (``fs_fastGridness.m``, Langston et al. 2011; see
  :mod:`barrylab.analysis.field_shuffle`). Not a duplicate of
  :func:`gridness`/:func:`grid_properties` — see its docstring for several
  confirmed algorithmic divergences from ``multiGridness.m``, faithfully
  preserved rather than reconciled.
- :func:`square_grid_props` — a 4-fold/square-symmetric analogue of
  :func:`grid_properties`/:func:`gridness` for non-hexagonal grids
  (``sacProps90deg.m``).

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
_ROT_ANGLES_90 = (90, 45, 135)  # square/4-fold-symmetric grid variant (sacProps90deg.m)


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
    sum1sq = _filter2(s1**2, ones2, mode)
    sum2sq = _filter2(ones1, s2**2, mode)
    n = _filter2(ones1, ones2, mode)
    n_sq = n**2

    with np.errstate(divide="ignore", invalid="ignore"):
        std1 = np.sqrt(sum1sq / n - sum1**2 / n_sq)
        std2 = np.sqrt(sum2sq / n - sum2**2 / n_sq)
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


def cross_correlogram(rate_map_a, rate_map_b, min_overlap=None):
    """Cross spatial-correlogram between two rate maps (``xPearson.m``, cross mode).

    Like :func:`spatial_autocorrelogram` but between two distinct maps —
    typically two cells from the same grid module, whose cross-correlogram's
    central-peak offset from the geometric centre encodes their relative
    spatial phase (see :func:`~barrylab.analysis.grid_phase.cross_sac_props`).
    Unlike an autocorrelogram, there is no guarantee the peak nearest the
    geometric centre sits exactly at it, and central-symmetry does not apply.

    Parameters
    ----------
    rate_map_a, rate_map_b : ndarray
        (Smoothed) firing-rate maps; NaN bins (unvisited) are ignored.
    min_overlap : int, optional
        If given, bins computed from fewer than this many overlapping map
        bins are set to NaN.

    Returns
    -------
    ndarray
        The cross spatial-correlogram, shape ``(ny_a+ny_b-1, nx_a+nx_b-1)``.
    """
    coef, n = _pearson_all_lags(rate_map_a, rate_map_b, "full")
    if min_overlap is not None:
        coef = np.where(n < min_overlap, np.nan, coef)
    return coef


# --------------------------------------------------------------------------
# Peak finding on the SAC (shared by sacProps.m and multiGridness.m)
# --------------------------------------------------------------------------
def _sac_peaks(sac, n_closest=7):
    """Find the central peak and its closest neighbours in a SAC.

    ``n_closest`` is the number of nearest-to-centre peaks to keep (7 =
    centre + 6 hexagonal neighbours, matching ``sacProps.m``/``multiGridness.m``;
    :func:`square_grid_props` passes 5 = centre + 4 for a square/4-fold grid,
    matching ``sacProps90deg.m``).

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
    dist = np.sqrt((xy_central**2).sum(axis=1))

    order = np.argsort(dist)
    closest = order[:n_closest] if len(order) >= n_closest else order
    return {
        "filled": filled,
        "central": (cy, cx),
        "label_ids": (closest + 1),  # 1-based label ids
        "dist": dist[closest],
        "xy_central": xy_central[closest],  # includes centre as row 0
        "peak_rc": np.rint(com[closest]).astype(int),  # 0-based (row, col)
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
    yy, xx = np.mgrid[1 : filled.shape[0] + 1, 1 : filled.shape[1] + 1]
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
    # Inner radius = half the mean distance to the central peaks. sum(axis=1)
    # is per-peak (each row is one peak's [x, y] offset) -- axis=0 would sum
    # down columns (all peaks' x^2 then all y^2), which was the multiGridness.m
    # bug fixed upstream (Universal_Matlab python-port-bugfixes commit bc4feb9).
    # MATLAB round() is half-away-from-zero; the value is >= 0 so floor(x+0.5).
    inner = int(np.floor(np.mean(np.sqrt((peaks["xy_central"] ** 2).sum(axis=1))) / 2.0 + 0.5))
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
        return {
            "scale": np.nan,
            "gridness": np.nan,
            "orientation": np.array([]),
            "peak_coord": np.empty((0, 2)),
            "mean_r": np.nan,
        }

    filled = peaks["filled"]
    scale = float(np.median(peaks["dist"][1:]))  # exclude central peak

    # Orientation: polar angle of peaks above the x-axis (y negated for image frame).
    xy = peaks["xy_central"]
    theta = np.arctan2(-xy[:, 1], xy[:, 0])
    keep = theta >= 0
    keep[0] = False  # drop the central point
    orientation = np.degrees(np.sort(theta[keep]))

    # Mean r at the peaks.
    mean_r = float(np.mean([filled[r, c] for r, c in peaks["peak_rc"]]))

    standard, _ = gridness(sac)
    # (x, y) with origin top-left = 0-based col, row + 1 back to MATLAB 1-based.
    peak_coord = np.column_stack([peaks["peak_rc"][:, 1] + 1, peaks["peak_rc"][:, 0] + 1])
    return {
        "scale": scale,
        "gridness": standard,
        "orientation": orientation,
        "peak_coord": peak_coord,
        "mean_r": mean_r,
    }


def grid_score(
    recording,
    electrode,
    cluster,
    bin_cm=2.0,
    smoothing="boxcar",
    smooth_size=None,
    min_overlap=None,
):
    """Convenience: rate map -> SAC -> standard gridness for one cluster.

    Returns
    -------
    (standard_gridness, scale_cm)
        Standard gridness and grid scale converted to cm (bins * ``bin_cm``).
    """
    rm = rate_map(
        recording, electrode, cluster, bin_cm=bin_cm, smoothing=smoothing, smooth_size=smooth_size
    )
    sac = spatial_autocorrelogram(rm, min_overlap=min_overlap)
    props = grid_properties(sac)
    return props["gridness"], props["scale"] * bin_cm


# --------------------------------------------------------------------------
# Square/4-fold-symmetric grid variant (sacProps90deg.m)
# --------------------------------------------------------------------------
def _rotational_gridness_90(sac_filled, mask):
    """90-degree-symmetry gridness over a mask: corr(90) - max(corr(45), corr(135)).

    The square-grid analogue of :func:`_rotational_gridness`'s 60/120-vs-
    30/90/150 hexagonal split.
    """
    sel_a = sac_filled[mask]
    corrs = []
    for angle in _ROT_ANGLES_90:
        rot = rotate(sac_filled, angle, reshape=False, order=1, mode="constant", cval=0.0)
        sel_b = rot[mask]
        valid = np.isfinite(sel_a) & np.isfinite(sel_b)
        if valid.sum() < 2:
            corrs.append(np.nan)
        else:
            corrs.append(np.corrcoef(sel_a[valid], sel_b[valid])[0, 1])
    # MATLAB: gridness = min(rotatedCorr(1)) - max(rotatedCorr([2,3])) --
    # min() of the single 90-degree value is a no-op; np.fmax matches
    # MATLAB max()'s NaN-ignoring semantics for the 45/135 pair (the same
    # recurring gotcha handled throughout this project, e.g. _brandon_gridness).
    return corrs[0] - np.fmax(corrs[1], corrs[2])


def square_grid_props(sac):
    """Square/4-fold-symmetric grid metrics from a SAC. Port of ``sacProps90deg.m``.

    A variant of :func:`grid_properties`/:func:`gridness` for grids with
    4-fold (square -- e.g. some trained-network grid cells) rather than
    6-fold (hexagonal) symmetry: looks for the 4 peaks closest to the SAC
    centre (not 6), and measures gridness as the correlation at a 90-degree
    rotation minus the maximum correlation at 45/135 degrees (instead of the
    usual min(60,120) minus max(30,90,150) hexagonal split).

    Uses the *same* peak-detection convention as :func:`grid_properties`/
    :func:`gridness` (``multiGridness.m``/``sacProps.m``: ``imregionalmax``
    restricted to positive SAC bins, one candidate per connected plateau) --
    confirmed by direct inspection of ``sacProps90deg.m``, NOT the
    "raw"/unrestricted convention used by :func:`fast_gridness`/
    :func:`~barrylab.analysis.grid_phase.cross_sac_props`/
    :func:`~barrylab.analysis.grid_ellipse.brandon_regularise`.

    Parameters
    ----------
    sac : ndarray
        A spatial autocorrelogram (from :func:`spatial_autocorrelogram`).

    Returns
    -------
    dict with keys:
        ``scale`` -- array (up to 5,): mean distance of the up-to-4 central
            peaks from the centre (element 0), then each peak's own distance
            (``[nan]`` if only the central peak was found).
        ``gridness`` -- scalar 90-degree gridness (``nan`` if fewer than 2
            peaks found).
        ``orientation`` -- sorted axis orientations (deg CCW from +x) of
            peaks above the x-axis (``[nan]`` if only the central peak was
            found).
        ``peak_coord`` -- (n, 2) peak positions as (x, y), origin top-left.
    """
    peaks = _sac_peaks(sac, n_closest=5)
    if peaks is None:
        return {
            "scale": np.array([np.nan]),
            "gridness": float("nan"),
            "orientation": np.array([np.nan]),
            "peak_coord": np.empty((0, 2)),
        }

    filled = peaks["filled"]
    dist = peaks["dist"]
    if dist.size > 1:
        scale = np.concatenate(([float(np.mean(dist[1:]))], dist[1:]))
    else:
        scale = np.array([np.nan])

    if len(peaks["label_ids"]) == 1:  # only the central peak found (numPeaks==1)
        orientation = np.array([np.nan])
    else:
        xy = peaks["xy_central"]
        theta = np.arctan2(-xy[:, 1], xy[:, 0])  # y negated: image y-down frame
        keep = theta >= 0
        keep[0] = False  # drop the central point
        orientation = np.degrees(np.sort(theta[keep]))

    if len(peaks["label_ids"]) < 2:
        gridness90 = float("nan")
    else:
        field = _peak_fields(filled, peaks)
        cy, cx = peaks["central"]
        yy, xx = np.mgrid[1 : filled.shape[0] + 1, 1 : filled.shape[1] + 1]
        centre_radi = np.sqrt((xx - cx) ** 2 + (yy - cy) ** 2)
        max_allow = int(np.min(np.floor(np.array(filled.shape) / 2.0)))
        central_id = peaks["label_ids"][0]

        # Disc mask out to the furthest edge of the (up-to-4) fields' extent
        # -- identical construction to gridness()'s "standard" mask.
        outer = 0.0
        for pid in peaks["label_ids"]:
            m = field == pid
            if m.any():
                outer = max(outer, centre_radi[m].max())
        outer = min(outer, max_allow)
        mask = (centre_radi <= outer) & (field != central_id)
        gridness90 = _rotational_gridness_90(filled, mask)

    peak_coord = np.column_stack([peaks["peak_rc"][:, 1] + 1, peaks["peak_rc"][:, 0] + 1])
    return {
        "scale": scale,
        "gridness": gridness90,
        "orientation": orientation,
        "peak_coord": peak_coord,
    }


# --------------------------------------------------------------------------
# Fast gridness for the FieldShuffle loop (fs_fastGridness.m)
# --------------------------------------------------------------------------
def _regional_minima_1d(x):
    """Positions of ``imregionalmin()``-style regional minima in a 1-D array.

    A regional minimum is a maximal flat run strictly lower than the values
    immediately bordering it on both sides (an array edge only needs the one
    existing neighbour to be higher).
    """
    x = np.asarray(x, dtype=float)
    n = x.size
    is_min = np.zeros(n, dtype=bool)
    i = 0
    while i < n:
        j = i
        while j + 1 < n and x[j + 1] == x[i]:
            j += 1
        left_ok = (i == 0) or (x[i - 1] > x[i])
        right_ok = (j == n - 1) or (x[j + 1] > x[j])
        if left_ok and right_ok:
            is_min[i : j + 1] = True
        i = j + 1
    return is_min


def _raw_regional_max_peaks(filled, positive_only=False):
    """0-based ``(row, col)`` of every ``imregionalmax`` peak pixel.

    Shared by :func:`fast_gridness`,
    :func:`~barrylab.analysis.grid_phase.cross_sac_props` and
    :func:`~barrylab.analysis.grid_ellipse.brandon_regularise`. Peaks are
    returned in MATLAB's column-major ``find()`` order, since that order
    matters for tie-breaking wherever a caller later truncates to the N
    closest peaks.

    ``positive_only`` restricts peaks to strictly-positive bins (the
    convention ``_sac_peaks`` uses for :func:`gridness`/:func:`grid_properties`
    via ``multiGridness.m``). ``fs_fastGridness.m`` and ``xSacProps.m`` were
    both fixed upstream to do this (Universal_Matlab python-port-bugfixes
    commits 7c256ae and a443363) -- without it, a large flat region of
    negative/NaN-substituted bins is treated as locally maximal everywhere,
    flooding the candidate list with spurious tied "peaks". ``brandonRegEcc.m``
    was not changed upstream, so :func:`brandon_regularise` keeps the raw
    (``positive_only=False``) behaviour.
    """
    is_max = filled == maximum_filter(filled, size=3, mode="nearest")
    if positive_only:
        is_max = is_max & (filled > 0)
    flat_idx = np.flatnonzero(is_max.ravel(order="F"))
    return np.unravel_index(flat_idx, is_max.shape, order="F")


def fast_gridness(sac):
    """Fast standard + expanding (Langston 2011) gridness, and grid scale.

    Port of ``fs_fastGridness.m`` — used inside the (not yet ported)
    FieldShuffle significance-testing loop, where it is called many times per
    cell and so trades a little precision for speed relative to
    :func:`gridness`/:func:`grid_properties`.

    .. warning::
       Despite the similar name and purpose, this is **not** numerically
       identical to :func:`gridness`/:func:`grid_properties`, even for
       standard gridness on the same SAC. Confirmed, faithfully-preserved
       divergences from ``multiGridness.m`` (see ``FINDINGS.md``):

       1. Peak detection here does not restrict ``imregionalmax`` to
          positive-valued SAC bins the way ``multiGridness.m`` does (it masks
          non-positive bins to -1 first) — so spurious peaks among negative
          or NaN-substituted (-1) SAC regions can in principle be picked up.
       2. It builds the candidate peak list directly from ``find()`` on the
          regional-max mask -- one candidate per *pixel* of any local-max
          plateau -- rather than ``bwlabel``+``regionprops`` centroids (one
          candidate per connected plateau).
       3. All radius/mask comparisons here use a distance-from-centre grid
          **rounded** to the nearest integer bin (needed for its precomputed
          disc-mask family); ``multiGridness.m`` keeps this distance
          continuous throughout.
       4. Grid scale is the **mean** distance to the closest ~6 peaks;
          :func:`grid_properties` uses the **median**.
       5. Expanding gridness here picks its central/inner radius from the
          SAC's own radial mean-correlation profile (first regional minimum
          or first negative crossing, Langston et al. 2011) rather than half
          the mean peak distance (Killian 2012, used by :func:`gridness`).

    Parameters
    ----------
    sac : ndarray
        A spatial autocorrelogram (from :func:`spatial_autocorrelogram`).

    Returns
    -------
    (std_g, exp_g, grid_scale)
        ``std_g``/``exp_g`` — standard / expanding gridness (``nan`` if fewer
        than 2 peaks are found). ``grid_scale`` — mean distance (in SAC bins)
        to the closest peaks, excluding the centre (``nan`` if none).
    """
    exp_g_min, exp_g_step, exp_g_max = 5, 1, 5

    sac = np.real(np.asarray(sac, dtype=float))
    filled = np.where(np.isnan(sac), -1.0, sac)
    n_rows, n_cols = filled.shape

    ys, xs = _raw_regional_max_peaks(filled, positive_only=True)
    if ys.size < 2:
        return float("nan"), float("nan"), float("nan")

    central_row = int(np.ceil(n_rows / 2.0))  # 1-based
    central_col = int(np.ceil(n_cols / 2.0))
    dx = (xs + 1.0) - central_col
    dy = (ys + 1.0) - central_row
    dist = np.sqrt(dx**2 + dy**2)

    order = np.argsort(dist, kind="stable")
    n_take = min(7, order.size)
    closest = order[:n_take]
    dist_from_cent = dist[closest]
    grid_scale = float(np.mean(dist_from_cent[1:])) if n_take > 1 else float("nan")

    peak_rc = np.column_stack([ys[closest], xs[closest]])  # 0-based (row, col)
    # Field label ids are the peaks' ORIGINAL (pre-distance-sort) find()-order
    # positions -- matching MATLAB's `closestPeaks` values exactly, not a
    # fresh 1..7 relabelling -- because _peak_fields resolves overlapping
    # half-height fields by keeping the numerically larger id, and that
    # tie-break must match MATLAB's (effectively arbitrary) id assignment.
    label_ids = closest + 1
    field = _peak_fields(filled, {"label_ids": label_ids, "peak_rc": peak_rc})

    # Rounded (integer-bin) distance-from-centre grid, used for every
    # radius/mask comparison below -- see divergence (3) in the docstring.
    yy, xx = np.mgrid[1 : n_rows + 1, 1 : n_cols + 1]
    dist_grid = np.floor(np.sqrt((xx - central_col) ** 2 + (yy - central_row) ** 2) + 0.5)
    half_min_sz = int(np.min(np.floor(np.array([n_rows, n_cols]) / 2.0)))

    # --- Standard gridness -------------------------------------------------
    field_mask = field != 0
    max_dist = dist_grid[field_mask].max() if field_mask.any() else 0.0
    max_dist = min(max_dist, half_min_sz)
    std_mask = (dist_grid <= max_dist) & (field != label_ids[0])
    std_g = _rotational_gridness(filled, std_mask)

    # --- Expanding (Langston 2011) gridness ---------------------------------
    # Radial mean-correlation profile: mean SAC value at each integer radius.
    radii = dist_grid[dist_grid > 0].astype(np.int64)
    vals = filled[dist_grid > 0]
    max_r = int(radii.max()) if radii.size else 0
    corr_vs_rad = np.zeros(max_r)
    if max_r:
        sums = np.bincount(radii - 1, weights=vals, minlength=max_r)
        counts = np.bincount(radii - 1, minlength=max_r)
        with np.errstate(invalid="ignore", divide="ignore"):
            corr_vs_rad = np.where(counts > 0, sums / np.where(counts > 0, counts, 1), 0.0)
    corr_vs_rad = np.where(np.isnan(corr_vs_rad), 0.0, corr_vs_rad)

    rad_of_cent = None
    if corr_vs_rad.size:
        candidates = []
        min_idx = np.flatnonzero(_regional_minima_1d(corr_vs_rad))
        if min_idx.size:
            candidates.append(min_idx[0] + 1)  # 1-based radius
        neg_idx = np.flatnonzero(corr_vs_rad < 0)
        if neg_idx.size:
            candidates.append(neg_idx[0] + 1)
        if candidates:
            rad_of_cent = min(candidates)

    def _candidate_radii(r0):
        return np.arange(r0 + exp_g_min, half_min_sz - exp_g_max + 1, exp_g_step)

    mask_ext_rad = _candidate_radii(rad_of_cent) if rad_of_cent is not None else np.array([])
    if mask_ext_rad.size == 0:
        rad_of_cent = 5
        mask_ext_rad = _candidate_radii(rad_of_cent)

    exp_g = float("nan")
    if mask_ext_rad.size:
        exp_gs = []
        for r in mask_ext_rad:
            annulus = (dist_grid <= r) & (dist_grid > rad_of_cent)
            g = _rotational_gridness(filled, annulus)
            if np.isfinite(g):
                exp_gs.append(g)
        if exp_gs:
            exp_g = max(exp_gs)

    return std_g, exp_g, grid_scale
