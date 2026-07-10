"""Elliptical grid correction: fit and remove elliptical distortion from a SAC.

Ports ``gridEllipse_fit.m`` (fit an ellipse to a grid cell's six central SAC
peaks, giving x/y scale and eccentricity), ``gridEllipse_correct.m``
(regularise an elliptical SAC back toward circular, so gridness can be
recomputed on the corrected map), and ``brandonRegEcc.m`` (an alternative,
simpler ellipse estimate from just the closest/furthest of the six central
peaks). ``regularise_eliptic_grid.m`` is a byte-identical duplicate of
``gridEllipse_correct.m`` in the original — ported once here.

The ellipse-fit itself (``_fit_ellipse``) is a least-squares conic fit
originally from the MathWorks File Exchange, embedded in ``gridEllipse_fit.m``
as a subfunction; ported faithfully for portability, as the original was.
"""

import numpy as np
from scipy.ndimage import map_coordinates, rotate

from .grid import _raw_regional_max_peaks, _sac_peaks


def _fit_ellipse(x, y):
    """Least-squares ellipse fit (MathWorks File Exchange ``fit_ellipse``).

    Fits ``a*x^2 + b*x*y + c*y^2 + d*x + e*y + f = 0`` by least squares, then
    extracts the standard ellipse parameters. Returns ``None`` if the points
    describe a parabola or hyperbola rather than an ellipse.

    Returns
    -------
    dict with ``a`` (major semi-axis), ``b`` (minor semi-axis), ``phi``
    (tilt, radians), ``x0``/``y0`` (centre, in the de-meaned frame).
    """
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    mean_x, mean_y = x.mean(), y.mean()
    x = x - mean_x
    y = y - mean_y

    design = np.column_stack([x**2, x * y, y**2, x, y])
    coeffs, *_ = np.linalg.lstsq(design.T @ design, design.T @ np.ones_like(x), rcond=None)
    a, b, c, d, e = coeffs

    orientation_tol = 1e-3
    if min(abs(b / a) if a else np.inf, abs(b / c) if c else np.inf) > orientation_tol:
        orientation_rad = 0.5 * np.arctan2(b, c - a)
        cos_phi, sin_phi = np.cos(orientation_rad), np.sin(orientation_rad)
        a, b, c, d, e = (
            a * cos_phi**2 - b * cos_phi * sin_phi + c * sin_phi**2,
            0.0,
            a * sin_phi**2 + b * cos_phi * sin_phi + c * cos_phi**2,
            d * cos_phi - e * sin_phi,
            d * sin_phi + e * cos_phi,
        )
        mean_x, mean_y = cos_phi * mean_x - sin_phi * mean_y, sin_phi * mean_x + cos_phi * mean_y
    else:
        orientation_rad = 0.0
        cos_phi, sin_phi = 1.0, 0.0

    if a * c <= 0:
        return None  # parabola or hyperbola, not an ellipse

    if a < 0:
        a, c, d, e = -a, -c, -d, -e

    x0 = mean_x - d / (2 * a)
    y0 = mean_y - e / (2 * c)
    f = 1 + (d**2) / (4 * a) + (e**2) / (4 * c)
    semi_a, semi_b = np.sqrt(f / a), np.sqrt(f / c)

    return {"a": semi_a, "b": semi_b, "phi": orientation_rad, "x0": x0, "y0": y0}


def fit_grid_ellipse(sac):
    """Fit an ellipse to a grid cell's six central SAC peaks.

    Port of ``gridEllipse_fit.m``. Requires at least 6 detected peaks around
    the centre; returns all-NaN if fewer are found or no ellipse fits (a
    parabola/hyperbola).

    Parameters
    ----------
    sac : ndarray
        Spatial autocorrelogram, ideally built from a smoothed rate map.

    Returns
    -------
    dict with:
        ``xy_scale`` — ``(x_scale, y_scale)`` in SAC bins: where the fitted
            ellipse crosses each axis.
        ``eccentricity`` — 0 for a circle.
        ``orientation`` — major-axis orientation, radians CCW from +x
            (image/ij convention).
        ``ab_scale`` — ``(a, b)`` major/minor semi-axis lengths, in bins.
    """
    peaks = _sac_peaks(sac)
    nan_result = {
        "xy_scale": (np.nan, np.nan),
        "eccentricity": np.nan,
        "orientation": np.nan,
        "ab_scale": (np.nan, np.nan),
    }
    if peaks is None or len(peaks["label_ids"]) < 7:
        return nan_result

    # peak_rc[0] is the central peak; drop it, keep the 6 closest neighbours.
    peak_rc = peaks["peak_rc"][1:7]
    px = peak_rc[:, 1] + 1.0  # 1-based x (col)
    py = peak_rc[:, 0] + 1.0  # 1-based y (row)

    fit = _fit_ellipse(px, py)
    if fit is None:
        return nan_result

    a, b, orient = fit["a"], fit["b"], fit["phi"]
    if a < b:
        a, b = b, a
        orient = np.mod(orient + np.pi / 2, 2 * np.pi)

    eccentricity = np.sqrt((a**2 - b**2) / a**2)

    theta_axes = np.array([orient, np.pi / 2 + orient])
    ellipse_x = a * np.cos(theta_axes)
    ellipse_y = b * np.sin(theta_axes)
    xy_scale = np.sqrt(ellipse_x**2 + ellipse_y**2)

    return {
        "xy_scale": tuple(xy_scale),
        "eccentricity": float(eccentricity),
        "orientation": float(orient),
        "ab_scale": (float(a), float(b)),
    }


def regularise_elliptic_sac(sac, ab_scale, orient_deg):
    """Reshape an elliptical SAC to be circular. Port of ``gridEllipse_correct.m``.

    Rotates the SAC so its major axis aligns with the x-axis, rescales the
    (now-)minor axis to match the major axis, then rotates back — regularising
    the grid so gridness can be recomputed on a corrected map. Also serves as
    ``regularise_eliptic_grid.m``, a byte-identical duplicate in the original.

    .. note::
       ``orient_deg`` is in **degrees**, unlike :func:`fit_grid_ellipse`'s
       ``orientation`` output, which is in **radians** — this mismatch exists
       in the original MATLAB too (``gridEllipse_fit.m`` returns radians,
       ``gridEllipse_correct.m`` expects degrees; callers such as
       ``brandonRegEcc.m`` explicitly `rad2deg`-convert between them).
       Convert with ``np.degrees(...)`` when chaining the two functions.

    Parameters
    ----------
    sac : ndarray
    ab_scale : (float, float)
        Major/minor semi-axis lengths (order doesn't matter — re-sorted).
    orient_deg : float
        Major-axis orientation in degrees, CCW from +x.

    Returns
    -------
    ndarray, same shape as ``sac``.
    """
    sac = np.asarray(sac, dtype=float)
    a, b = sorted(ab_scale, reverse=True)
    if a == b:
        return sac.copy()

    # Matches the rotation-angle convention already validated in grid.py's
    # gridness rotations (angle passed through to scipy.ndimage.rotate
    # unchanged); MATLAB instruction here is imrotate(sac, -orient, ...).
    reg = rotate(sac, -orient_deg, reshape=True, order=1, mode="constant", cval=0.0)

    n_rows, n_cols = reg.shape
    shrink = (n_rows / 2.0) - (n_rows / 2.0) * (b / a)
    y_start, y_end = shrink, n_rows - 1 - shrink  # 0-based (MATLAB: 1+tmp, sacSize(1)-tmp)
    x_coords = np.linspace(0, n_cols - 1, n_cols)
    y_coords = np.linspace(y_start, y_end, n_rows)
    yy, xx = np.meshgrid(y_coords, x_coords, indexing="ij")
    reg = map_coordinates(reg, [yy, xx], order=1, mode="constant", cval=np.nan)

    reg = rotate(reg, orient_deg, reshape=True, order=1, mode="constant", cval=np.nan)

    # Centre-crop back to the original size. MATLAB's own equivalent
    # (`sizeDif=round(...)`) can be off by one pixel for an odd size
    # difference (MATLAB's round-half-away-from-zero vs numpy's
    # round-half-to-even -- see FINDINGS.md); crop directly to the exact
    # target shape instead, matching the code's own stated intent ("keep
    # only the central portion... so that it matches the original size").
    start_row = (reg.shape[0] - sac.shape[0]) // 2
    start_col = (reg.shape[1] - sac.shape[1]) // 2
    reg = reg[start_row : start_row + sac.shape[0], start_col : start_col + sac.shape[1]]
    return reg


def brandon_regularise(sac):
    """Regularise an elliptical SAC via the Brandon (2011) method.

    Port of ``brandonRegEcc.m``. An alternative to :func:`fit_grid_ellipse` +
    :func:`regularise_elliptic_sac`'s least-squares approach: estimates the
    ellipse's major/minor axes from just the closest and furthest of the six
    central peaks, in two different ways, then regularises the SAC with each
    estimate independently.

    .. note::
       The MATLAB docstring calls this method "a bit clunky" and says it
       "doesn't actually work that well -- certainly not as well as" the
       least-squares fit, kept mainly because the least-squares fit can fail
       (multiple minima) in some cases.

       Peak detection here follows ``fs_fastGridness.m``/``xSacProps.m``'s
       convention (``imregionalmax`` with no positive-value restriction), NOT
       :func:`fit_grid_ellipse`/:func:`gridness`'s ``multiGridness.m``
       convention -- confirmed by reading ``brandonRegEcc.m`` directly, not
       assumed from the similar name. See ``FINDINGS.md``.

       The MATLAB source's own comment flags that its first axis-length
       formula "requires a sqrt of a negative value... I've just ignored the
       negative sign" -- i.e. it can legitimately go complex, and MATLAB's
       ``abs()`` of a complex number returns its magnitude. Faithfully
       reproduced with ``numpy.lib.scimath.sqrt`` (which, like MATLAB's
       ``sqrt``, promotes a negative real argument to a complex result
       rather than ``nan``). Method 1's axis lengths are wrapped in
       ``abs()`` before use, matching the MATLAB call
       (``abs([m1MajR, m1MinR])``) -- but **method 2's are not** (the MATLAB
       call is the bare ``[m2MajR, m2MinR]``), which is asymmetric with
       method 1's treatment of the identical failure mode and looks like an
       oversight in the original rather than a deliberate choice. Faithfully
       preserved: for method 2, if the axis length comes out complex, its
       real part is used, since :func:`regularise_elliptic_sac` sorts/compares
       its ``ab_scale`` argument (undefined for complex numbers in Python,
       unlike MATLAB's real-part-based complex comparison) -- this only
       changes what happens in an already-degenerate corner the original
       author's own comment says is not well understood.

    Parameters
    ----------
    sac : ndarray
        Spatial autocorrelogram (from a smoothed rate map).

    Returns
    -------
    (rsac1, rsac2)
        Two independently-regularised SACs (method 1: minor axis from the
        closest peak; method 2: major axis from the furthest peak). Returns
        the original ``sac`` unchanged (twice) if fewer than 4 peaks are
        found, or per-method if that method's major/minor axis ratio exceeds
        2x (Brandon 2011's own no-correction criterion).
    """
    sac = np.real(np.asarray(sac, dtype=float))
    filled = np.where(np.isnan(sac), -1.0, sac)

    peak_m, peak_n = _raw_regional_max_peaks(filled)  # 0-based
    n_rows, n_cols = filled.shape
    central_row = int(np.ceil(n_rows / 2.0))  # 1-based
    central_col = int(np.ceil(n_cols / 2.0))
    dx = (peak_n + 1.0) - central_col
    dy = (peak_m + 1.0) - central_row
    # "Not euclidian" (MATLAB's own comment) -- squared, unrooted distance,
    # used only to rank/sort peaks so the monotonic transform doesn't matter.
    dist_sq = dx**2 + dy**2
    order = np.argsort(dist_sq, kind="stable")
    # closestPeaks(2:min(7,length)) -- drop index 0, the central peak itself.
    closest = order[1 : min(7, order.size)]

    if closest.size < 4:
        return sac.copy(), sac.copy()

    # x,y pairs (1-based), y counting down from top, origin top-left.
    px = peak_n[closest] + 1.0
    py = peak_m[closest] + 1.0
    ang = np.arctan2(-py, px)
    r = np.hypot(px, py)
    close_ang, close_r = ang[0], r[0]
    dist_ang, dist_r = ang[-1], r[-1]

    theta = np.mod(close_ang - dist_ang, 2 * np.pi)
    cos_sq_th = np.cos(theta) ** 2
    sin_sq_th = np.sin(theta) ** 2
    csqrt = np.lib.scimath.sqrt

    # Method 1: minor axis set to the closest peak; estimate major axis.
    m1_min_r = close_r
    m1_maj_th = np.mod(close_ang + np.pi / 2, 2 * np.pi)
    m1_maj_r = csqrt(cos_sq_th / ((sin_sq_th / close_r**2) - (dist_r**-2)))

    # Method 2: major axis set to the furthest peak; estimate minor axis.
    m2_maj_r = dist_r
    m2_maj_th = np.mod(dist_ang, 2 * np.pi)
    m2_min_r = csqrt(sin_sq_th / ((close_r**-2) - (cos_sq_th / dist_r**2)))

    if abs(m1_maj_r) > 2 * abs(m1_min_r):
        rsac1 = sac.copy()
    else:
        rsac1 = regularise_elliptic_sac(
            sac, (abs(m1_maj_r), abs(m1_min_r)), float(np.degrees(m1_maj_th))
        )

    m2_min_r_real = m2_min_r.real if np.iscomplexobj(m2_min_r) else m2_min_r
    if m2_maj_r > 2 * m2_min_r_real:
        rsac2 = sac.copy()
    else:
        rsac2 = regularise_elliptic_sac(
            sac, (float(m2_maj_r), float(m2_min_r_real)), float(np.degrees(m2_maj_th))
        )

    return rsac1, rsac2
