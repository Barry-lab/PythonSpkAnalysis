"""Elliptical grid correction: fit and remove elliptical distortion from a SAC.

Ports ``gridEllipse_fit.m`` (fit an ellipse to a grid cell's six central SAC
peaks, giving x/y scale and eccentricity) and ``gridEllipse_correct.m``
(regularise an elliptical SAC back toward circular, so gridness can be
recomputed on the corrected map). ``regularise_eliptic_grid.m`` is a
byte-identical duplicate of ``gridEllipse_correct.m`` in the original — ported
once here.

The ellipse-fit itself (``_fit_ellipse``) is a least-squares conic fit
originally from the MathWorks File Exchange, embedded in ``gridEllipse_fit.m``
as a subfunction; ported faithfully for portability, as the original was.
"""

import numpy as np
from scipy.ndimage import map_coordinates, rotate

from .grid import _sac_peaks


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

    design = np.column_stack([x ** 2, x * y, y ** 2, x, y])
    coeffs, *_ = np.linalg.lstsq(design.T @ design, design.T @ np.ones_like(x), rcond=None)
    a, b, c, d, e = coeffs

    orientation_tol = 1e-3
    if min(abs(b / a) if a else np.inf, abs(b / c) if c else np.inf) > orientation_tol:
        orientation_rad = 0.5 * np.arctan2(b, c - a)
        cos_phi, sin_phi = np.cos(orientation_rad), np.sin(orientation_rad)
        a, b, c, d, e = (
            a * cos_phi ** 2 - b * cos_phi * sin_phi + c * sin_phi ** 2,
            0.0,
            a * sin_phi ** 2 + b * cos_phi * sin_phi + c * cos_phi ** 2,
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
    f = 1 + (d ** 2) / (4 * a) + (e ** 2) / (4 * c)
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
    nan_result = {"xy_scale": (np.nan, np.nan), "eccentricity": np.nan,
                  "orientation": np.nan, "ab_scale": (np.nan, np.nan)}
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

    eccentricity = np.sqrt((a ** 2 - b ** 2) / a ** 2)

    theta_axes = np.array([orient, np.pi / 2 + orient])
    ellipse_x = a * np.cos(theta_axes)
    ellipse_y = b * np.sin(theta_axes)
    xy_scale = np.sqrt(ellipse_x ** 2 + ellipse_y ** 2)

    return {"xy_scale": tuple(xy_scale), "eccentricity": float(eccentricity),
            "orientation": float(orient), "ab_scale": (float(a), float(b))}


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
    reg = reg[start_row: start_row + sac.shape[0], start_col: start_col + sac.shape[1]]
    return reg
