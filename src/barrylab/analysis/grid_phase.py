"""Cross spatial-autocorrelogram (xSAC) phase properties.

Port of ``xSacProps.m`` (``Analysis_Grid``): key metrics from the
cross-correlogram of two grid-cell rate maps (normally from the same grid
module), used to characterise their relative spatial phase offset. Feed it
the output of :func:`~barrylab.analysis.grid.cross_correlogram`.
"""

import numpy as np

from .grid import _raw_regional_max_peaks


def _cart2pol(x, y):
    """MATLAB ``cart2pol``: ``(theta, rho) = (atan2(y, x), hypot(x, y))``."""
    return np.arctan2(y, x), np.hypot(x, y)


def _wrap_to_2pi(theta):
    """MATLAB ``wrapTo2Pi``: wrap radians into ``[0, 2*pi)``.

    (MATLAB's version also special-cases exact positive multiples of
    ``2*pi`` to map to ``2*pi`` rather than ``0`` -- unreachable here since
    every angle passed through this module originates from ``atan2``, whose
    range ``(-pi, pi]`` never produces an exact multiple of ``2*pi`` other
    than ``0`` itself, so plain ``mod`` is exact for this module's inputs.)
    """
    return np.mod(theta, 2 * np.pi)


def cross_sac_props(xsac):
    """Phase offset and axis geometry from a cross spatial autocorrelogram.

    Port of ``xSacProps.m``. Finds the peak closest to the xSAC's geometric
    centre (the simplest measure of relative grid phase between the two
    source cells), then the six peaks surrounding *that* peak, to estimate
    candidate grid-axis directions.

    Requires at least 7 detected peaks in ``xsac`` (the central peak plus 6
    neighbours) -- exactly as intolerant of fewer as the original MATLAB,
    which indexes ``dist(2:7)`` unconditionally.

    .. warning::
       Confirmed via a hand-constructed, unambiguous synthetic xSAC that the
       core geometry here (angle/distance/peak-selection) is an exact
       translation of ``xSacProps.m``. But the MATLAB algorithm itself has a
       fragility worth knowing about, faithfully preserved rather than
       fixed (see ``FINDINGS.md``): it zero-clips every negative correlation
       bin *before* peak detection (``xsac(xsac<0)=0``), and
       ``imregionalmax`` treats a large flat region as trivially "locally
       maximal" everywhere within it. On a realistic correlogram -- where
       roughly half the bins are typically negative, away from genuine
       structure -- this floods the candidate-peak list with thousands of
       spuriously-tied, zero-valued "peaks". Because every returned field
       here is chosen by proximity (closest peak to the centre, then
       closest six peaks to *that* peak), a real correlogram's true
       hexagonal structure can be swamped by nearer, meaningless
       zero-plateau points -- ``peak_corr`` coming back as exactly ``0.0``
       is the tell. Check ``peak_corr`` before trusting ``phase_mag``/
       ``phase_ang``/``peak_coord``, and treat ``axis_ang_dist`` (which
       searches for neighbours even closer in than the already-fragile
       central peak) with real caution on non-idealised input.

    Parameters
    ----------
    xsac : ndarray
        A cross spatial autocorrelogram (e.g. from
        :func:`~barrylab.analysis.grid.cross_correlogram`), typically between
        two cells from the same grid module.

    Returns
    -------
    dict with keys:
        ``phase_mag`` -- distance (bins) of the xSAC's central peak (closest
            to the geometric centre) from that centre -- the simplest
            indication of relative grid phase offset.
        ``phase_ang`` -- angle (radians, CCW from +x, image-plot convention:
            +y is downward in array/row terms) of that peak from the centre.
        ``peak_coord`` -- ``(row, col)``, 1-based (matching this function's
            own MATLAB docstring, "m,n coordinates" -- unlike
            :func:`~barrylab.analysis.grid.grid_properties`'s ``peak_coord``,
            which is ``(x, y)``).
        ``peak_corr`` -- xsac value at that peak.
        ``axis_ang_dist`` -- ``(3, 2)`` array of ``(angle, distance)`` for the
            three smallest-angle of the six peaks surrounding the central
            peak (candidate grid axis directions).
    """
    xsac = np.real(np.asarray(xsac, dtype=float))
    xsac = np.where(np.isnan(xsac), -1.0, xsac)
    xsac = np.where(xsac < 0, 0.0, xsac)

    peak_m, peak_n = _raw_regional_max_peaks(xsac)  # 0-based
    if peak_m.size < 7:
        # MATLAB indexes dist(2:7) unconditionally and would error out here
        # too (an out-of-bounds index) -- raise explicitly rather than let
        # Python's silently-truncating slice semantics return a malformed
        # (<3, 2) axis_ang_dist.
        raise ValueError(
            "cross_sac_props needs at least 7 detected peaks (1 central + "
            "6 neighbours), found %d" % peak_m.size
        )

    n_rows, n_cols = xsac.shape
    # MATLAB centM=sizeM/2+0.5 is 1-based; the 0-based equivalent is -1 of that.
    cent_m = n_rows / 2.0 - 0.5
    cent_n = n_cols / 2.0 - 0.5

    ang, dist = _cart2pol(peak_n - cent_n, -(peak_m - cent_m))
    ang = _wrap_to_2pi(ang)
    order = np.argsort(dist, kind="stable")

    phase_mag = float(dist[order[0]])
    phase_ang = float(ang[order[0]])
    peak_row, peak_col = int(peak_m[order[0]]), int(peak_n[order[0]])
    peak_corr = float(xsac[peak_row, peak_col])

    # Six peaks around the central peak (not the xSAC's geometric centre).
    ang2, dist2 = _cart2pol(peak_n - peak_col, -(peak_m - peak_row))
    order2 = np.argsort(dist2, kind="stable")
    six = order2[1:7]  # drop index 0: the central peak itself (distance 0)
    dist6 = dist2[six]
    ang6 = _wrap_to_2pi(ang2[six])

    order3 = np.argsort(ang6, kind="stable")
    axis_ang_dist = np.column_stack([ang6[order3[:3]], dist6[order3[:3]]])

    return {
        "phase_mag": phase_mag,
        "phase_ang": phase_ang,
        "peak_coord": (peak_row + 1, peak_col + 1),
        "peak_corr": peak_corr,
        "axis_ang_dist": axis_ang_dist,
    }
