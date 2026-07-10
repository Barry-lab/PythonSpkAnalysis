"""Tests for cross-SAC phase properties (``xSacProps.m``)."""

import numpy as np

from barrylab.analysis.grid import cross_correlogram, spatial_autocorrelogram
from barrylab.analysis.grid_phase import cross_sac_props


def _hex_grid_map(n=61, spacing=12.0, sigma=2.5, offset=(0.0, 0.0)):
    yy, xx = np.mgrid[0:n, 0:n]
    a1 = np.array([spacing, 0.0])
    a2 = np.array([spacing * 0.5, spacing * np.sqrt(3) / 2])
    centre = np.array([n / 2, n / 2]) + np.array(offset)
    m = np.zeros((n, n))
    for i in range(-5, 6):
        for j in range(-5, 6):
            c = i * a1 + j * a2 + centre
            m += np.exp(-((xx - c[0]) ** 2 + (yy - c[1]) ** 2) / (2 * sigma**2))
    return m


def _clean_xsac_with_known_peaks(
    n=41, central_offset=(3, -2), radius=8, neighbour_angles_deg=(10, 70, 130, 190, 250, 310)
):
    """A hand-built xSAC with exactly-known peak positions and no flat/tied
    regions, isolating the geometry math from the zero-clipping fragility
    documented in cross_sac_props' docstring (see also FINDINGS.md)."""
    c = n // 2
    xsac = 0.001 * (np.arange(n)[:, None] + 1.3 * np.arange(n)[None, :])
    central_row, central_col = c + central_offset[0], c + central_offset[1]
    xsac[central_row, central_col] = 10.0
    neighbours = []
    for a in neighbour_angles_deg:
        rad = np.radians(a)
        dx, dy = radius * np.cos(rad), radius * np.sin(rad)
        row, col = int(round(central_row - dy)), int(round(central_col + dx))
        xsac[row, col] = 5.0
        neighbours.append((a, row, col))
    return xsac, c, (central_row, central_col), neighbours


def test_cross_correlogram_of_a_map_with_itself_matches_autocorrelogram():
    rm = _hex_grid_map()
    assert np.allclose(cross_correlogram(rm, rm), spatial_autocorrelogram(rm), equal_nan=True)


def test_cross_sac_props_recovers_known_peak_geometry_exactly():
    n = 41
    xsac, centre, (central_row, central_col), neighbours = _clean_xsac_with_known_peaks(n=n)

    props = cross_sac_props(xsac)

    expected_mag = np.hypot(central_row - centre, central_col - centre)
    expected_ang = np.arctan2(-(central_row - centre), central_col - centre) % (2 * np.pi)
    assert np.isclose(props["phase_mag"], expected_mag)
    assert np.isclose(props["phase_ang"], expected_ang)
    assert props["peak_coord"] == (central_row + 1, central_col + 1)
    assert np.isclose(props["peak_corr"], 10.0)

    # Exact expected (angle, distance) from each neighbour's actual (rounded
    # to integer pixel) placement -- not the pre-rounding ideal angle, since
    # placing a spike requires rounding to a pixel and that shifts its true
    # angle/distance from the central peak slightly.
    exact = []
    for _, row, col in neighbours:
        dist = np.hypot(row - central_row, col - central_col)
        ang = np.arctan2(-(row - central_row), col - central_col) % (2 * np.pi)
        exact.append((ang, dist))
    exact.sort(key=lambda t: t[0])
    expected_three = np.array(exact[:3])

    got = props["axis_ang_dist"][np.argsort(props["axis_ang_dist"][:, 0])]
    assert np.allclose(got, expected_three, atol=1e-9)


def test_cross_sac_props_requires_at_least_seven_peaks():
    n = 21
    xsac = 0.001 * (np.arange(n)[:, None] + 1.3 * np.arange(n)[None, :])
    xsac[10, 10] = 10.0  # a single unambiguous peak only
    import pytest

    with pytest.raises(ValueError):
        cross_sac_props(xsac)


def test_cross_sac_props_self_correlation_has_zero_phase_offset():
    # A self cross-correlation is the autocorrelogram: the true peak sits
    # exactly at the geometric centre (r=1), so it beats every zero-clipped
    # background pixel on distance regardless of the swamping issue below.
    rm = _hex_grid_map()
    xsac = cross_correlogram(rm, rm)
    props = cross_sac_props(xsac)
    assert np.isclose(props["phase_mag"], 0.0)
    assert np.isclose(props["peak_corr"], 1.0)


def test_cross_sac_props_swamped_by_zero_clipped_background_on_realistic_input():
    """Documents a genuine MATLAB fragility (xSacProps.m), faithfully
    replicated rather than fixed -- see FINDINGS.md and the cross_sac_props
    docstring. Zero-clipping every negative correlation bin before peak
    detection floods the candidate list with spuriously-tied zero-valued
    "peaks" on any realistic (non-self) correlogram, so the phase-offset
    result generally does NOT recover the true injected spatial offset
    between two grid modules -- peak_corr coming back as 0.0 is the tell.
    """
    rm = _hex_grid_map()
    rm2 = _hex_grid_map(offset=(1.5, 3.0))
    xsac = cross_correlogram(rm, rm2)
    props = cross_sac_props(xsac)
    assert props["peak_corr"] == 0.0
    assert not np.isclose(props["phase_mag"], np.hypot(3.0, 1.5), atol=0.5)
