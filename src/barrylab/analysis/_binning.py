"""General N-dimensional histogram binning. Port of ``histnd.m``.

This is a fresh, general implementation (not a refactor of the private 2-D/1-D
binning helpers already embedded in ``maps.py``/``direction.py`` — consolidating
those is a goal-2, after-parity item; see ``PORTING_PLAN.md``). It exists because
:mod:`pxd` is the first analysis that genuinely needs N-D (specifically 3-D:
direction x position-Y x position-X) binning that the existing 2-D-only helpers
can't express.
"""

import numpy as np

_EPS = np.finfo(float).eps


def histnd(data, ranges, bin_size):
    """N-dimensional histogram. Port of ``histnd.m``.

    Parameters
    ----------
    data : ndarray (n_points, n_dim)
        Points to bin.
    ranges : ndarray (n_dim,) or (2, n_dim)
        Row 0 (or the only row) is the maximum permitted value per dimension.
        An optional row 1 gives the minimum, subtracted from ``data`` before
        binning (matches ``histnd.m``'s optional second range row).
    bin_size : float or ndarray (n_dim,)
        Bin width per dimension.

    Returns
    -------
    (counts, bins)
        ``counts`` has shape ``ceil(ranges[0] / bin_size)`` per dimension.
        ``bins`` is ``(n_valid_points, n_dim)``, the 0-based bin index assigned
        to each point (rows containing NaN are dropped, matching ``histnd.m``'s
        NaN handling for 1-D/2-D inputs; for 3-D+ MATLAB's own NaN filter is
        actually skipped — a latent limitation not exercised by any current
        caller, since position/direction data reaching this function is always
        fully interpolated already. We drop NaN rows unconditionally here,
        which gives identical results for the well-formed inputs this is
        actually called with.)
    """
    data = np.asarray(data, dtype=float)
    ranges = np.atleast_2d(np.asarray(ranges, dtype=float))
    n_dim = data.shape[1]
    bin_size = np.broadcast_to(np.asarray(bin_size, dtype=float), (n_dim,))

    if ranges.shape[0] == 2:
        data = data - ranges[1]
    max_range = ranges[0]

    quotient = data / bin_size
    quotient = np.where(quotient == 0, quotient + _EPS, quotient)
    bins = np.ceil(quotient)

    max_bin = np.ceil(max_range / bin_size).astype(np.int64)

    valid = ~np.isnan(bins).any(axis=1)
    bins = bins[valid].astype(np.int64)
    bins = np.clip(bins, 1, None)
    bins = np.minimum(bins, max_bin)

    counts = np.zeros(tuple(max_bin))
    idx = tuple(bins[:, d] - 1 for d in range(n_dim))
    np.add.at(counts, idx, 1.0)
    return counts, bins - 1
