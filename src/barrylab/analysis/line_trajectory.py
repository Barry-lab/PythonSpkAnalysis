"""Best-fit straight-line trajectory through a decoded posterior matrix.

Ports the ``lineTraj_*`` group in ``Decoding/``, used to summarise replay
events: given a posterior probability matrix (space x time), find the
constant-speed straight line through it that captures the most probability
mass (e.g. for scoring candidate replay trajectories).

- :func:`best_line` -- exhaustive search over line gradient/offset for the
  best fit (``lineTraj_decode.m``).
- :func:`line_residual` -- score an arbitrary line-band mask against a pMat
  (``lineTraj_get_res.m``).
- :func:`reget_line_band` -- rebuild a line-band mask for a different-length
  pMat with the same spatial extent (``lineTraj_reget_lnBnd.m``).

The three functions are designed to round-trip: applying
:func:`reget_line_band` + :func:`line_residual` to the line found by
:func:`best_line` reproduces its own reported residual exactly (this is
literally what ``lineTraj_testScript.m`` checks, and is used here as the
validation oracle -- there is no MATLAB runtime on this machine to generate
independent reference numbers from the shipped ``pMat_data.mat``, see
``FINDINGS.md``).
"""

import numpy as np
from scipy.ndimage import convolve1d
from scipy.signal import convolve2d


def _matlab_round(x):
    """MATLAB ``round()``: half away from zero (numpy rounds half to even)."""
    return np.sign(x) * np.floor(np.abs(x) + 0.5)


def _define_line(sz_pm2, spd_to_test, y_range):
    """Build a line mask and its band. Port of ``lineTraj_define_line.m``.

    Parameters
    ----------
    sz_pm2 : int
        Width (time bins) of the target pMat.
    spd_to_test : float
        Gradient of the line, in spatial bins per time bin.
    y_range : int
        Half-width of the band around the line (bins above/below).

    Returns
    -------
    (tst_ln, tst_ln_bnd)
        Both ``(n_rows, sz_pm2)``, zero/one masks -- the bare line and the
        line plus its band.
    """
    if np.ceil(np.abs(spd_to_test)) == 0:
        tst_ln = np.ones((1, sz_pm2))
    else:
        n_rows = int(_matlab_round((sz_pm2 - 1) * abs(spd_to_test))) + 1
        tst_ln = np.zeros((n_rows, sz_pm2))
        tst_ln_y = _matlab_round(np.arange(sz_pm2) * abs(spd_to_test)).astype(np.int64)
        tst_ln[tst_ln_y, np.arange(sz_pm2)] = 1.0

    kern_len = 2 * y_range + 1
    tst_ln = np.pad(tst_ln, ((y_range, y_range), (0, 0)), mode="constant")
    # filter2(kern, tstLn) with a symmetric all-ones column kernel, 'same'
    # mode -- equivalent to a rolling sum along axis 0 (verified against a
    # brute-force conv2/rot90 translation before porting; see FINDINGS.md).
    tst_ln_bnd = convolve1d(tst_ln, np.ones(kern_len), axis=0, mode="constant", cval=0.0)

    if np.sign(spd_to_test) == -1:
        tst_ln = tst_ln[::-1, :]
        tst_ln_bnd = tst_ln_bnd[::-1, :]
    return tst_ln, tst_ln_bnd


def _il_filter2(a, b):
    """``conv2(b, rot90(a,2), 'valid')`` -- MATLAB's ``il_filter2`` inline fn."""
    return convolve2d(b, np.rot90(a, 2), mode="valid")


def best_line(pmat, t_bin_sz, s_bin_sz, y_range):
    """Find the best constant-speed straight-line trajectory through ``pmat``.

    Port of ``lineTraj_decode.m``. Exhaustively tests line gradients over a
    fixed speed range (100 to 5000 cm/s, positive and negative, in 50 cm/s
    steps -- hard-coded in the original MATLAB and kept as-is) and, for each
    gradient, all valid vertical offsets, scoring each candidate line by the
    mean probability mass of ``pmat`` falling within a band around it. Lines
    shorter than a third of ``pmat``'s width are discarded.

    Parameters
    ----------
    pmat : ndarray (n_pos_bin, n_time_bin)
        Posterior probability matrix (columns should sum to ~1; an all-zero
        column means "no data" and is fine). NaNs are treated as 0.
    t_bin_sz : float
        Temporal bin size in seconds.
    s_bin_sz : float
        Spatial bin size in cm.
    y_range : int
        Half-width (in spatial bins) of the band scored around each candidate
        line.

    Returns
    -------
    (best_res, best_spd, best_y_interc, best_line_band)
        ``best_res`` -- best line's mean band probability, normalised by the
        total probability mass in ``pmat`` (0 to 1); ``nan`` if no line of
        sufficient length was found.
        ``best_spd`` -- ``(2,)``: gradient in position-bins/time-bin, then in
        cm/s.
        ``best_y_interc`` -- y-intercept (position bin) of the best line.
        ``best_line_band`` -- ``(n_pos_bin, n_time_bin)`` mask of the bins
        contributing to the best residual.
    """
    pmat = np.where(np.isnan(pmat), 0.0, np.asarray(pmat, dtype=float))
    sz_pm1, sz_pm2 = pmat.shape

    spd_bnd = (100.0, 5000.0)
    spd_inc = 50.0
    min_lgth = sz_pm2 / 3.0

    n_pos = int(np.floor((spd_bnd[1] - (-spd_bnd[1])) / spd_inc)) + 1
    spd_to_test = -spd_bnd[1] + spd_inc * np.arange(n_pos)
    spd_to_test = spd_to_test[~((spd_to_test > -spd_bnd[0]) & (spd_to_test < spd_bnd[0]))]
    spd_to_test = (spd_to_test / s_bin_sz) * t_bin_sz
    n_spd_to_test = spd_to_test.size

    mean_p = np.zeros(n_spd_to_test)
    best_offset = np.zeros(n_spd_to_test, dtype=np.int64)
    pmat_ones = np.ones_like(pmat)

    for nn in range(n_spd_to_test):
        tst_ln, tst_ln_bnd = _define_line(sz_pm2, spd_to_test[nn], y_range)
        pad_rows = tst_ln.shape[0]

        pmat_padded = np.pad(pmat, ((pad_rows, pad_rows), (0, 0)), mode="constant")
        tmp_mean_p = _il_filter2(tst_ln_bnd, pmat_padded)

        ones_padded = np.pad(pmat_ones, ((pad_rows, pad_rows), (0, 0)), mode="constant")
        tmp_line_length = _il_filter2(tst_ln, ones_padded)

        # tst_ln_bnd/tst_ln always span pmat's full (unpadded) width, so the
        # 'valid' convolution collapses to a single column -- this is a plain
        # vector max, not a genuine 2-D tie-break.
        tmp_mean_p_col = tmp_mean_p[:, 0].copy()
        tmp_mean_p_col[tmp_line_length[:, 0] < min_lgth] = 0.0

        flat_idx = int(np.argmax(tmp_mean_p_col))
        mean_p[nn] = tmp_mean_p_col[flat_idx]
        best_offset[nn] = flat_idx

    tmp_ind = int(np.argmax(mean_p))
    best_res = mean_p[tmp_ind] / pmat.sum()

    best_spd = np.array([spd_to_test[tmp_ind], 0.0])
    best_spd[1] = best_spd[0] * s_bin_sz / t_bin_sz

    tst_ln, _ = _define_line(sz_pm2, spd_to_test[tmp_ind], y_range)

    # bstYInterc = bestOffSet(tmpInd) + find(tstLn(:,1)) - size(tstLn,1) - 2
    # (MATLAB, 1-based). Kept in explicit 1-based form here so it can be
    # audited line-by-line against the source rather than pre-simplified.
    r_1based = int(np.argmax(tst_ln[:, 0] != 0)) + 1  # find(tstLn(:,1)), first match
    best_offset_1based = int(best_offset[tmp_ind]) + 1
    best_y_interc = best_offset_1based + r_1based - tst_ln.shape[0] - 2

    best_line_band = reget_line_band(pmat.shape, best_spd[0], best_y_interc, y_range)

    if best_res == 0:
        best_res = np.nan
        best_y_interc = np.nan
        best_line_band = np.full_like(best_line_band, np.nan)
        best_spd = np.full(2, np.nan)

    return best_res, best_spd, best_y_interc, best_line_band


def reget_line_band(pmat_shape, spd, y_interc, y_range):
    """Rebuild a line-band mask for a (possibly different-length) pMat.

    Port of ``lineTraj_reget_lnBnd.m``. Lets a line found on one pMat (e.g.
    via :func:`best_line`) be re-applied to another pMat with the same
    number of spatial bins but a different number of time bins, to measure
    "chance" overlap.

    Parameters
    ----------
    pmat_shape : (int, int)
        Target ``(n_pos_bin, n_time_bin)``.
    spd : float
        Line gradient in position-bins per time-bin (``best_spd[0]`` from
        :func:`best_line`).
    y_interc : float
        Y-intercept (position bin) of the line.
    y_range : int
        Half-width (in spatial bins) of the band around the line.

    Returns
    -------
    ndarray ``pmat_shape``
        Zero/one band mask.
    """
    n_rows, n_cols = pmat_shape
    ln_bnd = np.zeros((n_rows * 3, n_cols))

    ln_bnd_x = np.arange(n_cols)
    ln_bnd_y = _matlab_round(ln_bnd_x * spd) + y_interc + n_rows
    in_range = (ln_bnd_y >= 0) & (ln_bnd_y < n_rows * 3)
    ln_bnd_x = ln_bnd_x[in_range]
    ln_bnd_y = ln_bnd_y[in_range].astype(np.int64)

    ln_bnd[ln_bnd_y, ln_bnd_x] = 1.0

    kern_len = 2 * y_range + 1
    ln_bnd = convolve1d(ln_bnd, np.ones(kern_len), axis=0, mode="constant", cval=0.0)

    return ln_bnd[n_rows: 2 * n_rows, :]


def line_residual(pmat, line_band):
    """Fraction of ``pmat``'s probability mass covered by a line-band mask.

    Port of ``lineTraj_get_res.m``.

    Parameters
    ----------
    pmat : ndarray (n_pos_bin, n_time_bin)
    line_band : ndarray (n_pos_bin, n_time_bin)
        A mask from :func:`best_line` or :func:`reget_line_band`.

    Returns
    -------
    float
        Between 0 and 1.
    """
    pmat = np.where(np.isnan(pmat), 0.0, np.asarray(pmat, dtype=float))
    res = pmat * line_band
    return res.sum() / pmat.sum()
