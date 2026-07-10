"""Circular-linear regression for theta phase precession.

Ports ``circ_regress_cb.m``: fits a line relating a linear variable (typically
in-field position) to a circular variable (theta phase, mod 2*pi) by
maximising the mean resultant length of the phase residuals, then assesses
significance by shuffling. Based on Kempter et al. (2012) and Fisher (1993)
"Statistical Analysis of Circular Data" section 6.3.3.
"""

import numpy as np
from scipy.optimize import minimize_scalar

from .direction import circular_mean


def _circ_r(alpha):
    """Mean resultant length of a set of angles (radians)."""
    return np.abs(np.mean(np.exp(1j * alpha)))


def _cost(slope, x, t):
    """Negative mean resultant length of the phase residual -- minimised
    (i.e. resultant length maximised) to find the best-fit slope."""
    return -_circ_r(t - slope * x)


def _kempter_corr(theta, t):
    """Circular-circular correlation (Kempter et al. 2012, eq. 3)."""
    t_bar = circular_mean(t)
    theta_bar = circular_mean(theta)
    num = np.sum(np.sin(t - t_bar) * np.sin(theta - theta_bar))
    denom = np.sqrt(np.sum(np.sin(t - t_bar) ** 2) * np.sum(np.sin(theta - theta_bar) ** 2))
    return num / denom


def _lcc(x, t, max_slope, bounds_frac):
    lo, hi = bounds_frac[0] * max_slope, bounds_frac[1] * max_slope
    result = minimize_scalar(
        _cost, args=(x, t), bounds=(lo, hi), method="bounded", options={"xatol": 1e-2}
    )
    slope_opt = result.x
    fval_opt = result.fun
    phase_opt = np.arctan2(np.sum(np.sin(t - slope_opt * x)), np.sum(np.cos(t - slope_opt * x)))

    theta = np.mod(np.abs(slope_opt) * x, 2 * np.pi)
    corr = _kempter_corr(theta, t)
    if np.sign(corr) != np.sign(slope_opt):
        corr = -corr

    return {"slope_opt": slope_opt, "fval_opt": fval_opt, "phase_opt": phase_opt, "corr": corr}


def circular_linear_regression(x, t, k=1000, pp_dir=0, random_state=None):
    """Circular-linear regression between a linear variable and a phase.

    Port of ``circ_regress_cb.m``. Finds the slope maximising the mean
    resultant length of ``t - slope*x`` (i.e. the best-fit line making the
    phase residuals as tightly clustered as possible), reports a
    circular-circular correlation for the fit, and estimates a p-value by
    shuffling ``t`` relative to ``x`` ``k`` times and refitting.

    Parameters
    ----------
    x : ndarray
        In-field positions (linear variable).
    t : ndarray
        Phases in radians, in ``[0, 2*pi)``.
    k : int
        Number of shuffles for the permutation p-value (MATLAB recommends
        1000 for a quick estimate, 10000 for a stable one).
    pp_dir : {0, 1, -1}
        Constrain the fitted slope's sign: ``0`` unconstrained, ``1`` forces a
        positive slope, ``-1`` forces a negative slope. The slope magnitude is
        always capped at ``4*pi / (max(x) - min(x))``.
    random_state : int or numpy.random.RandomState, optional
        Seed/generator for the shuffle (MATLAB's own ``randperm`` is likewise
        nondeterministic by default; passing a seed is purely for
        reproducibility, not a behavioural change from the original).

    Returns
    -------
    dict with ``slope_opt`` (rad/unit-x), ``fval_opt``, ``phase_opt`` (rad,
    intercept at x=0), ``corr`` (Kempter circular-circular R), ``psim``
    (permutation p-value). All ``nan`` if fewer than 3 valid (non-NaN) points
    are supplied.
    """
    x = np.asarray(x, dtype=float).ravel()
    t = np.asarray(t, dtype=float).ravel()
    if x.shape != t.shape:
        raise ValueError("x and t must be the same size")

    valid = ~(np.isnan(x) | np.isnan(t))
    x, t = x[valid], t[valid]

    if np.any(t > 2 * np.pi) or np.any(t < 0):
        raise ValueError("t must be in radians, in the range [0, 2*pi]")

    if len(x) < 3:
        return {
            "slope_opt": np.nan,
            "phase_opt": np.nan,
            "fval_opt": np.nan,
            "corr": np.nan,
            "psim": np.nan,
        }

    if pp_dir == -1:
        bounds_frac = (-1.0, np.finfo(float).eps)
    elif pp_dir == 1:
        bounds_frac = (np.finfo(float).eps, 1.0)
    else:
        bounds_frac = (-1.0, 1.0 - np.finfo(float).eps)

    max_slope = (4 * np.pi) / (x.max() - x.min())
    fit = _lcc(x, t, max_slope, bounds_frac)

    rng = (
        random_state
        if isinstance(random_state, np.random.RandomState)
        else np.random.RandomState(random_state)
    )
    n = len(t)
    corr_dist = np.empty(k)
    for i in range(k):
        shuffled = t[rng.permutation(n)]
        corr_dist[i] = _lcc(x, shuffled, max_slope, bounds_frac)["corr"]

    if np.isnan(fit["corr"]):
        fit["psim"] = np.nan
    elif fit["corr"] >= 0:
        fit["psim"] = float(np.mean(corr_dist >= fit["corr"]))
    else:
        fit["psim"] = float(np.mean(corr_dist <= fit["corr"]))

    return fit
