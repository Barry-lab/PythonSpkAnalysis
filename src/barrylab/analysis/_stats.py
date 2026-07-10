"""Small shared statistical utilities used across analyses."""

import numpy as np
from scipy import stats as _stats

_CORR_FUNCS = {
    "pearson": _stats.pearsonr,
    "spearman": _stats.spearmanr,
    "kendall": _stats.kendalltau,
}


def nan_corr(var1, var2, corr_type="pearson"):
    """Correlation that ignores NaN/Inf pairs. Port of ``nan_corr.m``.

    Any pair where either value is non-finite is dropped before correlating
    (mirrors the MATLAB, which additionally screens complex values — not
    applicable to real-valued Python arrays).

    Parameters
    ----------
    var1, var2 : ndarray
        Same-length arrays (flattened if not already 1-D).
    corr_type : {"pearson", "spearman", "kendall"}

    Returns
    -------
    (r, p)
        ``(nan, nan)`` if no valid pairs remain.
    """
    v1 = np.ravel(np.asarray(var1, dtype=float))
    v2 = np.ravel(np.asarray(var2, dtype=float))
    if v1.shape != v2.shape:
        raise ValueError("var1 and var2 must have the same number of elements")

    corr_fn = _CORR_FUNCS.get(corr_type.lower())
    if corr_fn is None:
        raise ValueError(f"corr_type must be one of {sorted(_CORR_FUNCS)}")

    valid = np.isfinite(v1) & np.isfinite(v2)
    if valid.sum() == 0:
        return float("nan"), float("nan")
    r, p = corr_fn(v1[valid], v2[valid])
    return float(r), float(p)
