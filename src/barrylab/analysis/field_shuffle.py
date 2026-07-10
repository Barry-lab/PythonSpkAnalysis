"""Watershed field segmentation for the field-shuffle grid-cell test.

Ports the watershed-dependent field-detection stage of ``FieldShuffle/``:

- :func:`watershed_fields` — segment a rate map into firing fields
  (``watershed(-rm)`` + ``fs_fixWatershed.m``'s ridge-reassignment).
- :func:`field_bin_locations` — for each field, its peak and every other
  bin's polar offset from that peak (``fs_binMembLoc.m``).

Requires the ``scikit-image`` extra (``pip install "barrylab[fieldshuffle]"``).

.. note::
   ``fs_shufFields.m`` (the actual randomised field-reallocation shuffle) and
   the ``fs_masterAnalysis.m``/``identify_grid_cell.m`` orchestration (which
   also need ``fs_fastGridness.m``'s Langston-2011 expanding-gridness variant
   and ``brandonRegEcc.m``'s ellipse-corrected gridness) are **not yet
   ported** — they are a substantially more intricate algorithm than the
   segmentation stage here, and were not given the same careful
   read-source/port/validate cycle as the rest of this package due to time
   constraints. See ``PORTING_PLAN.md``.
"""

import numpy as np
from scipy.ndimage import label
from skimage.segmentation import watershed

_FOUR_CONN = np.array([[0, 1, 0], [1, 1, 1], [0, 1, 0]])


def _fix_watershed(rate_map, labels):
    """Reassign unallocated (0) watershed ridge bins to the 4-connected
    neighbour with the closest rate-map value. Port of ``fs_fixWatershed.m``.

    All reassignments are computed against the original ``labels`` (not bins
    already reassigned in this same call), matching the source's simultaneous
    (vectorised) computation.
    """
    original = labels
    result = labels.copy()
    unassigned = np.argwhere(original == 0)
    if len(unassigned) == 0:
        return result

    ny, nx = rate_map.shape
    offsets = ((-1, 0), (1, 0), (0, -1), (0, 1))
    for r, c in unassigned:
        best_diff, best_label = np.inf, 0
        for dr, dc in offsets:
            rr, cc = r + dr, c + dc
            if 0 <= rr < ny and 0 <= cc < nx and original[rr, cc] > 0:
                diff = abs(rate_map[rr, cc] - rate_map[r, c])
                if diff < best_diff:
                    best_diff, best_label = diff, original[rr, cc]
        result[r, c] = best_label
    return result


def watershed_fields(rate_map, occupancy=None):
    """Segment a smoothed rate map into firing fields via watershed.

    Port of the field-segmentation stage of ``fs_masterAnalysis.m``
    (``watershed(-wsRm)`` + ``fs_fixWatershed.m``). Unvisited bins (NaN in
    ``rate_map``, or zero in ``occupancy`` if given) are excluded before
    watershedding, matching ``wsRm(binXY==0) = 0``.

    Parameters
    ----------
    rate_map : ndarray
        A smoothed rate map (Gaussian smoothing recommended, matching the
        original's field-segmentation-specific smoothing choice).
    occupancy : ndarray, optional
        Dwell-time map, same shape; bins with zero occupancy are treated as
        unvisited. If omitted, only NaN in ``rate_map`` marks unvisited bins.

    Returns
    -------
    ndarray of int labels, same shape as ``rate_map``. Unvisited bins get
    whatever label the zeroed-out landscape naturally assigns them (matching
    MATLAB); label 0 should not normally remain after the ridge-fix step,
    except in the rare case of a ridge pixel with no allocated neighbour.
    """
    rm = np.asarray(rate_map, dtype=float).copy()
    unvisited = np.isnan(rm)
    if occupancy is not None:
        unvisited = unvisited | (np.asarray(occupancy) == 0)
    rm[unvisited] = 0.0

    # Matches MATLAB exactly: zero out unvisited bins in the landscape BEFORE
    # watershedding (wsRm(binXY==0)=0) and let the transform itself decide
    # their fate, rather than forcing them to label 0 afterwards (which would
    # fight the watershed's own segmentation and corrupt the ridge-only
    # semantics fs_fixWatershed relies on).
    raw_labels = watershed(-rm, connectivity=2, watershed_line=True)
    return _fix_watershed(rm, raw_labels)


def field_bin_locations(labels, rate_map):
    """Each field's peak bin, and every other bin's polar offset from it.

    Port of ``fs_binMembLoc.m``.

    Parameters
    ----------
    labels : ndarray of int
        Field assignment from :func:`watershed_fields` (or any integer
        labelling, e.g. from ``scipy.ndimage.label``).
    rate_map : ndarray
        The rate map the labels were derived from.

    Returns
    -------
    list of dict, one per unique label, each with:
        ``label`` — the field id,
        ``peak_index`` — ``(row, col)`` of the field's peak bin,
        ``other_indices`` — list of ``(row, col)`` for the field's other bins,
        ``theta`` — ndarray, angle (rad, CCW from +x) from peak to each other bin,
        ``radius`` — ndarray, distance from peak to each other bin,
        ``n_other`` — ``len(other_indices)``.
    """
    ny, nx = rate_map.shape
    result = []
    for lbl in np.unique(labels):
        mask = labels == lbl
        masked_rm = np.where(mask, rate_map, np.nan)
        peak_flat = int(np.nanargmax(masked_rm))
        peak_r, peak_c = divmod(peak_flat, nx)

        other_flat = np.flatnonzero(mask.ravel())
        other_flat = other_flat[other_flat != peak_flat]
        other_r, other_c = np.divmod(other_flat, nx)

        off_r = other_r - peak_r
        off_c = other_c - peak_c
        theta = np.mod(np.arctan2(-off_r, off_c), 2 * np.pi)
        radius = np.hypot(off_c, off_r)

        result.append({
            "label": int(lbl),
            "peak_index": (peak_r, peak_c),
            "other_indices": list(zip(other_r.tolist(), other_c.tolist())),
            "theta": theta,
            "radius": radius,
            "n_other": len(other_flat),
        })
    return result
