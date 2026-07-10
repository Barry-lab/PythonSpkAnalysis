"""Watershed field segmentation and the field-shuffle grid-cell test.

Ports the whole of ``FieldShuffle/``:

- :func:`watershed_fields` — segment a rate map into firing fields
  (``watershed(-rm)`` + ``fs_fixWatershed.m``'s ridge-reassignment).
- :func:`field_bin_locations` — for each field, its peak and every other
  bin's polar offset from that peak (``fs_binMembLoc.m``).
- :func:`shuffle_fields` — randomly reallocate field centres and their
  member bins while preserving each field's internal shape as closely as
  possible (``fs_shufFields.m``).
- :func:`fast_gridness_with_brandon` — the four gridness measures
  ``fs_masterAnalysis.m`` computes per (real or shuffled) ratemap: standard
  and expanding, each with and without Brandon-method ellipse correction.
- :func:`field_shuffle_significance` — the Recording-level orchestration
  (``fs_masterAnalysis.m``): builds a ratemap, gets its four gridness
  measures, then compares them against ``n_shuffles`` iterations each of a
  temporal shuffle and this module's field shuffle.
- :func:`is_grid_cell` — compare a cell's gridness against its own shuffle
  distributions at a given percentile (``identify_grid_cell.m``).

Requires the ``scikit-image`` extra (``pip install "barrylab[fieldshuffle]"``).
"""

import numpy as np
from scipy.ndimage import label
from skimage.segmentation import watershed

from .grid import cross_correlogram, fast_gridness
from .grid_ellipse import brandon_regularise


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


def _as_random_state(rng):
    if isinstance(rng, np.random.RandomState):
        return rng
    return np.random.RandomState(rng)


def shuffle_fields(rate_map_shape, field_bin_loc, static_bin, rng=None):
    """Randomly reallocate field centres and their member bins.

    Port of ``fs_shufFields.m``. Randomly picks new locations for every
    field's centre, then works through each field's other bins in order of
    closeness to its own peak (with a random per-field rotational offset,
    since many bins tie on distance and would otherwise always be placed in
    the same angular order), placing each one at the nearest still-available
    bin to the field's *new* centre (also with a random per-field rotational
    offset applied to that ranking). Fields are processed round-robin (one
    bin from each field in turn, in a random field order) rather than one
    field at a time, so no single field systematically claims the best
    positions.

    .. note::
       Unlike every bit-exact-validated port in this package, this function's
       randomness cannot be matched against MATLAB's own (``rand``/
       ``randperm``, unseeded and running a different RNG algorithm) under
       any circumstances -- so faithfulness here means the *algorithm*, not
       specific numeric output. Given that, this operates throughout in
       0-based, **row-major** (C-order) flat indices -- matching
       :func:`field_bin_locations`'s own convention -- rather than MATLAB's
       1-based column-major linear indices, which would only add translation
       risk for zero validation benefit. See ``FINDINGS.md``.

    Parameters
    ----------
    rate_map_shape : (int, int)
        Shape of the ratemap ``field_bin_loc`` was derived from (only the
        size is used, matching ``smthRm`` in the original -- it never reads
        the ratemap's actual values).
    field_bin_loc : list of dict
        Output of :func:`field_bin_locations`.
    static_bin : array-like of int
        0-based row-major flat indices of bins that must not move (typically
        unvisited bins, ``np.flatnonzero(occupancy_map == 0)``).
    rng : int, numpy.random.RandomState, or None
        Seed or random state for reproducibility.

    Returns
    -------
    ndarray of int, shape ``rate_map_shape``
        A gather-index array: ``shuf_inds[i, j]`` is the flat index of the
        *original* bin whose value should now appear at ``(i, j)``. Apply to
        any same-shape array with
        ``values.ravel()[shuf_inds.ravel()].reshape(rate_map_shape)``.
    """
    rng = _as_random_state(rng)
    ny, nx = rate_map_shape
    n_total = ny * nx
    n_field = len(field_bin_loc)

    static_bin = np.asarray(static_bin, dtype=np.int64).ravel()
    all_bins = np.arange(n_total)
    bin_to_move = np.setdiff1d(all_bins, static_bin)
    n_bin = bin_to_move.size

    # Random field-processing order -- persists through every subsequent
    # field-indexed step below, matching MATLAB reordering fldBinLoc once
    # at the top and reusing that same variable throughout.
    fields = [field_bin_loc[i] for i in rng.permutation(n_field)]

    shuf_inds = np.zeros(n_total, dtype=np.int64)

    # --- reallocate field centres to random new bin locations ---
    fld_cent_ind = bin_to_move[rng.permutation(n_bin)[:n_field]]
    for k, fld in enumerate(fields):
        peak_r, peak_c = fld["peak_index"]
        shuf_inds[fld_cent_ind[k]] = peak_r * nx + peak_c
    shuf_inds[static_bin] = static_bin

    # --- each field's own bins, ordered by closeness to its OWN peak ---
    own_rnd_orient = rng.uniform(0.0, 2 * np.pi, size=n_field)
    processed = []
    for k, fld in enumerate(fields):
        other = fld["other_indices"]
        theta = np.asarray(fld["theta"], dtype=float)
        radius = np.asarray(fld["radius"], dtype=float)
        other_flat = np.array([r * nx + c for r, c in other], dtype=np.int64)

        if other_flat.size:
            keep = ~np.isin(other_flat, static_bin)
            other_flat, theta, radius = other_flat[keep], theta[keep], radius[keep]

        theta = np.mod(theta + own_rnd_orient[k], 2 * np.pi)
        # sortrows([theta, radius], [2, 1]): radius primary, theta secondary,
        # both ascending -- np.lexsort's LAST key is primary.
        if other_flat.size:
            other_flat = other_flat[np.lexsort((theta, radius))]
        processed.append(other_flat)

    # --- rank every bin in the ratemap by closeness to each NEW field centre ---
    all_rows, all_cols = np.divmod(all_bins, nx)
    new_rnd_orient = rng.uniform(0.0, 2 * np.pi, size=n_field)
    bin_ord_by_field = np.empty((n_field, n_total), dtype=np.int64)
    for k in range(n_field):
        cr, cc = divmod(int(fld_cent_ind[k]), nx)
        off_i = all_rows - cr
        off_j = all_cols - cc
        th = np.mod(np.arctan2(-off_i, off_j) + new_rnd_orient[k], 2 * np.pi)
        r = np.hypot(off_j, off_i)
        ranked = all_bins[np.lexsort((th, r))]
        taken = np.isin(ranked, fld_cent_ind) | np.isin(ranked, static_bin)
        bin_ord_by_field[k] = np.where(taken, -1, ranked)

    # --- round-robin bin-processing sequence: for each "slot" (rank), one
    # bin from each field in turn, skipping fields with fewer bins ---
    max_n_other = max((p.size for p in processed), default=0)
    bin_seq, fld_id_seq = [], []
    for slot in range(max_n_other):
        for k, p in enumerate(processed):
            if slot < p.size:
                bin_seq.append(p[slot])
                fld_id_seq.append(k)

    # --- reallocate: each bin goes to the closest still-available bin to
    # its (new, already-shuffled) field centre ---
    for bin_idx, fld_id in zip(bin_seq, fld_id_seq):
        row = bin_ord_by_field[fld_id]
        available = row[row != -1]
        if available.size == 0:
            # MATLAB's find(vector,1) on an all-taken row returns [], and
            # shufInds([])=x is a silent no-op there -- meaning MATLAB would
            # quietly return an incomplete/malformed permutation (some bin
            # never placed, its slot left at the preallocated 0) rather than
            # error. That can only happen if static_bin collectively excludes
            # so much of the map that the movable bins can't all be
            # reallocated -- raise clearly instead of silently corrupting the
            # output or crashing with a bare IndexError.
            raise ValueError(
                "shuffle_fields ran out of available destination bins -- "
                "static_bin excludes too much of the map for all movable "
                "field bins to be reallocated"
            )
        next_bin = available[0]
        shuf_inds[next_bin] = bin_idx
        bin_ord_by_field[bin_ord_by_field == next_bin] = -1

    return shuf_inds.reshape(rate_map_shape)
