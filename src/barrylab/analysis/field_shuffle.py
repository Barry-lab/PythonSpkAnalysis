"""Watershed field segmentation and the field-shuffle grid-cell test.

Ports the whole of ``FieldShuffle/``:

- :func:`watershed_fields` — segment a rate map into firing fields
  (``watershed(-rm)`` + ``fs_fixWatershed.m``'s ridge-reassignment).
- :func:`field_bin_locations` — for each field, its peak and every other
  bin's polar offset from that peak (``fs_binMembLoc.m``).
- :func:`shuffle_fields` — randomly reallocate field centres and their
  member bins while preserving each field's internal shape as closely as
  possible (``fs_shufFields.m``).
- :func:`field_shuffle_significance` — the Recording-level orchestration
  (``fs_masterAnalysis.m``): builds a ratemap, gets its four gridness
  measures (standard/expanding, each with and without Brandon-method
  ellipse correction), then compares them against ``n_shuffles`` iterations
  each of a temporal shuffle and this module's field shuffle.
- :func:`is_grid_cell` — compare a cell's gridness against its own shuffle
  distributions at a given percentile (``identify_grid_cell.m``).

Requires the ``scikit-image`` extra (``pip install "barrylab[fieldshuffle]"``).
"""

import numpy as np
from skimage.segmentation import watershed

from .grid import fast_gridness, spatial_autocorrelogram
from .grid_ellipse import brandon_regularise
from .maps import occupancy_map, smooth_ratemap, spike_map


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

       A field whose *peak* bin is itself in ``static_bin`` is treated as
       entirely static (none of its bins move) rather than relocated like a
       normal field. Unmarked watershed always tiles the *entire* ratemap
       into catchment basins, including never-visited regions -- a field
       whose single highest-rate bin happens to be one of those unvisited
       bins is a watershed artefact of a never-visited area, not a genuine
       firing field, and both MATLAB and a literal translation of it are
       vulnerable to a real, not-contrived capacity mismatch if it's
       relocated anyway: that field's peak would still consume a slot from
       the movable-bin pool via random reallocation (nothing in
       ``fs_binMembLoc.m``/:func:`field_bin_locations` checks whether a
       peak coincides with an unvisited bin) while contributing nothing back
       to supply, since it was excluded from ``static_bin`` accounting on
       the demand side. Confirmed via real data (the shipped example
       trial's non-grid cell t2c4 has exactly this: a tiny field whose peak
       sits on an unvisited bin, deterministically over-constraining the
       allocation by exactly one slot on every seed tried) that this is a
       likely, not rare, occurrence for cells with sparse/patchy firing --
       exactly the kind of cell this whole function exists to correctly
       call non-grid. See ``FINDINGS.md``.

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

    static_bin = np.asarray(static_bin, dtype=np.int64).ravel()

    # A field whose own peak is a static bin is a watershed artefact of a
    # never-visited region, not a genuine field -- exclude it entirely
    # (every one of its bins joins the static set) rather than relocate it.
    # See the docstring note above for why this matters for feasibility.
    static_set = set(static_bin.tolist())
    genuine_fields, extra_static = [], []
    for fld in field_bin_loc:
        peak_r, peak_c = fld["peak_index"]
        if (peak_r * nx + peak_c) in static_set:
            extra_static.append(peak_r * nx + peak_c)
            extra_static.extend(r * nx + c for r, c in fld["other_indices"])
        else:
            genuine_fields.append(fld)
    if extra_static:
        static_bin = np.union1d(static_bin, np.array(extra_static, dtype=np.int64))
    field_bin_loc = genuine_fields
    n_field = len(field_bin_loc)

    all_bins = np.arange(n_total)
    bin_to_move = np.setdiff1d(all_bins, static_bin)
    n_bin = bin_to_move.size
    if n_bin < n_field:
        # Should not be reachable via any self-consistent (field_bin_loc,
        # static_bin) pair post-exclusion above: a genuine field's peak is
        # never static by definition (that's exactly what triggers
        # exclusion), so every genuine field's peak is always itself part of
        # the movable pool, guaranteeing n_bin >= n_field. Kept as a
        # defensive guard for a caller passing mismatched/inconsistent data
        # (e.g. a static_bin from a different occupancy map) -- MATLAB's
        # fldCentInd(1:nField) would index past the end of a shorter array
        # and error there too; Python's slicing silently truncates instead,
        # which would otherwise surface later as a confusing IndexError once
        # the truncated fld_cent_ind ran out mid-loop.
        raise ValueError(
            "shuffle_fields has more genuine fields (%d) than movable bins "
            "(%d) to give them new centres" % (n_field, n_bin)
        )

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


def _brandon_gridness(sac):
    """Best-of-two-methods Brandon-corrected gridness. Port of the
    ``sf_getBrandonG`` subfunction embedded in ``fs_masterAnalysis.m``."""
    rsac1, rsac2 = brandon_regularise(sac)
    std1, exp1, _ = fast_gridness(rsac1)
    std2, exp2, _ = fast_gridness(rsac2)
    # np.fmax (not max()/np.maximum) to match MATLAB max()'s NaN-ignoring
    # behaviour: if only one method found enough peaks, use that one.
    return float(np.fmax(std1, std2)), float(np.fmax(exp1, exp2))


def field_shuffle_significance(recording, electrode, cluster, bin_cm=2.0, n_shuffles=100,
                                smooth_kern_size=1.5, smooth_kern_type="gaussian",
                                min_temporal_rotation_s=30.0, keep_ratemaps=False, rng=None):
    """Gridness of a cell vs. its own temporal- and field-shuffle distributions.

    A ``Recording``-level redesign of ``fs_masterAnalysis.m`` (like
    :func:`~barrylab.analysis.decoding.crunch_decode` -- validate the
    component functions individually; this end-to-end orchestrator is not
    bit-exact, and its field-shuffle step involves genuine randomness that
    can never be bit-exact regardless, see :func:`shuffle_fields`).

    For the real cell, builds a ratemap and gets four gridness measures
    (standard and expanding, each with and without Brandon-method ellipse
    correction — see :func:`~barrylab.analysis.grid.fast_gridness` and
    :func:`~barrylab.analysis.grid_ellipse.brandon_regularise`). Then repeats
    this ``n_shuffles`` times for two null distributions: a **temporal**
    shuffle (spikes' position-sample indices rotated by a random amount,
    with at least ``min_temporal_rotation_s`` seconds either side of zero)
    and a **field** shuffle (:func:`shuffle_fields`, on fields segmented once
    from the real, Gaussian-smoothed ratemap).

    .. note::
       The original MATLAB hardcodes a 50 Hz position sample rate in this
       rotation-range calculation (``trlLgth=trlLgthPosSamp/(50*60)``,
       ``rotRng=[...]*50``) rather than reading it from the data -- fine for
       its DACQ-only context (always 50 Hz) but not generalisable. This
       ``Recording``-based redesign uses ``recording.pos_sample_rate``
       directly instead, consistent with how ``crunch_decode`` already
       generalises the same kind of MATLAB DACQ-specific assumption
       elsewhere in this package.

       Also **does not reproduce** a genuine MATLAB bug found while reading
       ``fs_masterAnalysis.m``: ``rotSpkPos(rotSpkPos==0)=nLpForShuf`` -- the
       modulo-wrapped position index 0 should remap to ``trlLgthPosSamp``
       (MATLAB is 1-based, so 0 is never a valid position index and must
       wrap to the *last* one), but the code substitutes the unrelated
       shuffle-iteration count instead, corrupting a small number of spikes
       per shuffle iteration (whichever land exactly on a multiple of the
       trial length after rotation). In 0-based Python this whole
       remapping step is moot -- ``0`` is already a perfectly valid index --
       so ``(spk_idx + rotation) % n_pos`` is used directly with no special
       case, sidestepping the bug entirely rather than needing to replicate
       or explicitly avoid it. See ``FINDINGS.md``.

    Parameters
    ----------
    recording : Recording
    electrode, cluster : int
    bin_cm : float
        Spatial bin size in cm for both the ratemap and field segmentation.
    n_shuffles : int
        Number of temporal and field shuffle iterations (MATLAB default 100,
        via ``identify_grid_cell.m``).
    smooth_kern_size, smooth_kern_type : float, str
        Smoothing for the *gridness* ratemaps (field segmentation always
        uses Gaussian smoothing with sigma 3, matching the original's fixed
        choice for that specific step).
    min_temporal_rotation_s : float
        Minimum temporal-shuffle rotation magnitude, in seconds.
    keep_ratemaps : bool
        If ``True``, also return every shuffle's ratemap (slower; mirrors
        MATLAB's ``nargout==5`` conditional via an explicit parameter instead
        of output-count introspection).
    rng : int, numpy.random.RandomState, or None
        Seed or random state.

    Returns
    -------
    dict with:
        ``gridness`` -- ``(4,)`` ``[standard, standard_brandon, expanding,
            expanding_brandon]`` for the real cell.
        ``grid_scale`` -- grid scale (bins) of the real cell's SAC.
        ``p95`` -- ``(4, 2)`` 95th-percentile shuffle threshold, columns
            ``[temporal, field]``.
        ``gridness_from_shuffle`` -- ``(4, 2, n_shuffles)`` every shuffle's
            gridness.
        ``ratemaps`` -- (only if ``keep_ratemaps``) dict with ``real``,
            ``temporal_shuffle`` (list), ``field_shuffle`` (list).
    """
    rng = _as_random_state(rng)
    elec = recording.electrode(electrode)
    spk_idx = elec.spike_pos_index(cluster)
    n_pos = recording.n_pos

    bin_occ = occupancy_map(recording, bin_cm)
    bin_spk = spike_map(recording, spk_idx, bin_cm)

    rm, _, _ = smooth_ratemap(bin_occ, bin_spk, smooth_kern_size, smooth_kern_type, "norm")
    rm_gauss, _, _ = smooth_ratemap(bin_occ, bin_spk, 3.0, "gaussian", "norm")

    sac = spatial_autocorrelogram(rm)
    std_g, exp_g, grid_scale = fast_gridness(sac)
    std_bg, exp_bg = _brandon_gridness(sac)
    grid = np.array([std_g, std_bg, exp_g, exp_bg])

    # --- field segmentation, once, on the real (Gaussian-smoothed) ratemap ---
    # Matches MATLAB's wsRm(binXY==0)=0 literally: the SAME zeroed array is
    # used for both watershedding and fs_binMembLoc's peak-finding, not
    # watershed_fields' own (NaN-based) internal zeroing.
    unvisited = bin_occ == 0
    ws_rm = np.where(unvisited, 0.0, rm_gauss)
    labels = watershed_fields(ws_rm)
    field_bin_loc = field_bin_locations(labels, ws_rm)
    static_bin = np.flatnonzero(unvisited.ravel())

    min_rot_samp = int(round(min_temporal_rotation_s * recording.pos_sample_rate))
    rot_low, rot_high = min_rot_samp, n_pos - min_rot_samp

    grid_from_shuffle = np.zeros((4, 2, n_shuffles))
    time_rms = [] if keep_ratemaps else None
    space_rms = [] if keep_ratemaps else None

    for h in range(n_shuffles):
        # 1) temporal shuffle
        rot_amount = rng.randint(rot_low, rot_high + 1)
        rot_spk_idx = (spk_idx + rot_amount) % n_pos
        bin_spk_rot = spike_map(recording, rot_spk_idx, bin_cm)
        rm_t, _, _ = smooth_ratemap(bin_occ, bin_spk_rot, smooth_kern_size, smooth_kern_type, "norm")
        if keep_ratemaps:
            time_rms.append(rm_t)
        sac_t = spatial_autocorrelogram(rm_t)
        t_std, t_exp, _ = fast_gridness(sac_t)
        t_std_b, t_exp_b = _brandon_gridness(sac_t)

        # 2) field shuffle
        shuf_inds = shuffle_fields(rm.shape, field_bin_loc, static_bin, rng=rng)
        flat = shuf_inds.ravel()
        bin_occ_shuf = bin_occ.ravel()[flat].reshape(rm.shape)
        bin_spk_shuf = bin_spk.ravel()[flat].reshape(rm.shape)
        rm_s, _, _ = smooth_ratemap(bin_occ_shuf, bin_spk_shuf, smooth_kern_size, smooth_kern_type, "norm")
        if keep_ratemaps:
            space_rms.append(rm_s)
        sac_s = spatial_autocorrelogram(rm_s)
        s_std, s_exp, _ = fast_gridness(sac_s)
        s_std_b, s_exp_b = _brandon_gridness(sac_s)

        grid_from_shuffle[:, 0, h] = [t_std, t_std_b, t_exp, t_exp_b]
        grid_from_shuffle[:, 1, h] = [s_std, s_std_b, s_exp, s_exp_b]

    p95 = np.nanpercentile(grid_from_shuffle, 95, axis=2)

    result = {
        "gridness": grid,
        "grid_scale": grid_scale,
        "p95": p95,
        "gridness_from_shuffle": grid_from_shuffle,
    }
    if keep_ratemaps:
        result["ratemaps"] = {"real": rm, "temporal_shuffle": time_rms, "field_shuffle": space_rms}
    return result


def is_grid_cell(recording, electrode, cluster, n_shuffles=100, bin_cm=2.0, rng=None):
    """Is this cell a grid cell, by comparison to its own shuffle distributions?

    Port of ``identify_grid_cell.m`` (redesigned around ``Recording``, like
    :func:`field_shuffle_significance`). The original hardcodes its shuffle
    count (100) and bin size (2 cm) rather than exposing them as arguments;
    kept as the defaults here but overridable.

    .. note::
       The original MATLAB calls ``fs_masterAnalysis(data.pos, spkPos, ...)``,
       referencing an undefined variable ``data`` (the function's own
       parameter is named ``pos``, not ``data.pos``) -- this would raise
       "undefined variable" in MATLAB on every single call, meaning
       ``identify_grid_cell.m`` was very likely never actually run as a
       standalone function (only ever exercised, if at all, by pasting its
       body into a script/workspace context where a variable happened to be
       named ``data``). Not something that can be "faithfully reproduced" --
       there is no well-defined numerical behaviour to preserve, only a
       crash -- so this redesign simply calls
       :func:`field_shuffle_significance` with the actual ``recording``
       argument. See ``FINDINGS.md``.

    Parameters
    ----------
    recording : Recording
    electrode, cluster : int
    n_shuffles : int
    bin_cm : float
    rng : int, numpy.random.RandomState, or None

    Returns
    -------
    dict with:
        ``is_grid`` -- ``(4, 2)`` bool, rows ``[standard, standard_brandon,
            expanding, expanding_brandon]``, columns ``[temporal, field]``
            shuffle -- whether the real cell's gridness exceeds that
            shuffle's 95th percentile.
        ``gridness`` -- ``(4,)`` the cell's own gridness values.
        ``threshold`` -- ``(4, 2)`` the 95th-percentile threshold used.
    """
    result = field_shuffle_significance(
        recording, electrode, cluster, bin_cm=bin_cm, n_shuffles=n_shuffles, rng=rng
    )
    is_grid = result["gridness"][:, None] > result["p95"]
    return {"is_grid": is_grid, "gridness": result["gridness"], "threshold": result["p95"]}
