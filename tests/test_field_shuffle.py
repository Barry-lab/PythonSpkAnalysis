"""Tests for watershed field segmentation (field-shuffle grid-cell test, v0.7 partial)."""

import numpy as np
import pytest

skimage = pytest.importorskip("skimage", reason="scikit-image not installed (barrylab[fieldshuffle])")

from barrylab.analysis.field_shuffle import (  # noqa: E402
    _fix_watershed,
    field_bin_locations,
    shuffle_fields,
    watershed_fields,
)


def test_fix_watershed_assigns_all_ridge_bins():
    rm = np.array([[1.0, 2.0, 1.0], [2.0, 3.0, 2.0], [1.0, 2.0, 1.0]])
    labels = np.array([[0, 1, 0], [0, 1, 0], [0, 2, 0]])
    fixed = _fix_watershed(rm, labels)
    assert np.all(fixed != 0)
    assert fixed.shape == labels.shape


def test_fix_watershed_uses_frozen_original_not_in_progress_reassignments():
    # A single unassigned bin (1,1) sits between two allocated regions with
    # very different values; it should go to whichever is CLOSER in value,
    # using the ORIGINAL (pre-reassignment) neighbour labels.
    rm = np.array([[10.0, 10.0, 0.0], [10.0, 10.5, 0.0], [10.0, 10.0, 0.0]])
    labels = np.array([[1, 1, 2], [1, 0, 2], [1, 1, 2]])
    fixed = _fix_watershed(rm, labels)
    # (1,1)=10.5 is closer to the left region's 10.0 values than the right's 0.0
    assert fixed[1, 1] == 1


def test_fix_watershed_no_ridge_bins_is_a_noop():
    labels = np.array([[1, 1], [2, 2]])
    rm = np.array([[1.0, 2.0], [3.0, 4.0]])
    fixed = _fix_watershed(rm, labels)
    assert np.array_equal(fixed, labels)


def test_watershed_fields_produces_valid_labelling_on_synthetic_fields():
    yy, xx = np.mgrid[0:30, 0:30]
    rm = (5.0 * np.exp(-((xx - 8) ** 2 + (yy - 8) ** 2) / 10)
          + 5.0 * np.exp(-((xx - 22) ** 2 + (yy - 22) ** 2) / 10))
    labels = watershed_fields(rm)
    assert labels.shape == rm.shape
    assert labels.min() >= 1  # no unresolved ridge bins on this clean synthetic case
    # two well-separated peaks should produce at least two distinct fields
    assert len(np.unique(labels)) >= 2


def test_watershed_fields_respects_occupancy():
    rm = np.ones((10, 10)) * 5.0
    occ = np.ones((10, 10))
    occ[0:3, 0:3] = 0  # a corner never visited
    labels = watershed_fields(rm, occupancy=occ)
    assert labels.shape == rm.shape  # runs without error on partially-unvisited input


def test_field_bin_locations_peaks_match_independent_check():
    yy, xx = np.mgrid[0:20, 0:20]
    rm = 5.0 * np.exp(-((xx - 5) ** 2 + (yy - 5) ** 2) / 8) + 5.0 * np.exp(-((xx - 15) ** 2 + (yy - 15) ** 2) / 8)
    labels = watershed_fields(rm)
    locs = field_bin_locations(labels, rm)

    for f in locs:
        mask = labels == f["label"]
        expected_peak = np.unravel_index(np.nanargmax(np.where(mask, rm, np.nan)), rm.shape)
        assert expected_peak == f["peak_index"]
        # field size should equal n_other + the peak itself
        assert int(mask.sum()) == f["n_other"] + 1


def test_field_bin_locations_polar_offsets_are_consistent():
    yy, xx = np.mgrid[0:15, 0:15]
    rm = 5.0 * np.exp(-((xx - 7) ** 2 + (yy - 7) ** 2) / 6)
    labels = watershed_fields(rm)
    locs = field_bin_locations(labels, rm)
    for f in locs:
        assert len(f["theta"]) == f["n_other"]
        assert len(f["radius"]) == f["n_other"]
        assert np.all(f["radius"] >= 0)
        assert np.all((f["theta"] >= 0) & (f["theta"] < 2 * np.pi))


def test_watershed_fields_on_example_grid_cell(example_trial):
    import barrylab as bl

    rec = bl.load_dacq(example_trial)
    rm = bl.rate_map(rec, 2, 3, smoothing="gaussian", smooth_size=3)
    occ = bl.occupancy_map(rec)
    labels = watershed_fields(rm, occupancy=occ)
    locs = field_bin_locations(labels, rm)

    assert labels.shape == rm.shape
    assert len(locs) >= 1
    total_bins = sum(f["n_other"] + 1 for f in locs)
    assert total_bins == labels.size  # every bin belongs to exactly one field


# --------------------------------------------------------------------------
# shuffle_fields (fs_shufFields.m) -- involves genuine per-call randomness
# (unseeded rand()/randperm() in MATLAB), so bit-exact comparison is
# impossible regardless of implementation care; validated via structural
# invariants the algorithm must satisfy by construction, plus a
# closeness-preservation property check (its whole stated purpose).
# --------------------------------------------------------------------------
def _two_field_setup(n=20):
    yy, xx = np.mgrid[0:n, 0:n]
    rm = (5.0 * np.exp(-((xx - 5) ** 2 + (yy - 5) ** 2) / 6)
          + 5.0 * np.exp(-((xx - 14) ** 2 + (yy - 14) ** 2) / 6))
    labels = watershed_fields(rm)
    locs = field_bin_locations(labels, rm)
    static_bin = np.array([], dtype=np.int64)
    return rm, locs, static_bin


def test_shuffle_fields_is_a_valid_permutation_of_all_bins():
    rm, locs, static_bin = _two_field_setup()
    shuf_inds = shuffle_fields(rm.shape, locs, static_bin, rng=0)
    assert shuf_inds.shape == rm.shape
    assert sorted(shuf_inds.ravel().tolist()) == list(range(rm.size))


def test_shuffle_fields_leaves_static_bins_in_place():
    rm, locs, _ = _two_field_setup()
    n = rm.shape[1]
    # A compact corner, like a genuinely-unvisited patch of arena -- unlike a
    # bin scattered maximally evenly across the whole map (e.g. every 7th
    # bin), which leaves zero slack when two fields' watershed catchment
    # basins already tile the entire map with no "background" area, and can
    # make the (greedy, faithfully-ported) reallocation genuinely
    # unsatisfiable -- see test_shuffle_fields_raises_when_over_constrained.
    static_bin = np.array([r * n + c for r in range(3) for c in range(3)], dtype=np.int64)
    shuf_inds = shuffle_fields(rm.shape, locs, static_bin, rng=1)
    assert np.array_equal(shuf_inds.ravel()[static_bin], static_bin)


def test_shuffle_fields_raises_when_over_constrained():
    # Documents a genuine MATLAB fragility (fs_shufFields.m), not fixed: if
    # static_bin excludes too much of the map relative to how tightly the
    # fields already tile it (e.g. two watershed catchment basins covering
    # every bin, no unlabelled "background" slack at all), the greedy
    # round-robin reallocation can run out of destinations for a field
    # partway through. MATLAB's find(vector,1) returns [] in that case and
    # shufInds([])=x is a silent no-op there, meaning MATLAB would quietly
    # return an incomplete/malformed permutation rather than fail loudly.
    # Raising clearly here is safer than reproducing that silent corruption
    # or crashing with a bare IndexError.
    rm, locs, _ = _two_field_setup()
    static_bin = np.arange(0, rm.size, 7)  # scattered maximally -- zero slack
    import pytest
    with pytest.raises(ValueError):
        shuffle_fields(rm.shape, locs, static_bin, rng=1)


def test_shuffle_fields_is_reproducible_with_a_seed():
    rm, locs, static_bin = _two_field_setup()
    a = shuffle_fields(rm.shape, locs, static_bin, rng=42)
    b = shuffle_fields(rm.shape, locs, static_bin, rng=42)
    assert np.array_equal(a, b)


def test_shuffle_fields_typically_moves_field_centres():
    rm, locs, static_bin = _two_field_setup()
    shuf_inds = shuffle_fields(rm.shape, locs, static_bin, rng=7)
    nx = rm.shape[1]
    for f in locs:
        peak_r, peak_c = f["peak_index"]
        peak_flat = peak_r * nx + peak_c
        # the bin that now HOLDS this field's original peak value is where
        # shuf_inds equals peak_flat -- find it and confirm it moved.
        new_home = np.flatnonzero(shuf_inds.ravel() == peak_flat)
        assert new_home.size == 1
        assert new_home[0] != peak_flat


def test_shuffle_fields_preserves_relative_field_compactness():
    # The algorithm's whole stated purpose: bins belonging to a field should
    # remain clustered around the field's NEW centre after reallocation, not
    # scattered randomly across the map. Check that the reallocated bins
    # belonging to each field are, on average, much closer to that field's
    # new centre than a uniformly-random placement would be.
    rm, locs, static_bin = _two_field_setup()
    shuf_inds = shuffle_fields(rm.shape, locs, static_bin, rng=3)
    ny, nx = rm.shape
    new_pos_of = {}
    for flat_old in range(rm.size):
        new_flat = int(np.flatnonzero(shuf_inds.ravel() == flat_old)[0])
        new_pos_of[flat_old] = divmod(new_flat, nx)

    max_possible_dist = np.hypot(ny, nx)
    for f in locs:
        peak_r, peak_c = f["peak_index"]
        peak_flat = peak_r * nx + peak_c
        new_centre = new_pos_of[peak_flat]
        dists = []
        for r, c in f["other_indices"]:
            old_flat = r * nx + c
            new_r, new_c = new_pos_of[old_flat]
            dists.append(np.hypot(new_r - new_centre[0], new_c - new_centre[1]))
        if dists:
            assert np.mean(dists) < 0.5 * max_possible_dist


def test_shuffle_fields_runs_on_example_grid_cell(example_trial):
    import barrylab as bl

    rec = bl.load_dacq(example_trial)
    rm = bl.rate_map(rec, 2, 3, smoothing="gaussian", smooth_size=3)
    occ = bl.occupancy_map(rec)
    labels = watershed_fields(rm, occupancy=occ)
    locs = field_bin_locations(labels, rm)
    static_bin = np.flatnonzero(occ.ravel() == 0)

    shuf_inds = shuffle_fields(rm.shape, locs, static_bin, rng=0)
    assert shuf_inds.shape == rm.shape
    assert sorted(shuf_inds.ravel().tolist()) == list(range(rm.size))
    assert np.array_equal(shuf_inds.ravel()[static_bin], static_bin)
