"""Tests for watershed field segmentation and the field-shuffle grid-cell test (v0.7)."""

import numpy as np
import pytest

skimage = pytest.importorskip(
    "skimage", reason="scikit-image not installed (barrylab[fieldshuffle])"
)

from barrylab.analysis.field_shuffle import (  # noqa: E402
    _brandon_gridness,
    _fix_watershed,
    field_bin_locations,
    field_shuffle_significance,
    is_grid_cell,
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
    rm = 5.0 * np.exp(-((xx - 8) ** 2 + (yy - 8) ** 2) / 10) + 5.0 * np.exp(
        -((xx - 22) ** 2 + (yy - 22) ** 2) / 10
    )
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
    rm = 5.0 * np.exp(-((xx - 5) ** 2 + (yy - 5) ** 2) / 8) + 5.0 * np.exp(
        -((xx - 15) ** 2 + (yy - 15) ** 2) / 8
    )
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
    rm = 5.0 * np.exp(-((xx - 5) ** 2 + (yy - 5) ** 2) / 6) + 5.0 * np.exp(
        -((xx - 14) ** 2 + (yy - 14) ** 2) / 6
    )
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
    static_bin = np.array([r * n + c for r in range(3) for c in range(3)], dtype=np.int64)
    shuf_inds = shuffle_fields(rm.shape, locs, static_bin, rng=1)
    assert np.array_equal(shuf_inds.ravel()[static_bin], static_bin)


def test_shuffle_fields_excludes_a_field_whose_peak_is_static():
    # Documents a real, likely-common (not contrived) MATLAB fragility found
    # via the shipped example trial's non-grid cell t2c4: unmarked watershed
    # tiles the ENTIRE ratemap into catchment basins, including never-visited
    # regions, so a tiny/marginal field's peak can itself land on an
    # unvisited (static) bin -- a plausible outcome for any cell with
    # sparse/patchy firing, exactly the kind of cell this whole module
    # exists to help classify as non-grid. If relocated like a normal field
    # anyway, that field's peak still consumes a slot from the movable pool
    # while contributing nothing to declared demand (having been excluded
    # from static_bin's accounting) -- deterministically over-constraining
    # the allocation by exactly one slot every time, confirmed by hand on the
    # real fixture: 0/30 seeds succeeded before this field-exclusion handling
    # was added, matching what would be MATLAB's own silent-corruption
    # failure mode there (see FINDINGS.md). Fixed by excluding any field
    # whose peak is itself static entirely (it's watershed noise from a
    # never-visited area, not a genuine field) rather than relocating it.
    rm = np.zeros((1, 6))
    locs = [
        {  # a genuine field: peak=0, others=1,2
            "label": 1,
            "peak_index": (0, 0),
            "n_other": 2,
            "other_indices": [(0, 1), (0, 2)],
            "theta": np.array([0.0, 0.0]),
            "radius": np.array([1.0, 2.0]),
        },
        {  # a degenerate "field" whose peak (bin 3) is itself unvisited
            "label": 2,
            "peak_index": (0, 3),
            "n_other": 2,
            "other_indices": [(0, 4), (0, 5)],
            "theta": np.array([0.0, 0.0]),
            "radius": np.array([1.0, 2.0]),
        },
    ]
    static_bin = np.array([3], dtype=np.int64)  # exactly field 2's peak

    shuf_inds = shuffle_fields(rm.shape, locs, static_bin, rng=0)
    assert sorted(shuf_inds.ravel().tolist()) == list(range(rm.size))
    # the whole degenerate field (peak + its declared others) stays in place,
    # not just the single bin explicitly marked static
    for flat in (3, 4, 5):
        assert shuf_inds.ravel()[flat] == flat


def test_shuffle_fields_runs_on_real_grid_and_non_grid_cells(example_trial):
    # The exact real-data scenario that motivated the peak-exclusion fix
    # above: t2c4 (a documented non-grid cell) has a tiny watershed field
    # whose peak sits on an unvisited bin, which reliably (0/30 seeds)
    # over-constrained the allocation before that fix. Confirms it now
    # succeeds across several seeds for every example cluster, grid or not.
    import barrylab as bl

    rec = bl.load_dacq(example_trial)
    occ = bl.occupancy_map(rec, 2.0)
    unvisited = occ == 0
    static_bin = np.flatnonzero(unvisited.ravel())
    for cluster in (1, 3, 4, 5):
        rm = bl.rate_map(rec, 2, cluster, smoothing="gaussian", smooth_size=3)
        # Zero (not NaN) unvisited bins before segmenting/locating, matching
        # field_shuffle_significance's own wsRm(binXY==0)=0 handling -- an
        # all-NaN field slice would otherwise crash field_bin_locations'
        # nanargmax.
        ws_rm = np.where(unvisited, 0.0, rm)
        labels = watershed_fields(ws_rm)
        locs = field_bin_locations(labels, ws_rm)
        for seed in range(5):
            shuf_inds = shuffle_fields(ws_rm.shape, locs, static_bin, rng=seed)
            assert sorted(shuf_inds.ravel().tolist()) == list(range(rm.size))


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


# --------------------------------------------------------------------------
# field_shuffle_significance / is_grid_cell (fs_masterAnalysis.m /
# identify_grid_cell.m) -- Recording-level redesigns (like crunch_decode),
# not bit-exact, and inheriting shuffle_fields' genuine randomness. See
# their docstrings/FINDINGS.md for two real MATLAB bugs found while reading
# fs_masterAnalysis.m/identify_grid_cell.m that are NOT reproduced here
# (a rotSpkPos(rotSpkPos==0)=nLpForShuf variable mix-up that's moot in
# 0-based indexing, and identify_grid_cell.m referencing an undefined
# `data.pos` that would have errored on every real MATLAB call).
# --------------------------------------------------------------------------
def test_brandon_gridness_matches_best_of_both_methods():
    from barrylab.analysis.grid import fast_gridness, spatial_autocorrelogram
    from barrylab.analysis.grid_ellipse import brandon_regularise

    def hex_grid_map(n=81, spacing=10.0, sigma=2.0, stretch_x=1.6):
        yy, xx = np.mgrid[0:n, 0:n]
        a1 = np.array([spacing * stretch_x, 0.0])
        a2 = np.array([spacing * 0.5 * stretch_x, spacing * np.sqrt(3) / 2])
        centre = np.array([n / 2.0, n / 2.0])
        m = np.zeros((n, n))
        for i in range(-5, 6):
            for j in range(-5, 6):
                c = i * a1 + j * a2 + centre
                m += np.exp(-((xx - c[0]) ** 2 + (yy - c[1]) ** 2) / (2 * sigma**2))
        return m

    sac = spatial_autocorrelogram(hex_grid_map())
    rsac1, rsac2 = brandon_regularise(sac)
    std1, exp1, _ = fast_gridness(rsac1)
    std2, exp2, _ = fast_gridness(rsac2)
    bg_std, bg_exp = _brandon_gridness(sac)
    assert bg_std == max(std1, std2)
    assert bg_exp == max(exp1, exp2)


def test_field_shuffle_significance_runs_end_to_end(example_trial):
    import barrylab as bl

    rec = bl.load_dacq(example_trial)
    result = field_shuffle_significance(rec, 2, 3, n_shuffles=4, rng=0)
    assert result["gridness"].shape == (4,)
    assert result["p95"].shape == (4, 2)
    assert result["gridness_from_shuffle"].shape == (4, 2, 4)
    assert np.isfinite(result["grid_scale"])
    # standard gridness (no correction) should match fast_gridness directly
    from barrylab.analysis.grid import fast_gridness

    sac = bl.spatial_autocorrelogram(bl.rate_map(rec, 2, 3, smoothing="gaussian", smooth_size=1.5))
    expected_std, expected_exp, _ = fast_gridness(sac)
    assert np.isclose(result["gridness"][0], expected_std)
    assert np.isclose(result["gridness"][2], expected_exp)


def test_field_shuffle_significance_keep_ratemaps(example_trial):
    import barrylab as bl

    rec = bl.load_dacq(example_trial)
    result = field_shuffle_significance(rec, 2, 3, n_shuffles=2, keep_ratemaps=True, rng=1)
    assert "ratemaps" in result
    assert result["ratemaps"]["real"].ndim == 2
    assert len(result["ratemaps"]["temporal_shuffle"]) == 2
    assert len(result["ratemaps"]["field_shuffle"]) == 2
    assert all(
        rm.shape == result["ratemaps"]["real"].shape
        for rm in result["ratemaps"]["temporal_shuffle"]
    )


def test_is_grid_cell_discriminates_example_cells(example_trial):
    import barrylab as bl

    rec = bl.load_dacq(example_trial)
    for cluster in (1, 3, 4, 5):
        result = is_grid_cell(rec, 2, cluster, n_shuffles=6, rng=0)
        assert result["is_grid"].shape == (4, 2)
        assert result["gridness"].shape == (4,)
        assert result["threshold"].shape == (4, 2)
