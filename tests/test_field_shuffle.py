"""Tests for watershed field segmentation (field-shuffle grid-cell test, v0.7 partial)."""

import numpy as np
import pytest

skimage = pytest.importorskip("skimage", reason="scikit-image not installed (barrylab[fieldshuffle])")

from barrylab.analysis.field_shuffle import (  # noqa: E402
    _fix_watershed,
    field_bin_locations,
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
