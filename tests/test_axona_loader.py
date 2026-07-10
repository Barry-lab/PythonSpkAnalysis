"""Golden-master tests for the DACQ loader (skip if example data absent)."""

import numpy as np


def test_load_dacq_basic(example_trial):
    import barrylab as bl

    rec = bl.load_dacq(example_trial)
    assert rec.n_pos == 60050
    assert rec.pos_sample_rate == 50.0
    assert rec.pixels_per_metre == 365.0
    assert sorted(rec.electrodes) == [1, 2, 3, 4, 5, 6, 7, 8]
    assert abs(rec.duration - 1201.0) < 1.0


def test_position_pipeline_matches_matlab(example_trial, golden_mat):
    """xy / direction / speed reproduce MATLAB to machine precision (midpoint)."""
    import barrylab as bl

    g = golden_mat["data"].pos
    rec = bl.load_dacq(example_trial, head_pos_weighting="midpoint")

    assert np.nanmax(np.abs(rec.xy - g.xy)) < 1e-9
    dd = np.abs(np.mod(rec.direction - g.dir + 180, 360) - 180)
    assert np.nanmax(dd) < 1e-8
    assert np.nanmax(np.abs(rec.speed - g.speed)) < 1e-9


def test_spike_mapping_matches_matlab(example_trial, golden_mat):
    import barrylab as bl

    rec = bl.load_dacq(example_trial, head_pos_weighting="midpoint")
    elec = rec.electrode(2)
    g_tet = next(t for t in golden_mat["data"].tetrode if int(t.id) == 2)

    # pos_sample: port is 0-based, MATLAB 1-based.
    assert np.array_equal(elec.pos_index + 1, np.asarray(g_tet.pos_sample))
    assert np.array_equal(elec.clusters, golden_mat["exactCut"])
    assert np.array_equal(elec.spike_pos_index(3) + 1, golden_mat["spikePosSamp"])


def test_binning_totals(example_trial, golden_mat):
    """Occupancy sums to duration; spikes sum to the cell's spike count."""
    import barrylab as bl

    rec = bl.load_dacq(example_trial, head_pos_weighting="midpoint")
    occ = bl.occupancy_map(rec)
    spk = bl.spike_map(rec, rec.electrode(2).spike_pos_index(3))
    assert abs(occ.sum() - rec.duration) < 1e-6
    assert spk.sum() == golden_mat["spikePosSamp"].size
    # rate map is a strong sanity correlation vs the (precursor-smoothed) fixture
    rm = bl.rate_map(rec, 2, 3)
    both = ~np.isnan(rm) & ~np.isnan(golden_mat["smthRm"])
    assert np.corrcoef(rm[both], golden_mat["smthRm"][both])[0, 1] > 0.9


def test_rate_map_smoothing_options(example_trial):
    """boxcar (default), gaussian and no-smoothing all produce valid maps."""
    import barrylab as bl

    rec = bl.load_dacq(example_trial)
    box = bl.rate_map(rec, 2, 3)                                   # default boxcar-5
    gauss = bl.rate_map(rec, 2, 3, smoothing="gaussian", smooth_size=2)
    raw = bl.rate_map(rec, 2, 3, smoothing=None)

    for rm in (box, gauss, raw):
        assert rm.shape == (62, 62)
        assert np.array_equal(np.isnan(rm), np.isnan(raw))  # same visited support
    # smoothing reduces bin-to-bin variance relative to the raw map
    assert np.nanstd(box) < np.nanstd(raw)
    assert np.nanstd(gauss) < np.nanstd(raw)
