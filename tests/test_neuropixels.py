"""Tests for the (experimental) Neuropixels front-end loader.

No real Neuropixels data or SpikeInterface is needed: a ``Sorting`` is just an
object exposing ``get_unit_ids`` / ``get_unit_spike_train`` / a sampling rate,
so a small duck-typed mock stands in for one. These check that the loader
correctly turns (sorted spikes + a position track) into a ``Recording`` that
the existing spatial analyses run on unchanged -- which is the whole point of
the one-backend-many-frontends design.
"""

import numpy as np

import barrylab as bl
from barrylab.io.neuropixels import load_neuropixels


class _MockSorting:
    """Minimal stand-in for a SpikeInterface Sorting object."""

    def __init__(self, spike_trains, sampling_frequency):
        self._trains = spike_trains  # {unit_id: np.array of spike frames}
        self.sampling_frequency = float(sampling_frequency)

    def get_unit_ids(self):
        return list(self._trains.keys())

    def get_unit_spike_train(self, unit_id):
        return self._trains[unit_id]


def _synthetic_place_cell_sorting(fs=30000.0, n_pos=6000, pos_rate=50.0):
    """A place cell that fires in one corner, plus a silent-ish second unit,
    on a synthetic random-walk position track."""
    rng = np.random.RandomState(0)
    # random-walk position in a 100x100 (px) box
    steps = rng.normal(0, 2.0, size=(n_pos, 2))
    xy = np.clip(np.cumsum(steps, axis=0) % 100.0, 0, 100.0)

    duration_s = n_pos / pos_rate
    # unit 7: fires when the animal is near (20, 20)
    near = np.hypot(xy[:, 0] - 20, xy[:, 1] - 20) < 12
    fire_times = []
    for i in np.flatnonzero(near):
        t0 = i / pos_rate
        fire_times.extend(t0 + rng.uniform(0, 1.0 / pos_rate, size=rng.poisson(3)))
    unit7 = np.sort(np.array(fire_times)) * fs  # -> sample frames
    # unit 3: sparse, uniform in time (no place field)
    unit3 = np.sort(rng.uniform(0, duration_s, size=200)) * fs
    trains = {7: unit7.astype(np.int64), 3: unit3.astype(np.int64)}
    return _MockSorting(trains, fs), xy, pos_rate


def test_load_neuropixels_builds_a_valid_recording():
    sorting, xy, pos_rate = _synthetic_place_cell_sorting()
    rec = load_neuropixels(sorting, positions=xy, pos_sample_rate=pos_rate, pixels_per_metre=100.0)

    assert isinstance(rec, bl.Recording)
    assert rec.n_pos == xy.shape[0]
    assert rec.pos_sample_rate == pos_rate
    assert set(rec.electrodes) == {1}  # default probe id
    elec = rec.electrode(1)
    assert set(elec.cluster_ids().tolist()) == {3, 7}
    # every spike maps to a valid position sample
    assert elec.pos_index.min() >= 0
    assert elec.pos_index.max() < rec.n_pos
    # direction and speed were derived from the track
    assert rec.direction.shape == (rec.n_pos,)
    assert rec.speed.shape == (rec.n_pos,)


def test_neuropixels_recording_works_with_rate_map_and_spatial_info():
    sorting, xy, pos_rate = _synthetic_place_cell_sorting()
    rec = load_neuropixels(sorting, positions=xy, pos_sample_rate=pos_rate, pixels_per_metre=100.0)

    # The whole point: existing analyses run unchanged on a NPX Recording.
    rmap = bl.rate_map(rec, electrode=1, cluster=7, bin_cm=4.0)
    occ = bl.occupancy_map(rec, bin_cm=4.0)
    assert rmap.shape == occ.shape
    bits_per_spike, _ = bl.spatial_information(rmap, occ)
    assert np.isfinite(bits_per_spike)

    # the place cell (unit 7) should carry more spatial information than the
    # spatially-uniform unit 3
    rmap3 = bl.rate_map(rec, electrode=1, cluster=3, bin_cm=4.0)
    bps3, _ = bl.spatial_information(rmap3, occ)
    assert bits_per_spike > bps3


def test_load_neuropixels_restricts_to_requested_units():
    sorting, xy, pos_rate = _synthetic_place_cell_sorting()
    rec = load_neuropixels(
        sorting, positions=xy, pos_sample_rate=pos_rate, pixels_per_metre=100.0, unit_ids=[7]
    )
    assert rec.electrode(1).cluster_ids().tolist() == [7]


def test_load_neuropixels_accepts_explicit_window_and_probe_id():
    sorting, xy, pos_rate = _synthetic_place_cell_sorting()
    rec = load_neuropixels(
        sorting,
        positions=xy,
        pos_sample_rate=pos_rate,
        pixels_per_metre=100.0,
        window=(0.0, 100.0, 0.0, 100.0),
        probe_id=42,
    )
    assert set(rec.electrodes) == {42}
    assert rec.window == (0.0, 100.0, 0.0, 100.0)
    assert rec.window_extent == (100.0, 100.0)


def test_load_neuropixels_rejects_bad_positions():
    sorting, _, pos_rate = _synthetic_place_cell_sorting()
    import pytest

    with pytest.raises(ValueError):
        load_neuropixels(sorting, positions=np.zeros((10, 3)), pos_sample_rate=pos_rate,
                         pixels_per_metre=100.0)


def test_load_neuropixels_path_without_spikeinterface_raises_clear_error(tmp_path):
    # Passing a path triggers the SpikeInterface import; if it's absent the
    # error must name the optional extra rather than a bare ImportError.
    import pytest

    try:
        import spikeinterface  # noqa: F401
        has_si = True
    except ImportError:
        has_si = False
    if has_si:
        import pytest as _pytest
        _pytest.skip("spikeinterface is installed; can't test the missing-dep path")
    with pytest.raises(ImportError, match="barrylab\\[neuropixels\\]"):
        load_neuropixels(str(tmp_path), positions=np.zeros((5, 2)), pos_sample_rate=50.0,
                         pixels_per_metre=100.0)


def test_load_neuropixels_puts_positions_in_window_frame():
    """Positions away from the camera origin must not pile into the edge bins.

    Before the fix the default window was the positions' bounding box but xy
    was not shifted by its minimum, so every sample above (max - min) was
    clamped into the last row/column of the occupancy map.
    """
    sorting, xy, pos_rate = _synthetic_place_cell_sorting()
    offset = np.array([150.0, 80.0])
    rec_origin = load_neuropixels(
        sorting, positions=xy, pos_sample_rate=pos_rate, pixels_per_metre=100.0
    )
    rec_offset = load_neuropixels(
        sorting, positions=xy + offset, pos_sample_rate=pos_rate, pixels_per_metre=100.0
    )
    np.testing.assert_allclose(rec_offset.xy, rec_origin.xy)
    np.testing.assert_allclose(
        bl.analysis.maps.occupancy_map(rec_offset), bl.analysis.maps.occupancy_map(rec_origin)
    )
    np.testing.assert_allclose(
        bl.rate_map(rec_offset, electrode=1, cluster=7),
        bl.rate_map(rec_origin, electrode=1, cluster=7),
        equal_nan=True,
    )
    assert rec_offset.window[0] == xy[:, 0].min() + offset[0]
