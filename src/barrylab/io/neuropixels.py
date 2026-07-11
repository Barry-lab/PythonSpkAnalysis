"""Neuropixels front-end loader (EXPERIMENTAL, lightly tested).

Turns a *spike-sorted* Neuropixels recording plus a *separate position track*
into a :class:`~barrylab.core.recording.Recording`, so every spatial analysis
in :mod:`barrylab` (rate maps, gridness, head direction, decoding, ...) works
on Neuropixels data unchanged.

Two things a Neuropixels experiment gives you that this loader consumes:

1. **Sorted spikes.** Neuropixels raw data must first be spike-sorted (e.g.
   with Kilosort). This loader does NOT sort -- it reads *already-sorted*
   output. Point ``sorting`` at a Kilosort/phy output folder and it is read via
   `SpikeInterface <https://spikeinterface.readthedocs.io>`_
   (``pip install "barrylab[neuropixels]"``); or pass any SpikeInterface
   ``Sorting`` object, or in fact any object exposing ``get_unit_ids()``,
   ``get_unit_spike_train(unit_id)`` (spike sample indices) and a
   ``sampling_frequency``.

2. **Position.** A Neuropixels probe records neural data only -- it has NO idea
   where the animal is. Position comes from a *separate* behavioural source
   (a VR rig's position log, an overhead-camera tracker, DeepLabCut, ...),
   sampled on its own clock. You pass it in as ``positions`` (an ``(n_pos, 2)``
   array of ``x, y``) together with its ``pos_sample_rate``. The neural and
   behavioural clocks must already be aligned to a common ``t=0`` (typically
   via a shared sync/TTL pulse) -- this loader assumes they are.

.. warning::
   This front end is **experimental and far less tested than the DACQ loader**.
   It has no golden-master validation (unlike the bit-exact DACQ chain) and is
   exercised only against synthetic/mock sortings, because open Neuropixels
   spatial datasets are large (10s-100s of GB) and require GPU spike-sorting
   that this package doesn't do. Treat results with care, check them against a
   method you trust, and please improve/fix this loader -- it is a starting
   point, not a finished, validated pipeline. See the README.
"""

import numpy as np

from ..core.recording import Electrode, Recording
from . import _position


def _resolve_sorting(sorting, sorting_format="auto"):
    """Return ``(unit_ids, {unit_id: spike_frames}, sampling_frequency)``.

    ``sorting`` may be a path to a Kilosort/phy output folder (read via
    SpikeInterface), a SpikeInterface ``Sorting`` object, or any duck-typed
    object exposing ``get_unit_ids()``, ``get_unit_spike_train(unit_id)`` and
    ``sampling_frequency`` (attribute or ``get_sampling_frequency()``).
    """
    if isinstance(sorting, (str, bytes)) or hasattr(sorting, "__fspath__"):
        try:
            from spikeinterface import extractors as se  # noqa: PLC0415
        except ImportError as exc:  # pragma: no cover - env-dependent
            raise ImportError(
                "reading a Kilosort/phy folder needs SpikeInterface -- install "
                'the optional extra: pip install "barrylab[neuropixels]". '
                "Alternatively pass an already-loaded Sorting object."
            ) from exc
        path = str(sorting)
        readers = {
            "phy": se.read_phy,
            "kilosort": se.read_kilosort,
        }
        if sorting_format == "auto":
            # phy's reader also handles a bare Kilosort folder; try it first.
            sorting = readers["phy"](path)
        elif sorting_format in readers:
            sorting = readers[sorting_format](path)
        else:
            raise ValueError(
                "sorting_format must be 'auto', 'phy' or 'kilosort', got %r" % (sorting_format,)
            )

    unit_ids = list(sorting.get_unit_ids())
    if hasattr(sorting, "sampling_frequency"):
        fs = float(sorting.sampling_frequency)
    elif hasattr(sorting, "get_sampling_frequency"):
        fs = float(sorting.get_sampling_frequency())
    else:
        raise AttributeError("sorting object exposes no sampling_frequency")

    spike_trains = {uid: np.asarray(sorting.get_unit_spike_train(uid)) for uid in unit_ids}
    return unit_ids, spike_trains, fs


def load_neuropixels(
    sorting,
    positions,
    pos_sample_rate,
    pixels_per_metre,
    window=None,
    direction=None,
    speed=None,
    name="neuropixels",
    probe_id=1,
    unit_ids=None,
    sorting_format="auto",
    metadata=None,
):
    """Build a :class:`~barrylab.core.recording.Recording` from sorted
    Neuropixels spikes and a behavioural position track.

    Parameters
    ----------
    sorting
        A spike-sorting result: a path to a Kilosort/phy output folder, a
        SpikeInterface ``Sorting`` object, or any object exposing
        ``get_unit_ids()`` / ``get_unit_spike_train(unit_id)`` /
        ``sampling_frequency``. Spike trains are spike *sample indices* on the
        neural (Neuropixels) clock.
    positions : ndarray (n_pos, 2)
        Behavioural position (``x, y``) sampled at ``pos_sample_rate``, on a
        clock already aligned to the neural clock's ``t=0``. Units are whatever
        ``pixels_per_metre`` is expressed in (pixels if from a camera; if your
        track is already in cm, pass it in cm and set ``pixels_per_metre=100``).
    pos_sample_rate : float
        Position sampling rate in Hz.
    pixels_per_metre : float
        Spatial scale, for converting position units to cm.
    window : (float, float, float, float), optional
        ``(min_x, max_x, min_y, max_y)`` tracked window. Defaults to the
        bounding box of ``positions`` (ignoring NaNs).
    direction : ndarray (n_pos,), optional
        Head direction in degrees. If ``None``, it is derived from the
        frame-to-frame displacement of ``positions`` (the 1-LED / movement
        heading, degrees anticlockwise from +x) -- note this is movement
        direction, a proxy for head direction, not a true 2-marker heading.
    speed : ndarray (n_pos,), optional
        Running speed in cm/s. If ``None``, derived from ``positions``.
    name : str
        Session name.
    probe_id : int
        The id the probe is stored under in ``Recording.electrodes`` (reach a
        unit with ``rate_map(rec, electrode=probe_id, cluster=<unit_id>)``).
    unit_ids : sequence, optional
        Restrict to these sorted units (e.g. only "good" clusters). Defaults to
        all units in the sorting.
    sorting_format : {"auto", "phy", "kilosort"}
        Which SpikeInterface reader to use when ``sorting`` is a path.
    metadata : dict, optional
        Extra provenance stored on the ``Recording``.

    Returns
    -------
    Recording
        With one :class:`~barrylab.core.recording.Electrode` (``probe_id``)
        whose ``clusters`` are the sorted unit ids and whose ``pos_index`` maps
        each spike to the position sample it co-occurred with.
    """
    all_units, spike_trains, fs = _resolve_sorting(sorting, sorting_format)
    if unit_ids is not None:
        keep = set(unit_ids)
        all_units = [u for u in all_units if u in keep]

    positions = np.asarray(positions, dtype=float)
    if positions.ndim != 2 or positions.shape[1] != 2:
        raise ValueError("positions must be (n_pos, 2)")
    n_pos = positions.shape[0]
    pos_times = np.arange(n_pos) / pos_sample_rate

    if window is None:
        window = (
            float(np.nanmin(positions[:, 0])),
            float(np.nanmax(positions[:, 0])),
            float(np.nanmin(positions[:, 1])),
            float(np.nanmax(positions[:, 1])),
        )

    if direction is None:
        direction = _position._heading_from_displacement(positions, pos_sample_rate)
    else:
        direction = np.asarray(direction, dtype=float)

    if speed is None:
        step = np.hypot(np.diff(positions[:, 0]), np.diff(positions[:, 1]))
        speed = np.empty(n_pos)
        speed[: n_pos - 1] = step
        speed[n_pos - 1] = speed[n_pos - 2] if n_pos > 1 else 0.0
        speed = speed * (100.0 * pos_sample_rate / pixels_per_metre)  # px/sample -> cm/s
    else:
        speed = np.asarray(speed, dtype=float)

    # Flatten all units into one Electrode: per-spike time (s), cluster id, and
    # the position sample each spike co-occurred with. Spikes are sorted by
    # time (units interleaved) so the arrays read naturally; ordering does not
    # affect any downstream binning (which selects by cluster).
    times_list, clusters_list = [], []
    for uid in all_units:
        frames = spike_trains[uid]
        if frames.size == 0:
            continue
        times_list.append(frames.astype(np.float64) / fs)
        clusters_list.append(np.full(frames.size, uid))
    if times_list:
        spike_times = np.concatenate(times_list)
        clusters = np.concatenate(clusters_list)
        order = np.argsort(spike_times, kind="stable")
        spike_times, clusters = spike_times[order], clusters[order]
    else:
        spike_times = np.array([], dtype=np.float64)
        clusters = np.array([], dtype=object)

    # Position sample co-occurring with each spike (0-based), clamped to range.
    pos_index = np.clip(np.floor(spike_times * pos_sample_rate).astype(np.int64), 0, n_pos - 1)

    electrode = Electrode(
        identifier=probe_id,
        spike_times=spike_times,
        pos_index=pos_index,
        clusters=clusters,
        metadata={"source": "neuropixels", "sampling_frequency": fs, "n_units": len(all_units)},
    )

    meta = {"format": "neuropixels", "sorting_sample_rate": fs}
    if metadata:
        meta.update(metadata)

    return Recording(
        name=name,
        xy=positions,
        direction=direction,
        speed=speed,
        pos_times=pos_times,
        pos_sample_rate=float(pos_sample_rate),
        pixels_per_metre=float(pixels_per_metre),
        window=window,
        electrodes={probe_id: electrode},
        metadata=meta,
    )
