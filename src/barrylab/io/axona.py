"""Axona / DACQ front-end loader.

``load_dacq`` reads a raw DACQ trial (a set of files sharing one basename) and
returns a :class:`~barrylab.core.recording.Recording`. It ports the master
``read_DACQ.m`` and the per-file readers ``read_pos_file.m``,
``read_tetrode_file.m`` and ``read_cut_file.m``.

The public entry point is :func:`load_dacq`; the ``_read_*`` helpers are exposed
for validation and testing but are not part of the stable API.
"""

import os
from typing import Optional

import numpy as np

from ..core.recording import Electrode, Recording
from ..params import POS_MAX_SPEED_MS, POS_SMOOTH_BOXCAR_S
from . import _position
from ._axona_files import big_endian_view, key_value, read_axona_file

# Tetrode files store 4 channels x (4-byte timestamp + 50 int8 waveform samples)
# = 54 bytes/channel = 216 bytes/spike. Timestamps are identical across the 4
# channels, so only the first is read.
_BYTES_PER_SPIKE = 216
_SPIKE_TIMESTAMP_BYTES = 4


def _read_pos(path):
    """Read a ``.pos`` file into raw per-LED positions.

    Returns ``(led_pos, led_pix, header)`` where ``led_pos`` is
    ``(n_pos, n_leds, 2)`` in pixels (untracked flagged as NaN) and ``led_pix``
    is ``(n_pos, n_leds)`` pixel counts. Port of ``read_pos_file.m``.
    """
    header, data = read_axona_file(path, has_binary=True)
    n_pos = int(key_value(header, "num_pos_samples", "num"))
    n_colours = int(key_value(header, "num_colours", "num"))

    # Each record: 4-byte timestamp + n_colours * (x, y) as big-endian uint16.
    record_len = _SPIKE_TIMESTAMP_BYTES + n_colours * 2 * 2
    buf = np.frombuffer(data, dtype=np.uint8)
    buf = buf[: n_pos * record_len].reshape(n_pos, record_len)
    payload = buf[:, _SPIKE_TIMESTAMP_BYTES:]  # drop the (unused) timestamp
    shorts = big_endian_view(payload, np.uint16)  # (n_pos, n_colours * 2)

    # Current 2-spot format is 't,x1,y1,x2,y2,numpix1,numpix2'.
    led_pos = np.zeros((n_pos, 2, 2))
    led_pos[:, 0, 0], led_pos[:, 0, 1] = shorts[:, 0], shorts[:, 1]
    led_pos[:, 1, 0], led_pos[:, 1, 1] = shorts[:, 2], shorts[:, 3]
    led_pix = np.zeros((n_pos, 2))
    led_pix[:, 0], led_pix[:, 1] = shorts[:, 4], shorts[:, 5]

    led_pos, led_pix = _drop_unused_leds(led_pos, led_pix)
    led_pos[led_pos == 1023] = np.nan
    led_pix[led_pix == 1023] = np.nan
    return led_pos, led_pix, header


def _drop_unused_leds(led_pos, led_pix):
    """Remove LED columns that are entirely absent (all 0 or all 1023)."""
    bad = (led_pos == 0) | (led_pos == 1023)
    used = ~np.all(bad.reshape(led_pos.shape[0], led_pos.shape[1], 2).all(axis=2), axis=0)
    if not used.all():
        led_pos = led_pos[:, used, :]
        led_pix = led_pix[:, used]
    return led_pos, led_pix


def _read_tetrode(path):
    """Read spike times (seconds) from a tetrode file. Port of ``read_tetrode_file.m``.

    Returns ``(spike_times, header)``. Waveforms are not loaded (not needed for
    spatial analyses); this mirrors the fast ``loadWave=0`` MATLAB path.
    """
    header, data = read_axona_file(path, has_binary=True)
    timebase = key_value(header, "timebase", "num")
    n_spikes = int(key_value(header, "num_spikes", "num"))

    buf = np.frombuffer(data, dtype=np.uint8)
    buf = buf[: n_spikes * _BYTES_PER_SPIKE].reshape(n_spikes, _BYTES_PER_SPIKE)
    ts_raw = big_endian_view(buf[:, :_SPIKE_TIMESTAMP_BYTES], np.int32).reshape(-1)
    spike_times = ts_raw.astype(np.float64) / timebase
    return spike_times, header


def _read_cut(path):
    """Read cluster ids per spike from a ``.cut`` (Tint) or ``.clu`` file.

    Port of ``read_cut_file.m``: after the ``Exact_cut_for: <file> spikes: <n>``
    line, ``n`` cluster ids follow (possibly across several lines). A bare
    ``.clu`` file is just the list of ids.
    """
    with open(path, "r") as fh:
        tokens = fh.read().split()
    if "Exact_cut_for:" in tokens:
        start = tokens.index("spikes:", tokens.index("Exact_cut_for:"))
        n_spikes = int(tokens[start + 1])
        ids = tokens[start + 2: start + 2 + n_spikes]
        return np.array([int(x) for x in ids], dtype=int)
    # Bare .clu: an optional leading count then one id per spike.
    ids = [int(t) for t in tokens if t.lstrip("-").isdigit()]
    return np.array(ids, dtype=int)


def _tetrodes_in_use(settings):
    """Tetrode numbers (1-16) marked as collected in the set file."""
    used = []
    for t in range(1, 17):
        mask = key_value(settings, "collectMask_%d" % t, "num", exact=True)
        if mask:
            used.append(t)
    return used


def load_dacq(path, load_position=True,
              max_speed=POS_MAX_SPEED_MS, box_car=POS_SMOOTH_BOXCAR_S,
              head_pos_weighting="tracking"):
    """Load a raw Axona/DACQ trial into a :class:`Recording`.

    Parameters
    ----------
    path
        Path to any file of the trial, or the shared basename (extension is
        stripped). E.g. ``".../83812141"`` or ``".../83812141.set"``.
    load_position
        If ``True`` read and post-process the ``.pos`` file.
    max_speed, box_car
        Position-processing parameters (see :mod:`barrylab.params`).
    head_pos_weighting : {"tracking", "midpoint"}
        Two-LED head-position convention; see
        :func:`barrylab.io._position.postprocess_pos_data`.

    Returns
    -------
    Recording
    """
    stem = os.path.splitext(path)[0]
    name = os.path.basename(stem)
    settings, _ = read_axona_file(stem + ".set", has_binary=False)

    rec = Recording(name=name, metadata={"source": stem, "format": "dacq"})

    # --- position -------------------------------------------------------
    sample_rate = None
    if load_position and os.path.exists(stem + ".pos"):
        led_pos, led_pix, pos_header = _read_pos(stem + ".pos")
        pos = _position.postprocess_pos_data(
            led_pos, led_pix, pos_header, settings,
            max_speed=max_speed, box_car=box_car,
            head_pos_weighting=head_pos_weighting,
        )
        sample_rate = key_value(pos_header, "sample_rate", "num")
        rec.xy = pos["xy"]
        rec.direction = pos["direction"]
        rec.speed = pos["speed"]
        rec.pos_times = pos["times"]
        rec.pos_sample_rate = sample_rate
        rec.pixels_per_metre = key_value(pos_header, "pixels_per_metre", "num")
        rec.window = (
            key_value(pos_header, "window_min_x", "num"),
            key_value(pos_header, "window_max_x", "num"),
            key_value(pos_header, "window_min_y", "num"),
            key_value(pos_header, "window_max_y", "num"),
        )
        rec.metadata["pos_header"] = pos_header

    # --- tetrodes (+ cut) ----------------------------------------------
    for t in _tetrodes_in_use(settings):
        tet_path = "%s.%d" % (stem, t)
        if not os.path.exists(tet_path):
            continue
        spike_times, tet_header = _read_tetrode(tet_path)

        pos_index = None
        if sample_rate is not None:
            # Position sample co-occurring with each spike. MATLAB uses a
            # 1-based ceil; we keep the same rounding then convert to 0-based.
            one_based = np.ceil(spike_times * sample_rate).astype(np.int64)
            in_range = one_based <= rec.n_pos
            spike_times = spike_times[in_range]
            pos_index = one_based[in_range] - 1

        cut_path = "%s_%d.cut" % (stem, t)
        clusters = _read_cut(cut_path) if os.path.exists(cut_path) else None

        rec.electrodes[t] = Electrode(
            identifier=t,
            spike_times=spike_times,
            pos_index=pos_index,
            clusters=clusters,
            metadata={"header": tet_header},
        )

    return rec
