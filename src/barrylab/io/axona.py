"""Axona / DACQ front-end loader.

``load_dacq`` reads a raw DACQ trial (a set of files sharing one basename) and
returns a :class:`~barrylab.core.recording.Recording`. It ports the master
``read_DACQ.m`` and the per-file readers ``read_pos_file.m``,
``read_tetrode_file.m`` and ``read_cut_file.m``. Also ports the standalone
auxiliary readers/writer :func:`read_inp` (``read_INP_file.m``), :func:`read_stm`
(``read_stm_file.m``) and :func:`write_cut` (``write_cut_file.m``).

The public entry point is :func:`load_dacq`; the ``_read_*`` helpers are exposed
for validation and testing but are not part of the stable API.
"""

import os
import warnings

import numpy as np

from ..core.recording import Electrode, Recording
from ..params import POS_MAX_SPEED_MS, POS_SMOOTH_BOXCAR_S
from . import _position
from ._axona_files import (
    big_endian_view,
    key_value,
    little_endian_view,
    read_axona_file,
)

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
    with open(path) as fh:
        tokens = fh.read().split()
    if "Exact_cut_for:" in tokens:
        start = tokens.index("spikes:", tokens.index("Exact_cut_for:"))
        n_spikes = int(tokens[start + 1])
        ids = tokens[start + 2 : start + 2 + n_spikes]
        return np.array([int(x) for x in ids], dtype=int)
    # Bare .clu: an optional leading count then one id per spike.
    ids = [int(t) for t in tokens if t.lstrip("-").isdigit()]
    return np.array(ids, dtype=int)


def read_inp(path):
    """Read a ``.inp`` file: DACQ digital-timestamp events (button presses, DIO).

    Port of ``read_INP_file.m``. Each event is a fixed 7-byte record: a 4-byte
    big-endian timestamp, a 1-byte ASCII event type (``'I'``/``'O'``/``'K'`` for
    digital input / digital output / keypress respectively), and a 2-byte
    encoded channel bitmask. The file logs both the onset *and* offset of
    digital events, so a positive-channel ``'I'``/``'O'`` event is typically
    followed by a channel-0 "off" event (see the MATLAB source's own header
    comment for the full field semantics); keypress (``'K'``) channels are not
    decoded (the two-number DACQ keypress code is left as NaN), matching the
    original.

    No example ``.inp`` file ships with ``Universal_Matlab``'s ``Example_data``,
    so this is validated structurally (a synthetic file built byte-for-byte
    from the documented header/record layout) rather than against real data --
    see ``tests/test_inp_stm_cut.py`` and ``FINDINGS.md``.

    Parameters
    ----------
    path
        Full path to the ``.inp`` file.

    Returns
    -------
    (time, type, chan)
        ``time`` -- ndarray (n_events,) of event times in seconds.
        ``type`` -- ndarray (n_events,) of single-character event-type strings.
        ``chan`` -- ndarray (n_events,) of channel numbers (1-16 for digital
        events, 0 for a digital "off" event, NaN for keypresses).
    """
    header, data = read_axona_file(path, has_binary=True)
    timebase = key_value(header, "timebase", "num")

    # Assumes the binary segment is an exact multiple of 7 bytes (one event
    # each); any trailing partial bytes are silently dropped, matching the
    # truncating convention already used by _read_pos/_read_tetrode above.
    n_events = len(data) // 7
    buf = np.frombuffer(data, dtype=np.uint8)[: n_events * 7].reshape(n_events, 7)

    time = big_endian_view(buf[:, :4], np.int32).reshape(-1).astype(np.float64) / timebase
    event_type = np.array([chr(b) for b in buf[:, 4]])

    # Channel: byte 6 (0-indexed column 0) encodes channels 9-16 (bit 0 = ch9),
    # byte 7 (column 1) encodes channels 1-8 (bit 0 = ch1). At most one bit is
    # expected to be set per byte, so log2 of the byte value recovers its bit
    # position; a byte of 0 (no bit set) -> log2(0) = -inf -> nan, and the
    # nansum below then contributes 0 for that byte -- exactly reproducing the
    # "off" (channel 0) signal the MATLAB header comment describes, and MATLAB's
    # own nansum-of-two-NaNs-is-0 behaviour when neither byte has a bit set.
    chan_bytes = buf[:, 5:7].astype(np.float64)
    with np.errstate(divide="ignore"):
        bits = np.log2(chan_bytes)
    bits[np.isinf(bits)] = np.nan
    bits[:, 0] += 9  # byte 6: channels 9-16
    bits[:, 1] += 1  # byte 7: channels 1-8
    chan = np.nansum(bits, axis=1)
    chan[event_type == "K"] = np.nan

    return time, event_type, chan


def read_stm(path):
    """Read a ``.stm`` file: DACQ stimulator event timestamps (via the DIO port).

    Port of ``read_stm_file.m``. The MATLAB source hard-codes a fixed 11-line
    text header followed directly by a big-endian ``uint32`` binary stream (one
    timestamp per 4 bytes), ending just before a trailing ``data_end`` marker.

    Ported using the same robust ``data_start``/``data_end`` marker search
    already used by every other reader in this module, rather than
    transcribing the original's fragile fixed-line-count parsing (which only
    works if the header always has exactly 11 lines before the marker) --
    the identical design choice already made for ``.egf`` (see
    :func:`_read_egf`'s docstring: ``getegf.m``'s hard-coded "skip exactly 10
    header lines" approach was replaced with the marker search for the same
    reason).

    No example ``.stm`` file ships with ``Universal_Matlab``'s ``Example_data``,
    so this is validated structurally rather than against real data -- see
    ``tests/test_inp_stm_cut.py`` and ``FINDINGS.md``.

    Parameters
    ----------
    path
        Full path to the ``.stm`` file.

    Returns
    -------
    (stim_time, header)
        ``stim_time`` -- ndarray of stimulator event times in seconds.
        ``header`` -- the parsed header (list of ``(key, value)`` pairs).
    """
    header, data = read_axona_file(path, has_binary=True)
    timebase = key_value(header, "timebase", "num")

    n_values = len(data) // 4
    buf = np.frombuffer(data, dtype=np.uint8)[: n_values * 4]
    stim_time = big_endian_view(buf, np.uint32).reshape(-1).astype(np.float64) / timebase
    return stim_time, header


def write_cut(spike_assignment, dir_path, file_stem, tetrode_no, max_cluster_warn=30):
    """Write cluster assignments to a Tint-compatible ``.cut`` file.

    Port of ``write_cut_file.m``. Writes (overwriting without warning, matching
    the original) ``<dir_path>/<file_stem>_<tetrode_no>.cut`` -- the filename
    convention Tint requires to associate a cut file with its tetrode file.
    Cluster 0 is reserved for "junk" (unlike KlustaKwik's ``.clu`` convention).

    **Two MATLAB bugs found while reading the source, both fixed rather than
    faithfully reproduced -- see ``FINDINGS.md`` for the full detail:**

    1. ``write_cut_file.m``'s data-header line hard-codes the literal string
       ``Exact_cut_for: 1500728a`` -- one specific historic session's name --
       completely ignoring its own ``fileStem`` argument (a clear copy-paste
       bug, not an intentional convention: a self-consistent cut-file writer
       must embed the filename it was actually given). Fixed: this port writes
       the real ``file_stem``.
    2. A stale-loop-variable bug silently drops the *entire* spike-assignment
       body for any cut file with <=25 total spikes. MATLAB's data-row loop
       (``for n=0:ceil(nSpikes/25)-2``) reuses the loop variable ``n`` from an
       earlier, unrelated loop (writing the per-cluster header blocks); when
       the data-row loop's range is empty (true whenever ``nSpikes<=25``),
       MATLAB does not execute the loop body at all and leaves ``n`` at its
       *stale* value from the cluster-header loop, so the subsequent
       ``startLastRow=25*(n+1)+1`` computes a nonsensical, far-too-large
       starting index -- the final tail loop then also finds its range empty
       and writes *no* spike-assignment values at all, producing a cut file
       with a correct header but a completely empty (and, worse, silently so)
       body. This is a genuine, previously-undocumented bug (distinct from the
       one flagged in ``PORTING_PLAN.md``), confirmed by tracing the loop
       variable's scope by hand. Not preserved: this is silent data loss with
       no well-defined behaviour worth keeping, and faithfully reproducing it
       would break the very round-trip validation this port is checked with
       (a small synthetic cluster-assignment vector -- the natural choice for
       a unit test -- has well under 25 spikes). Fixed by computing the number
       of full 25-wide rows directly (``max(0, ceil(n_spikes/25) - 1)``)
       instead of relying on a loop variable's post-loop value.

    Parameters
    ----------
    spike_assignment : array-like of int
        One integer cluster id per spike (0 = junk).
    dir_path : str
        Directory to write into.
    file_stem : str
        Session filename stem (matches the ``.set``/tetrode files' basename);
        embedded in the ``Exact_cut_for:`` header line.
    tetrode_no : int
        Tetrode number this cut file corresponds to.
    max_cluster_warn : int, optional
        Warn (not raise, matching the original) if the maximum cluster id
        exceeds this (default 30, matching ``write_cut_file.m``'s ``maxCluster``).

    Returns
    -------
    str
        The full path written to.
    """
    spike_assignment = np.asarray(spike_assignment, dtype=int).ravel()
    if spike_assignment.size == 0:
        raise ValueError("spike_assignment must not be empty")
    n_spikes = spike_assignment.size
    max_cluster = int(spike_assignment.max())

    if max_cluster > max_cluster_warn:
        warnings.warn(
            "Cluster assignments exceed the max allowed (%d) -- might not be "
            "compatible with Tint." % max_cluster_warn,
            stacklevel=2,
        )

    full_path = os.path.join(dir_path, f"{file_stem}_{tetrode_no}.cut")
    with open(full_path, "w", newline="") as fh:
        fh.write(f"n_clusters: {max_cluster + 1}\r\n")
        fh.write("n_channels: 4\r\n")
        fh.write("n_params: 2\r\n")
        fh.write("times_used_in_Vt:   15   12   12   16\r\n")
        for n in range(max_cluster + 1):
            fh.write(f" cluster: {n:2d} centre:    0    0    0    0    0    0    0    0\r\n")
            fh.write("                min:    0    0    0    0    0    0    0    0\r\n")
            fh.write("                max:    256  256  256  256  256  256  256  256\r\n")
        fh.write("\r\n")

        fh.write(f"Exact_cut_for: {file_stem} spikes: {n_spikes:6d}\r\n")

        # Body: rows of 25 space-padded cluster ids, then any remaining tail
        # written with no full row and (matching the original exactly) no
        # trailing line terminator after the very last value.
        n_full_rows = max(0, -(-n_spikes // 25) - 1)  # ceil(n_spikes / 25) - 1
        for row in range(n_full_rows):
            chunk = spike_assignment[row * 25 : row * 25 + 25]
            fh.write("".join(f"  {v}" for v in chunk) + "\r\n")
        tail = spike_assignment[n_full_rows * 25 :]
        fh.write("".join(f"  {v}" for v in tail))

    return full_path


def _tetrodes_in_use(settings):
    """Tetrode numbers (1-16) marked as collected in the set file."""
    used = []
    for t in range(1, 17):
        mask = key_value(settings, f"collectMask_{t}", "num", exact=True)
        if mask:
            used.append(t)
    return used


def _read_eeg(path):
    """Read a ``.eeg`` file (250 Hz, signed 8-bit). Port of ``read_eeg_file.m``.

    Returns ``(samples, header)``; ``samples`` are raw signed int8 counts
    (not yet volts -- see :func:`eeg_convert_to_volts`).
    """
    header, data = read_axona_file(path, has_binary=True)
    n_samples = int(key_value(header, "num_EEG_samples", "num"))
    samples = np.frombuffer(data, dtype=np.int8)[:n_samples].astype(np.float64)
    return samples, header


def _read_egf(path):
    """Read a ``.egf`` file (4.8 kHz, signed 16-bit, little-endian).

    Port of ``getegf.m``, minus its optional decimation step (kept out of
    this low-level reader; see :mod:`barrylab.analysis.lfp` for an explicit,
    clearly-caveated resampling helper -- MATLAB's ``decimate()`` uses a
    specific IIR filter design that is not guaranteed to reproduce bit-exact
    in another language/library, so the raw undecimated signal is the
    faithful, unambiguous thing to return here). Unlike every other DACQ
    file, ``.egf`` is little-endian.

    Returns ``(samples, header)``; ``samples`` are raw signed int16 counts.
    """
    header, data = read_axona_file(path, has_binary=True)
    n_samples = int(key_value(header, "num_EGF_samples", "num"))
    buf = np.frombuffer(data, dtype=np.uint8)
    samples = little_endian_view(buf, np.int16).reshape(-1)[:n_samples].astype(np.float64)
    return samples, header


def eeg_convert_to_volts(samples, eeg_channel, settings, kind="eeg"):
    """Convert raw DACQ EEG/EGF counts to volts. Port of ``eeg_convert_to_volts.m``.

    Parameters
    ----------
    samples : ndarray
        Raw signed counts from :func:`_read_eeg` (``kind="eeg"``, int8) or
        :func:`_read_egf` (``kind="egf"``, int16).
    eeg_channel : int
        The 1-based EEG channel number (matches ``EEG_ch_N``/``saveEEG_ch_N``
        in the ``.set`` file).
    settings : Header
        The parsed ``.set`` header.
    kind : {"eeg", "egf"}
        ``"eeg"`` (default) for 8-bit 250 Hz data, ``"egf"`` for 16-bit 4.8 kHz.

    Returns
    -------
    ndarray of volts.
    """
    adc_fullscale = key_value(settings, "ADC_fullscale_mv", "num")
    conv_factor = {1500.0: 1.5, 3680.0: 3.68, 3000.0: 3.00}.get(adc_fullscale)
    if conv_factor is None:
        raise ValueError(f"unrecognised ADC_fullscale_mv {adc_fullscale!r} in settings")

    samples = np.asarray(samples, dtype=np.float64)
    if kind == "eeg":
        samples = samples / 128.0
    elif kind == "egf":
        samples = samples / 32768.0
    else:
        raise ValueError("kind must be 'eeg' or 'egf'")

    # eeg_channel is 1-based; the .set file's gain_ch_N keys are 0-based.
    gain = key_value(settings, f"gain_ch_{eeg_channel - 1}", "num", exact=True)
    return (conv_factor / gain) * samples


def _active_eeg_channels(settings, max_channels=16):
    """EEG channel numbers (1-based) flagged as saved in the ``.set`` file."""
    active = []
    for n in range(1, max_channels + 1):
        flag = key_value(settings, f"saveEEG_ch_{n}", "num", exact=True)
        if flag:
            active.append(n)
    return active


def _load_lfp_channel(stem, channel, settings):
    """Load one EEG channel's data (prefers the high-rate .egf if present)."""
    suffix = "" if channel == 1 else str(channel)
    egf_path = f"{stem}.egf{suffix}"
    eeg_path = f"{stem}.eeg{suffix}"

    if os.path.exists(egf_path):
        raw, header = _read_egf(egf_path)
        kind = "egf"
    elif os.path.exists(eeg_path):
        raw, header = _read_eeg(eeg_path)
        kind = "eeg"
    else:
        return None

    volts = eeg_convert_to_volts(raw, channel, settings, kind=kind)
    return {
        "samples": volts,
        "sample_rate": key_value(header, "sample_rate", "num"),
        "channel": channel,
        "kind": kind,
    }


def load_dacq(
    path,
    load_position=True,
    load_lfp=False,
    max_speed=POS_MAX_SPEED_MS,
    box_car=POS_SMOOTH_BOXCAR_S,
    head_pos_weighting="tracking",
):
    """Load a raw Axona/DACQ trial into a :class:`Recording`.

    Parameters
    ----------
    path
        Path to any file of the trial, or the shared basename (extension is
        stripped). E.g. ``".../83812141"`` or ``".../83812141.set"``.
    load_position
        If ``True`` read and post-process the ``.pos`` file.
    load_lfp
        If ``True``, load every EEG channel flagged as saved in the ``.set``
        file into ``Recording.lfp`` (preferring the high-rate ``.egf`` over
        ``.eeg`` where both exist). If a list of channel numbers, load only
        those. ``False`` (default) skips LFP loading.
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
            led_pos,
            led_pix,
            pos_header,
            settings,
            max_speed=max_speed,
            box_car=box_car,
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
        tet_path = f"{stem}.{t}"
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

        cut_path = f"{stem}_{t}.cut"
        clusters = _read_cut(cut_path) if os.path.exists(cut_path) else None

        rec.electrodes[t] = Electrode(
            identifier=t,
            spike_times=spike_times,
            pos_index=pos_index,
            clusters=clusters,
            metadata={"header": tet_header},
        )

    # --- LFP (.eeg / .egf) ----------------------------------------------
    if load_lfp:
        channels = (
            load_lfp if isinstance(load_lfp, (list, tuple)) else _active_eeg_channels(settings)
        )
        for ch in channels:
            channel_data = _load_lfp_channel(stem, ch, settings)
            if channel_data is not None:
                rec.lfp[ch] = channel_data

    return rec
