"""Tests for read_inp / read_stm / write_cut.

No example ``.inp``/``.stm`` file ships with ``Universal_Matlab``'s
``Example_data`` (confirmed by directly checking that directory), so these two
readers are validated structurally: build a byte-exact synthetic file from the
documented header/record layout and check the parsed values against known
ground truth. ``write_cut`` is validated by a genuine round trip against the
already-validated ``_read_cut``.
"""

import os
import struct

import numpy as np
import pytest

from barrylab.io.axona import _read_cut, read_inp, read_stm, write_cut


def _write_synthetic_axona_file(path, header_lines, binary):
    """Build a minimal Axona-style file: text header + data_start + binary + data_end."""
    header_text = "".join(line + "\r\n" for line in header_lines)
    with open(path, "wb") as fh:
        fh.write(header_text.encode("ascii"))
        fh.write(b"data_start")
        fh.write(binary)
        fh.write(b"\r\ndata_end\r\n")


# --------------------------------------------------------------------------
# read_inp
# --------------------------------------------------------------------------
def _inp_event(t_samples, type_char, byte6, byte7):
    """Pack one 7-byte .inp record: 4-byte big-endian time, 1-byte type, 2 chan bytes."""
    return struct.pack(">i", t_samples) + type_char.encode("ascii") + bytes([byte6, byte7])


def test_read_inp_parses_digital_onset_and_offset(tmp_path):
    timebase = 1000  # Hz
    # Event 1: digital input onset on channel 12 (byte6 bit for ch9-16: ch12 is
    # bit index 12-9=3 -> value 2**3=8). Byte7 = 0 (no ch1-8 bit set).
    ev1 = _inp_event(500, "I", 8, 0)
    # Event 2: the matching "off" signal -- both channel bytes 0.
    ev2 = _inp_event(750, "I", 0, 0)
    # Event 3: digital input onset on channel 3 (byte7 bit index 3-1=2 -> 4).
    ev3 = _inp_event(1000, "I", 0, 4)
    # Event 4: a keypress -- channel bytes are meaningless/not decoded.
    ev4 = _inp_event(2000, "K", 5, 9)

    header = [
        "trial_date Wednesday, 15 Jan 2020",
        "trial_time 12:00:00",
        "num_INP_samples 4",
        f"timebase {timebase} hz",
        "bytes_per_timestamp 4",
        "bytes_per_type 1",
        "bytes_per_channel 2",
    ]
    path = os.path.join(str(tmp_path), "test.inp")
    _write_synthetic_axona_file(path, header, ev1 + ev2 + ev3 + ev4)

    time, event_type, chan = read_inp(path)

    assert time.shape == (4,)
    np.testing.assert_allclose(time, np.array([500, 750, 1000, 2000]) / timebase)
    assert list(event_type) == ["I", "I", "I", "K"]
    assert chan[0] == 12  # onset, channel 12
    assert chan[1] == 0  # off signal
    assert chan[2] == 3  # onset, channel 3
    assert np.isnan(chan[3])  # keypress channel not decoded


def test_read_inp_truncates_trailing_partial_record(tmp_path):
    """A trailing incomplete 7-byte record is silently dropped, like _read_pos."""
    ev1 = _inp_event(100, "O", 1, 0)  # channel 9 (bit0 of byte6 -> ch9)
    path = os.path.join(str(tmp_path), "test.inp")
    _write_synthetic_axona_file(path, ["timebase 1000 hz"], ev1 + b"\x01\x02\x03")

    time, event_type, chan = read_inp(path)
    assert time.shape == (1,)
    assert event_type[0] == "O"
    assert chan[0] == 9


# --------------------------------------------------------------------------
# read_stm
# --------------------------------------------------------------------------
def test_read_stm_parses_uint32_timestamps(tmp_path):
    timebase = 1000
    values = [0, 1500, 3000, 12345]
    binary = b"".join(struct.pack(">I", v) for v in values)
    header = [
        "ClientVersion 1.0",
        "trial_date Wednesday, 15 Jan 2020",
        "trial_time 12:00:00",
        "data_path C:\\data",
        "set_file_name test",
        "experimenter CB",
        "comments none",
        f"timebase {timebase} hz",
        "duration 100",
        "spike_format none",
        f"num_stm_samples {len(values)}",
    ]
    path = os.path.join(str(tmp_path), "test.stm")
    _write_synthetic_axona_file(path, header, binary)

    stim_time, header_out = read_stm(path)
    np.testing.assert_allclose(stim_time, np.array(values) / timebase)
    assert any(k == "timebase" for k, _ in header_out)


def test_read_stm_no_stimulator_samples(tmp_path):
    path = os.path.join(str(tmp_path), "empty.stm")
    _write_synthetic_axona_file(path, ["timebase 1000 hz"], b"")
    stim_time, _ = read_stm(path)
    assert stim_time.shape == (0,)


# --------------------------------------------------------------------------
# write_cut
# --------------------------------------------------------------------------
def test_write_cut_round_trips_small_assignment(tmp_path):
    """The exact case the original's stale-loop-variable bug would break:
    fewer than 25 spikes."""
    assignment = np.array([4, 1, 2, 7, 4, 2, 12, 4, 3])
    path = write_cut(assignment, str(tmp_path), "2140612a", 1)
    assert os.path.basename(path) == "2140612a_1.cut"

    recovered = _read_cut(path)
    np.testing.assert_array_equal(recovered, assignment)


def test_write_cut_round_trips_multi_row_assignment(tmp_path):
    """>25 spikes exercises the full-row / tail split."""
    rng = np.random.RandomState(0)
    assignment = rng.randint(0, 6, size=137)
    path = write_cut(assignment, str(tmp_path), "83812141", 2)

    recovered = _read_cut(path)
    np.testing.assert_array_equal(recovered, assignment)


def test_write_cut_round_trips_exact_multiple_of_25(tmp_path):
    assignment = np.arange(25) % 4
    path = write_cut(assignment, str(tmp_path), "stem", 3)
    recovered = _read_cut(path)
    np.testing.assert_array_equal(recovered, assignment)


def test_write_cut_embeds_actual_file_stem_not_hardcoded_literal(tmp_path):
    """The confirmed PORTING_PLAN.md bug: MATLAB hard-codes 'Exact_cut_for:
    1500728a' regardless of the fileStem argument. This port must not."""
    assignment = np.array([0, 1, 0, 1])
    path = write_cut(assignment, str(tmp_path), "a_totally_different_stem", 5)
    with open(path) as fh:
        text = fh.read()
    assert "Exact_cut_for: a_totally_different_stem" in text
    assert "1500728a" not in text


def test_write_cut_header_declares_n_clusters_and_reserves_cluster_0(tmp_path):
    assignment = np.array([0, 2, 2, 5, 0])
    path = write_cut(assignment, str(tmp_path), "stem", 1)
    with open(path) as fh:
        text = fh.read()
    assert "n_clusters: 6" in text  # clusters 0..5 declared (max id 5 + 1)
    assert " cluster:  0 centre:" in text
    assert " cluster:  5 centre:" in text


def test_write_cut_warns_above_max_cluster(tmp_path):
    assignment = np.array([0, 40])
    with pytest.warns(UserWarning):
        write_cut(assignment, str(tmp_path), "stem", 1)


def test_write_cut_rejects_empty_assignment(tmp_path):
    with pytest.raises(ValueError):
        write_cut(np.array([], dtype=int), str(tmp_path), "stem", 1)
