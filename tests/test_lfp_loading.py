"""Tests for EEG/EGF (LFP) loading."""

import numpy as np


def test_load_lfp_all_channels(example_trial):
    import barrylab as bl

    rec = bl.load_dacq(example_trial, load_lfp=True)
    assert len(rec.lfp) >= 1
    for ch, d in rec.lfp.items():
        assert d["channel"] == ch
        assert d["kind"] in ("eeg", "egf")
        assert d["samples"].ndim == 1
        assert np.all(np.isfinite(d["samples"]))
        # duration implied by n_samples/sample_rate should match the session
        implied_duration = d["samples"].size / d["sample_rate"]
        assert abs(implied_duration - rec.duration) < 1.0


def test_load_lfp_selective_channels(example_trial):
    import barrylab as bl

    rec = bl.load_dacq(example_trial, load_lfp=[1])
    assert list(rec.lfp.keys()) == [1]


def test_load_lfp_false_loads_nothing(example_trial):
    import barrylab as bl

    rec = bl.load_dacq(example_trial, load_lfp=False)
    assert rec.lfp == {}


def test_egf_channel_prefers_high_rate_over_eeg(example_trial):
    import barrylab as bl

    rec = bl.load_dacq(example_trial, load_lfp=[1])
    # the example trial ships both .eeg (250Hz) and .egf (4.8kHz) for ch1;
    # the loader should prefer the higher-rate .egf.
    assert rec.lfp[1]["kind"] == "egf"
    assert rec.lfp[1]["sample_rate"] > 1000


def test_lfp_shows_plausible_theta_peak(example_trial):
    """A crude physiological sanity check: the dominant 4-12Hz frequency in a
    foraging rat's hippocampal LFP should fall in the classic theta band."""
    import barrylab as bl

    rec = bl.load_dacq(example_trial, load_lfp=[1])
    sig = rec.lfp[1]["samples"]
    fs = rec.lfp[1]["sample_rate"]
    n = min(len(sig), int(fs * 60))
    spec = np.abs(np.fft.rfft(sig[:n] - sig[:n].mean())) ** 2
    freqs = np.fft.rfftfreq(n, d=1.0 / fs)
    theta_mask = (freqs >= 4) & (freqs <= 12)
    peak_freq = freqs[theta_mask][np.argmax(spec[theta_mask])]
    assert 5.0 <= peak_freq <= 11.0


def test_eeg_convert_to_volts_scales_by_adc_and_gain():
    from barrylab.io.axona import eeg_convert_to_volts

    settings = [
        ("ADC_fullscale_mv", "1500"),
        ("gain_ch_0", "1000"),
    ]
    raw = np.array([64.0, -64.0, 0.0])  # half-scale int8 counts
    volts = eeg_convert_to_volts(raw, eeg_channel=1, settings=settings, kind="eeg")
    expected = (1.5 / 1000.0) * (raw / 128.0)
    assert np.allclose(volts, expected)


def test_eeg_convert_to_volts_rejects_unknown_kind():
    import pytest

    from barrylab.io.axona import eeg_convert_to_volts

    settings = [("ADC_fullscale_mv", "1500"), ("gain_ch_0", "1000")]
    with pytest.raises(ValueError):
        eeg_convert_to_volts(np.array([1.0]), 1, settings, kind="bogus")


def test_load_lfp_gain_uses_actual_hardware_channel_not_slot(example_trial):
    """Regression for the slot-vs-actual-channel gain bug: the per-channel gain
    must be keyed by the hardware channel (value of EEG_ch_{slot}), not the EEG
    slot number. Exercises the full load path (the direct eeg_convert_to_volts
    unit test bypasses it). On r1838, EEG slot 1 -> hardware channel 4 ->
    gain_ch_3, NOT gain_ch_0."""
    import barrylab as bl
    from barrylab.io._axona_files import key_value, read_axona_file
    from barrylab.io.axona import _read_egf

    rec = bl.load_dacq(example_trial, load_lfp=[1])
    settings, _ = read_axona_file(example_trial + ".set", has_binary=False)
    raw, _ = _read_egf(example_trial + ".egf")

    actual = int(key_value(settings, "EEG_ch_1", "num", exact=True))
    assert rec.lfp[1]["hardware_channel"] == actual  # slot 1 resolved to its hw channel
    gain = key_value(settings, f"gain_ch_{actual - 1}", "num", exact=True)
    expected = (1.5 / gain) * (raw / 32768.0)
    np.testing.assert_allclose(rec.lfp[1]["samples"], expected)
