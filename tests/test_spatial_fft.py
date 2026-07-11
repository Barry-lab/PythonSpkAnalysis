"""Tests for the 2-D FFT grid-analysis subpackage (``2DFFT_Anl/``)."""

import numpy as np
import pytest

from barrylab.analysis.spatial_fft import (
    FFTParams,
    _fspecial_gaussian_square,
    _matlab_round,
    _smooth_polar_fft,
    grid_fft_components,
    power_spectrum_2d,
    to_polar_fft,
)
from barrylab.core.recording import Electrode, Recording


# --------------------------------------------------------------------------
# FFTParams (fftVars.m)
# --------------------------------------------------------------------------
def test_fft_params_max_rad_matches_matlab_formula_and_tracks_field_changes():
    p = FFTParams()
    assert p.max_rad == int(np.ceil((1.0 / p.min_wave_cm) * (p.bin_size_cm * p.fft_size)))
    assert p.max_rad == 86  # defaults: ceil((1/12)*(2*512)) = ceil(85.33) = 86

    # derived, not cached -- changing a dependency must be reflected immediately
    p.fft_size = 256
    assert p.max_rad == int(np.ceil((1.0 / p.min_wave_cm) * (p.bin_size_cm * p.fft_size)))
    assert p.max_rad == 43


# --------------------------------------------------------------------------
# _fspecial_gaussian_square helper
# --------------------------------------------------------------------------
def test_fspecial_gaussian_square_normalised_symmetric():
    kern = _fspecial_gaussian_square(9, 2.0)
    assert kern.shape == (9, 9)
    assert np.isclose(kern.sum(), 1.0)
    assert np.allclose(kern, kern[::-1, ::-1])  # point-symmetric about centre
    assert kern[4, 4] == kern.max()  # peak at the centre


def test_fspecial_gaussian_square_even_hsize_matches_matlab_half_integer_grid():
    # MATLAB fspecial('gaussian', 4, 1): siz=(4-1)/2=1.5, coords=[-1.5,-0.5,0.5,1.5]
    kern = _fspecial_gaussian_square(4, 1.0)
    assert kern.shape == (4, 4)
    assert np.isclose(kern.sum(), 1.0)
    coords = np.arange(4) - 1.5
    xx, yy = np.meshgrid(coords, coords)
    expected = np.exp(-(xx**2 + yy**2) / 2.0)
    expected = expected / expected.sum()
    assert np.allclose(kern, expected)


# --------------------------------------------------------------------------
# power_spectrum_2d (get_2d_fft.m)
# --------------------------------------------------------------------------
def test_power_spectrum_2d_removes_dc_component():
    rng = np.random.RandomState(0)
    rm = rng.uniform(1.0, 5.0, size=(20, 20))
    ps2d = power_spectrum_2d(rm, 64, 0, 0)
    cent = 64 // 2
    assert abs(ps2d[cent, cent]) < 1e-8  # 0th (DC) component removed by mean-normalisation


def test_power_spectrum_2d_recovers_known_spatial_frequency():
    # A pure grating along x with an exact integer number of cycles across the
    # (unpadded) ratemap gives a clean, unambiguous single-frequency peak.
    n, wavelength = 64, 8
    grating = 1.0 + np.cos(2 * np.pi * np.arange(n) / wavelength)
    rm = np.tile(grating, (n, 1))  # varies along columns (x) only -> angle 0/180
    ps2d = power_spectrum_2d(rm, n, 0, 0)
    cent = n // 2
    ps2d_no_dc = ps2d.copy()
    ps2d_no_dc[cent, cent] = -1.0  # exclude the (already ~0) DC bin from argmax
    row, col = np.unravel_index(np.argmax(ps2d_no_dc), ps2d.shape)
    assert row == cent  # all energy on the centre row -> pure x-direction grating
    assert abs(col - cent) == n // wavelength  # peak at +/- fft_size/wavelength bins


def test_power_spectrum_2d_max_rad_masks_outside_radius_only():
    n = 40
    rm = np.random.RandomState(1).uniform(size=(15, 15))
    max_rad = 10
    ps2d = power_spectrum_2d(rm, n, 0, max_rad)
    cent = n // 2
    yy, xx = np.mgrid[0:n, 0:n]
    dist = np.hypot(xx - cent, yy - cent)
    assert np.all(np.isnan(ps2d[dist > max_rad]))
    assert np.all(np.isfinite(ps2d[dist <= max_rad]))


def test_power_spectrum_2d_smoothing_preserves_peak_location_and_stays_finite():
    n, wavelength = 64, 8
    grating = 1.0 + np.cos(2 * np.pi * np.arange(n) / wavelength)
    rm = np.tile(grating, (n, 1))
    ps2d_smooth = power_spectrum_2d(rm, n, 1.5, 0)
    assert np.all(np.isfinite(ps2d_smooth))  # edge-corrected -- no NaN/inf from the division
    cent = n // 2
    ps2d_smooth_no_dc = ps2d_smooth.copy()
    ps2d_smooth_no_dc[cent, cent] = -1.0
    row, col = np.unravel_index(np.argmax(ps2d_smooth_no_dc), ps2d_smooth.shape)
    # ceil(smth_sig*4) with smth_sig=1.5 gives an EVEN kernel size (6, not
    # forced odd -- unlike _smooth_polar_fft's kernel, get_2d_fft.m's own
    # fspecial('gaussian', ceil(smthSig.*4), smthSig) call has no such
    # odd-forcing step either). An even-sized kernel has no single centre
    # tap (its coordinates straddle 0, e.g. [-2.5,...,2.5]), so correlating
    # it with a peak sitting exactly on an integer bin genuinely, correctly
    # ties the two straddling output bins -- confirmed directly (both rows
    # 31 and 32 come back bit-identical here) -- so allow either.
    assert abs(row - cent) <= 1
    assert abs(col - (cent + n // wavelength)) <= 1 or abs(col - (cent - n // wavelength)) <= 1


# --------------------------------------------------------------------------
# to_polar_fft (get_polar_fft.m)
# --------------------------------------------------------------------------
def test_to_polar_fft_places_single_peak_at_expected_angle_and_radius():
    n = 41
    cent = n // 2
    ps2d = np.zeros((n, n))
    dx, dy = 6, -4  # a peak above-and-right of centre, no ties
    row, col = cent + dy, cent + dx
    ps2d[row, col] = 10.0

    ps_pol, th, frq = to_polar_fft(ps2d, max_rad=15, th_bin=1.0, r_bin=1.0, sigma=None)

    expected_r = int(round(np.hypot(dx, dy)))
    expected_th = int(round(np.degrees(np.arctan2(-dy, dx)))) % 180
    assert ps_pol[expected_th, expected_r] == 10.0
    assert np.nanargmax(ps_pol) == np.ravel_multi_index((expected_th, expected_r), ps_pol.shape)
    assert th[expected_th] == expected_th
    assert frq[expected_r] == expected_r


def test_to_polar_fft_rejects_non_square_input():
    with pytest.raises(ValueError):
        to_polar_fft(np.zeros((10, 12)), max_rad=3)


def test_to_polar_fft_rejects_max_rad_larger_than_half_size():
    with pytest.raises(ValueError):
        to_polar_fft(np.zeros((21, 21)), max_rad=11)  # cent=10, max_rad=11 > cent


def test_to_polar_fft_scalar_sigma_matches_equal_tuple_sigma():
    n = 31
    ps2d = np.random.RandomState(2).uniform(size=(n, n))
    a, _, _ = to_polar_fft(ps2d, max_rad=10, sigma=3.0)
    b, _, _ = to_polar_fft(ps2d, max_rad=10, sigma=(3.0, 3.0))
    assert np.allclose(a, b, equal_nan=True)


def test_to_polar_fft_binning_matches_brute_force_reference_including_nan_propagation():
    """The decisive test for get_polar_fft.m's rebinning step: an independent,
    literal re-implementation of the (round, mod, mean-with-explicit-nan-fill)
    accumarray logic, cross-checked bin-for-bin against ``to_polar_fft``.

    Deliberately uses coarse bins (so many pixels collide into the same
    (angle, radius) cell -- exercising the mean aggregation for real, not just
    singleton bins) and injects NaN pixels (simulating power_spectrum_2d's own
    max_rad masking) to confirm MATLAB's plain ``mean()`` NaN-propagation
    semantics (NOT ``nanmean``) -- a genuinely different convention from the
    default-0-fill radial-profile aggregation in grid.py::fast_gridness, see
    FINDINGS.md.
    """
    rng = np.random.RandomState(3)
    n = 31
    cent = n // 2
    ps2d = rng.uniform(0.0, 5.0, size=(n, n))
    ps2d[rng.rand(n, n) < 0.1] = np.nan

    max_rad = 12
    th_bin, r_bin = 7.0, 3.0  # coarse -> guarantees genuine multi-pixel bins
    ps_pol, th, frq = to_polar_fft(ps2d, max_rad, th_bin=th_bin, r_bin=r_bin, sigma=None)

    n_th_bins = int(round(180.0 / th_bin))
    x = np.arange(-max_rad, max_rad + 1)
    y = np.arange(-max_rad, 1)
    sums, counts = {}, {}
    for dy in y:
        for dx in x:
            val = ps2d[cent + dy, cent + dx]
            th_ang = abs(np.degrees(np.arctan2(-dy, dx)))
            r = np.hypot(dx, dy)
            r_sub = int(_matlab_round(np.array(r / r_bin)))
            th_sub = int(_matlab_round(np.array(th_ang / th_bin))) % n_th_bins
            key = (th_sub, r_sub)
            sums[key] = sums.get(key, 0.0) + val
            counts[key] = counts.get(key, 0) + 1

    max_r_sub = max(k[1] for k in counts)
    ref = np.full((n_th_bins, max_r_sub + 1), np.nan)
    for key, c in counts.items():
        ref[key] = sums[key] / c
    ref = ref[:, : max_rad + 1]

    assert ps_pol.shape == ref.shape
    assert np.array_equal(np.isnan(ps_pol), np.isnan(ref))  # exact NaN-mask match
    finite = ~np.isnan(ref)
    assert np.allclose(ps_pol[finite], ref[finite], atol=1e-9)


def test_to_polar_fft_empty_bins_are_nan_not_zero():
    # accumarray(..., @(x) mean(x), nan) -- an explicit NaN fill, unlike the
    # default-0 fill grid.py::fast_gridness's own (different) accumarray-style
    # aggregation relies on. With max_rad well beyond what a (th_bin, r_bin)
    # grid this coarse can ever be fully populated by, some bins MUST be
    # genuinely untouched (count==0) -- distinct from a bin whose mean
    # legitimately computes to 0 from real (zero-valued) contributing pixels.
    # (dy=-5, dx=3): deliberately off the y=0 row, which would otherwise
    # collide with its own mod-180 mirror at (dy=0, -dx) -- see the
    # dedicated wraparound-merging test below -- and dilute this bin's mean.
    n = 81
    ps2d = np.zeros((n, n))
    ps2d[n // 2 - 5, n // 2 + 3] = 1.0
    ps_pol, th, frq = to_polar_fft(ps2d, max_rad=35, th_bin=1.0, r_bin=1.0, sigma=None)
    hit = ~np.isnan(ps_pol)
    assert hit.sum() > 0
    assert (~hit).sum() > 0  # some bins are genuinely never touched -> NaN, not 0
    assert ps_pol[hit].max() == 1.0  # the one real peak survived the (mean, NaN-fill) aggregation


def test_to_polar_fft_wraps_0_and_180_degrees_into_the_same_bin():
    # th=0 and th=180 represent the same undirected axis of a point-symmetric
    # power spectrum -- get_polar_fft.m deliberately merges them
    # (`mod(round(th/thBin), 180/thBin)+1`). A point on the y==0 row at +dx
    # and its mirror at -dx therefore land in the SAME (angle, radius) bin
    # and get averaged together, not kept separate.
    n = 61
    cent = n // 2
    ps2d = np.zeros((n, n))
    ps2d[cent, cent + 5] = 1.0
    ps2d[cent, cent - 5] = 3.0
    ps_pol, th, frq = to_polar_fft(ps2d, max_rad=15, th_bin=1.0, r_bin=1.0, sigma=None)
    assert ps_pol[0, 5] == 2.0  # mean of the two merged mirror points


# --------------------------------------------------------------------------
# _smooth_polar_fft (the smoothing block inside get_polar_fft.m)
# --------------------------------------------------------------------------
def test_smooth_polar_fft_edge_corrected_constant_field_invariant():
    const = np.full((180, 20), 3.0)
    sm = _smooth_polar_fft(const, r_sig=2.0, th_sig=5.0)
    assert np.allclose(sm, 3.0, atol=1e-8)


def test_smooth_polar_fft_wraps_the_angular_axis():
    ps_pol = np.zeros((180, 20))
    ps_pol[0, 10] = 100.0
    sm = _smooth_polar_fft(ps_pol, r_sig=2.0, th_sig=10.0)
    assert sm[179, 10] > 0.0  # row 179 is the angular neighbour of row 0 (circular)
    assert sm[90, 10] == 0.0  # far from the impulse on the circle -> untouched


def test_smooth_polar_fft_does_not_wrap_the_radial_axis():
    ps_pol = np.zeros((180, 20))
    ps_pol[90, 19] = 100.0  # impulse at the far (max) radius edge
    sm = _smooth_polar_fft(ps_pol, r_sig=3.0, th_sig=2.0)
    assert sm[90, 18] > 0.0  # spreads towards lower radius...
    assert sm[90, 0] == 0.0  # ...but never wraps all the way to r=0


# --------------------------------------------------------------------------
# grid_fft_components (main_spatial_fft2.m) -- Recording-level orchestrator.
# Redesigned around Recording like decode_crunch_data/field_shuffle_significance;
# not bit-exact. Validated structurally + a synthetic known-orientation/
# known-wavelength grating (the 2-D analogue of test_grid_phase.py's
# hand-built known-peak fixture), plus a real example grid cell.
# --------------------------------------------------------------------------
def _synthetic_grating_recording(wavelength_cm=24.0, amplitude=5.0, theta_deg=60.0, seed=0):
    """A Recording whose unsmoothed ratemap is (up to integer spike-count
    rounding) an exact sinusoidal grating of known wavelength/orientation.

    Positions sit exactly at bin centres (R repeats per bin, in a RANDOM
    time-order -- not a raster scan -- so the temporal shuffle's circular
    rotation genuinely decorrelates spike-to-position assignment; an ordered
    raster trajectory would make even a "random" rotation land back on a
    highly self-similar position sequence, defeating the shuffle-null
    assumption). Each bin gets an identical, deterministic spike-count
    (proportional to the target rate shape) attached to every one of its
    visits, so the resulting ratemap tracks the analytic grating closely
    (no Poisson sampling noise) -- a deterministic-style fixture in the
    spirit of test_grid_phase.py's hand-built known-peak SAC, adapted to a
    full Recording pipeline.

    NOTE on angle convention: ``theta_deg`` is the grating's orientation in
    the ordinary (y-up) mathematical sense. Recording.xy's y increases
    downward (image convention, like the ratemap array itself), so the
    Fourier component this produces appears at angle ``(-theta_deg) % 180``
    under this package's cart2pol(dx, -dy) convention (see
    get_polar_fft.m's own "y increases down... need to invert" comment,
    already established in grid_phase.py/field_shuffle.py) -- confirmed
    empirically below, not just asserted.
    """
    n_side, bin_cm, ppm, r_visits, pos_sample_rate = 24, 4.0, 100.0, 20, 50.0
    bin_px = ppm / 100.0 * bin_cm
    xs_bin = (np.arange(n_side) + 0.5) * bin_px
    ys_bin = (np.arange(n_side) + 0.5) * bin_px
    gx, gy = np.meshgrid(xs_bin, ys_bin)
    bin_x, bin_y = gx.ravel(), gy.ravel()
    n_bins_total = n_side * n_side

    if amplitude:
        th_r = np.radians(theta_deg)
        proj = bin_x * np.cos(th_r) + bin_y * np.sin(th_r)
        rate_shape = 5.0 + amplitude * np.cos(2 * np.pi * proj / wavelength_cm)
    else:
        rate_shape = np.full(n_bins_total, 5.0)
    k_per_bin = np.clip(np.round(rate_shape / 10.0 * 16).astype(int), 0, None)

    bin_id_per_sample = np.tile(np.arange(n_bins_total), r_visits)
    rng = np.random.RandomState(seed)
    bin_id_per_sample = bin_id_per_sample[rng.permutation(bin_id_per_sample.size)]
    pos_xy = np.column_stack([bin_x, bin_y])[bin_id_per_sample]
    n_pos = pos_xy.shape[0]
    spk_idx = np.repeat(np.arange(n_pos), k_per_bin[bin_id_per_sample])

    pos_times = np.arange(n_pos) / pos_sample_rate
    window = (0.0, n_side * bin_px, 0.0, n_side * bin_px)
    rec = Recording(
        name="synthetic_grating",
        xy=pos_xy,
        direction=np.zeros(n_pos),
        speed=np.zeros(n_pos),
        pos_times=pos_times,
        pos_sample_rate=pos_sample_rate,
        pixels_per_metre=ppm,
        window=window,
        electrodes={
            1: Electrode(
                1,
                spike_times=pos_times[spk_idx],
                pos_index=spk_idx,
                clusters=np.zeros(len(spk_idx), dtype=int),
            )
        },
    )
    return rec, bin_cm


def test_grid_fft_components_recovers_known_grating_angle_and_wavelength():
    rec, bin_cm = _synthetic_grating_recording(
        wavelength_cm=24.0, amplitude=5.0, theta_deg=60.0, seed=0
    )
    params = FFTParams()
    params.bin_size_cm = bin_cm
    params.fft_size = 64  # small for test speed; n_shuffles small per porting-brief convention

    result = grid_fft_components(rec, 1, 0, params=params, n_shuffles=8, rng=0)

    assert result["ps2d"].shape == (64, 64)
    assert result["ps_polar"].shape[1] == params.max_rad + 1
    assert result["shuffle_max_power"].shape == (8,)
    assert result["significant"]
    assert result["peak_angle_deg"].size > 0

    expected_angle = (-60.0) % 180  # see _synthetic_grating_recording's angle-convention note
    angle_err = np.min(np.abs(result["peak_angle_deg"] - expected_angle))
    assert angle_err <= 10.0

    wavelength_err = np.min(np.abs(result["peak_wavelength_cm"] - 24.0))
    assert wavelength_err <= 10.0  # zero-padding leakage/discretisation broadens this a bit


def test_grid_fft_components_flat_ratemap_is_not_significant():
    # A no-modulation ("flat") control: with no genuine spatial structure the
    # real polar FFT should not exceed its own temporal-shuffle 95th
    # percentile. Empirically stable across many seeds (checked 0-15 while
    # developing this test); a low-probability false positive is possible in
    # principle (this is a genuine significance test) but not expected here.
    rec, bin_cm = _synthetic_grating_recording(amplitude=0.0, seed=11)
    params = FFTParams()
    params.bin_size_cm = bin_cm
    params.fft_size = 64

    result = grid_fft_components(rec, 1, 0, params=params, n_shuffles=8, rng=11)
    assert not result["significant"]
    assert result["peak_angle_deg"].size == 0
    assert result["peak_wavelength_cm"].size == 0


def test_grid_fft_components_runs_on_real_example_grid_cell(example_trial):
    import barrylab as bl

    rec = bl.load_dacq(example_trial)
    params = FFTParams()  # default fft_size=512 -- fine at n_shuffles=5 (fast, see below)
    result = grid_fft_components(rec, 2, 3, params=params, n_shuffles=5, rng=0)

    assert result["ps2d"].shape == (params.fft_size, params.fft_size)
    assert result["ps_polar"].shape[1] == params.max_rad + 1
    assert result["shuffle_max_power"].shape == (5,)
    assert np.isfinite(result["crit95"])
    if result["significant"]:
        # documented grid scale for t2c3 is ~40cm (FINDINGS.md); a Fourier
        # component's wavelength should be shorter than the grid scale
        # (cos(30deg)*scale for a perfectly regular grid) and not absurd.
        assert np.all(result["peak_wavelength_cm"] > 0)
        assert np.any(result["peak_wavelength_cm"] < 40.0)
