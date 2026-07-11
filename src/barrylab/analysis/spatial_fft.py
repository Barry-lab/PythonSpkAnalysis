"""2-D FFT grid-cell analysis.

Ports ``Analysis_Grid/2DFFT_Anl/`` -- Julia Krupic-style spatial-frequency
analysis of grid-cell ratemaps ("bandy cell" detection), rewritten by the
lab "for understanding, to test rigour and for speed":

- :class:`FFTParams` -- the common parameter struct (``fftVars.m``).
- :func:`power_spectrum_2d` -- mean-normalised, zero-padded 2-D FFT power
  spectrum of a ratemap, with optional Gaussian smoothing and a radius mask
  (``get_2d_fft.m``).
- :func:`to_polar_fft` -- reorganise a 2-D FFT power spectrum into polar
  (angle, radius) coordinates, mean-aggregated per bin, with optional
  circular-in-angle Gaussian smoothing (``get_polar_fft.m``).
- :func:`grid_fft_components` -- the ``Recording``-level orchestrator (a
  redesign of ``main_spatial_fft2.m``, like :func:`~barrylab.analysis.decoding.decode_crunch_data`
  and :func:`~barrylab.analysis.field_shuffle.field_shuffle_significance`):
  builds an unsmoothed ratemap, gets its 2-D and polar FFT, and tests each
  angular Fourier component for significance against a temporal-shuffle null
  distribution.

``main_spatial_fft.m`` (the v1 predecessor) is deliberately **not** ported --
its own header says to use v2, it bins with the superseded ``bin_pos_data``,
and it starts parallel workers with the obsolete ``goPar`` -- see
``PORTING_PLAN.md``.
"""

import warnings
from dataclasses import dataclass

import numpy as np
from scipy.ndimage import convolve
from scipy.signal import correlate2d, find_peaks

from .maps import occupancy_map, smooth_ratemap, spike_map


def _matlab_round(x):
    """MATLAB ``round()``: half away from zero (numpy's default is half to even)."""
    return np.sign(x) * np.floor(np.abs(x) + 0.5)


def _fspecial_gaussian_square(hsize, sigma):
    """MATLAB ``fspecial('gaussian', hsize, sigma)``: a normalised ``hsize x
    hsize`` kernel. ``hsize`` need not be odd -- MATLAB's ``-siz:siz`` colon
    range (``siz=(hsize-1)/2``, possibly a half-integer) always produces
    exactly ``hsize`` samples spaced 1 apart, symmetric about 0 (not
    necessarily centred on an integer coordinate when ``hsize`` is even).
    """
    siz = (hsize - 1) / 2.0
    coords = np.arange(hsize, dtype=float) - siz
    xx, yy = np.meshgrid(coords, coords)
    h = np.exp(-(xx**2 + yy**2) / (2.0 * sigma**2))
    h[h < np.finfo(float).eps * h.max()] = 0.0
    total = h.sum()
    if total != 0:
        h = h / total
    return h


def _as_random_state(rng):
    if isinstance(rng, np.random.RandomState):
        return rng
    return np.random.RandomState(rng)


# --------------------------------------------------------------------------
# fftVars.m -- the common parameter struct
# --------------------------------------------------------------------------
@dataclass
class FFTParams:
    """Parameters shared by the 2-D FFT grid-analysis family. Port of ``fftVars.m``.

    Kept as a plain dataclass of defaults (matching ``fftVars.m``'s role as a
    parameter struct, not an algorithm) rather than folded into
    ``params.py``'s shared constants -- this whole family of parameters is
    self-contained to this one analysis, unlike e.g. ``RATEMAP_BIN_CM`` which
    many unrelated functions share.

    ``max_rad`` is a derived property (not a stored field), exactly like
    MATLAB's ``f.maxRad=ceil((1/f.minWaveCm)*(f.binSizeCm*f.fftSize))`` --
    recomputed from ``min_wave_cm``/``bin_size_cm``/``fft_size`` so it can
    never go stale if any of those three are changed after construction.
    """

    #: Spatial bin size in cm (converted to pixels via the recording's ppm).
    bin_size_cm = 2.0
    #: Minimum temporal-shuffle rotation magnitude, in seconds.
    min_shuf_s = 20.0
    #: FFT is zero-padded to be square with sides of this length (ideally a
    #: power of 2). JK used 256; the MATLAB author's own comment recommends
    #: 512-1024 for accuracy in the 1-4 cm bin-size range.
    fft_size = 512
    #: Sigma (FFT bins) to smooth the 2-D power spectrum with (0 = no smooth;
    #: the MATLAB author's own comment recommends NOT smoothing in 2-D coords).
    smth_2d_sig = 0.0
    #: Minimum wavelength (cm) that should survive the radius mask; converted
    #: to ``max_rad`` (see the class docstring).
    min_wave_cm = 12.0
    #: Angular bin size in degrees for the polar FFT.
    th_bin = 1.0
    #: Radial bin size (FFT components) for the polar FFT -- best left at 1
    #: to match the 2-D FFT's own bin spacing (see :func:`to_polar_fft`).
    r_bin = 1.0
    #: Sigma (degrees) for angular (circular) Gaussian smoothing of the polar FFT.
    th_sig = 10.0
    #: Sigma (FFT components) for radial Gaussian smoothing of the polar FFT.
    r_sig = 4.0
    #: Number of temporal-shuffle iterations for the significance threshold.
    n_it = 1000
    #: Percentile of the shuffled max-power distribution used as the
    #: significance threshold for the real polar FFT.
    p_thresh = 50.0
    #: Secondary angular-peak threshold, as a fraction of the tallest peak.
    snd_peak_coef = 0.2
    #: Minimum separation between angular peaks, in degrees.
    min_peak_dist = 10.0

    @property
    def max_rad(self):
        """Radius (FFT bins) beyond which spatial frequencies are discarded."""
        return int(np.ceil((1.0 / self.min_wave_cm) * (self.bin_size_cm * self.fft_size)))


# --------------------------------------------------------------------------
# get_2d_fft.m
# --------------------------------------------------------------------------
def power_spectrum_2d(rate_map_2d, fft_size, smth_sig=0.0, max_rad=0):
    """Mean-normalised, zero-padded 2-D FFT power spectrum of a ratemap.

    Port of ``get_2d_fft.m``. Ideally takes an *unsmoothed* ratemap (see the
    MATLAB docstring), mean-normalises it (removing the 0th/DC component),
    zero-pads to a square ``fft_size x fft_size`` FFT, and returns the
    fftshifted magnitude spectrum (low frequency at the centre).

    Parameters
    ----------
    rate_map_2d : ndarray
        A ratemap; NaN (unvisited) bins are treated as 0 after mean-removal.
    fft_size : int
        Zero-pad to this square size (ideally a power of 2).
    smth_sig : float
        Sigma (FFT bins) of an edge-corrected Gaussian smooth applied to the
        power spectrum (0 = no smoothing; not recommended by the MATLAB
        author's own comment -- best not to smooth in 2-D FFT coordinates).
    max_rad : int
        Set bins further than this many FFT-bins from the centre to NaN
        (0 = return every frequency, unmasked).

    Returns
    -------
    ndarray (fft_size, fft_size)
        The (optionally smoothed/masked) power spectrum.
    """
    rm = np.asarray(rate_map_2d, dtype=float)
    rm = rm - np.nanmean(rm)  # mean-normalise -- removes the 0th FFT component
    rm = np.where(np.isnan(rm), 0.0, rm)  # replace remaining NaN with 0
    ps2d = np.abs(np.fft.fft2(rm, s=(int(fft_size), int(fft_size))))  # FFT w/ zero-padding
    ps2d = np.fft.fftshift(ps2d)  # low freq to centre

    if smth_sig:
        kern = _fspecial_gaussian_square(int(np.ceil(smth_sig * 4)), smth_sig)
        # Edge-corrected smoothing: divide by the smoothed all-ones mask so
        # edge bins (which see fewer kernel taps) aren't artificially dimmed
        # -- the identical trick used in maps.py::smooth_ratemap. filter2's
        # default MATLAB mode ('same', correlation not convolution) is
        # replicated with correlate2d(..., mode="same"); the Gaussian kernel
        # is symmetric so correlation and convolution coincide here anyway.
        ones_mask = np.ones_like(ps2d)
        denom = correlate2d(ones_mask, kern, mode="same")
        with np.errstate(invalid="ignore", divide="ignore"):
            ps2d = correlate2d(ps2d, kern, mode="same") / denom

    if max_rad:
        # Even-N centre after fftshift is at 0-based index N//2 -- matches
        # MATLAB's 1-based `floor(fftSize/2)+1` exactly (one less, 0-based).
        centre = int(fft_size) // 2
        yy, xx = np.mgrid[0 : int(fft_size), 0 : int(fft_size)]
        dist = np.hypot(xx - centre, yy - centre)
        ps2d = np.where(dist > max_rad, np.nan, ps2d)

    return ps2d


# --------------------------------------------------------------------------
# get_polar_fft.m
# --------------------------------------------------------------------------
def _smooth_polar_fft(ps_pol, r_sig, th_sig):
    """Edge-corrected 2-D Gaussian smoothing of a polar FFT: circular
    (wraps) in the angular dimension, non-wrapping in the radial dimension.

    Port of the smoothing block inside ``get_polar_fft.m``. The angular axis
    genuinely is periodic here (rows 0 and 179 are adjacent directions, both
    representing the same undirected axis of a point-symmetric power
    spectrum), but the radial axis is not (r=0 has no meaningful relationship
    to r=max). ``imfilter(..., 'circular')`` wraps *both* axes uniformly, so
    the MATLAB source works around this by padding extra NaN columns before
    the (never-wrapping-in-practice) low-radius edge, wide enough that the
    kernel can never reach past the padding to see genuine data at the other
    end -- this reproduces "no wrap in radius" using only a uniform-mode
    circular filter, without needing a mode that wraps one axis and not the
    other (which plain ``scipy.ndimage.convolve`` cannot express per-axis in
    the scipy version this project targets -- confirmed directly, a per-axis
    ``mode`` tuple raises "boundary mode not supported").

    Parameters
    ----------
    ps_pol : ndarray (n_th, n_r)
        Polar FFT (angle x radius); NaN marks unobserved bins.
    r_sig, th_sig : float
        Gaussian sigma, in radial FFT components / angular bins respectively.

    Returns
    -------
    ndarray (n_th, n_r), same shape as ``ps_pol``.
    """
    # Kernel is 4*sigma wide in each dimension, forced odd (MATLAB:
    # `kernSize=kernSize+mod(kernSize+1,2)` -- rounds an even size up by 1).
    kern_size_th = int(np.ceil(th_sig * 4))
    kern_size_r = int(np.ceil(r_sig * 4))
    if kern_size_th % 2 == 0:
        kern_size_th += 1
    if kern_size_r % 2 == 0:
        kern_size_r += 1

    r_off = np.arange(kern_size_r) - kern_size_r // 2
    th_off = np.arange(kern_size_th) - kern_size_th // 2
    xx, yy = np.meshgrid(r_off, th_off)  # xx: radius offset (cols); yy: angle offset (rows)
    kern = np.exp(-((xx**2) / (2.0 * r_sig**2) + (yy**2) / (2.0 * th_sig**2)))

    # Pad the radius axis (columns) with NaN on the low-r ("pre") side only --
    # wide enough (kern_size_r + 1, one more than the kernel's own half-width
    # margin) that the kernel can never wrap from one real edge to the other.
    pad_width = kern_size_r + 1
    padded = np.pad(
        np.asarray(ps_pol, dtype=float),
        ((0, 0), (pad_width, 0)),
        mode="constant",
        constant_values=np.nan,
    )
    norm = (~np.isnan(padded)).astype(float)
    padded = np.where(np.isnan(padded), 0.0, padded)

    smoothed = convolve(padded, kern, mode="wrap")
    norm_smoothed = convolve(norm, kern, mode="wrap")
    with np.errstate(invalid="ignore", divide="ignore"):
        smoothed = smoothed / norm_smoothed
    return smoothed[:, pad_width:]


def to_polar_fft(ps2d, max_rad, th_bin=1.0, r_bin=1.0, sigma=None):
    """Reorganise a 2-D FFT power spectrum into polar (angle, radius) coords.

    Port of ``get_polar_fft.m``. Because a real ratemap's 2-D FFT magnitude
    spectrum is point-symmetric about its centre, only the "top half" (angles
    0-180 deg) is computed, then rebinned by (rounded) angle and radius with a
    *mean* aggregation, and optionally smoothed.

    Parameters
    ----------
    ps2d : ndarray (N, N)
        A square 2-D FFT power spectrum, e.g. from :func:`power_spectrum_2d`
        (with the SAME ``max_rad`` -- see the note on NaN propagation below).
    max_rad : int
        Ignore/discard frequencies further than this many FFT bins from the
        centre (0 = use the whole array -- see the warning below).
    th_bin : float
        Angular bin size in degrees.
    r_bin : float
        Radial bin size in FFT components (best left at 1 to match ``ps2d``'s
        own bin spacing -- see the note on ``frq`` below).
    sigma : float, (r_sig, th_sig), or None
        Gaussian smoothing sigma for :func:`_smooth_polar_fft` (0/None skips
        smoothing). A scalar applies the same sigma to both dimensions,
        matching MATLAB's ``[rSig,thSig]=deal(sigma)`` for a 1-element input.

    Returns
    -------
    (ps_pol, th, frq)
        ``ps_pol`` -- ``(180/th_bin, max_rad+1)`` polar FFT.
        ``th`` -- angles (deg) labelling ``ps_pol``'s rows.
        ``frq`` -- FFT components labelling ``ps_pol``'s columns (``frq * (
        bin_size_cm/fft_size)`` is *not* wavelength directly -- see
        :func:`grid_fft_components` for the wavelength conversion actually
        used).

    .. warning::
       If ``max_rad == 0``, the final crop line (``psPol(:,1:maxRad+1)`` in
       the MATLAB, unconditional -- it does not participate in the
       ``if maxRad`` branch above it) collapses the output to radius bin 0
       only, rather than "no limit" as ``max_rad=0`` means in
       :func:`power_spectrum_2d`. This is a real inconsistency in the
       original between the two functions' treatment of the same 0 sentinel,
       faithfully preserved (well-defined, not a crash) rather than fixed --
       see ``FINDINGS.md``. Always call with a genuine positive ``max_rad``
       in practice (as :func:`grid_fft_components`/``fftVars.m`` do).

    .. note::
       Mean-aggregation uses an *explicit NaN fill* for bins no pixel ever
       lands in (``accumarray(..., @(x) mean(x), nan)`` -- NOT accumarray's
       default-0 fill). This is a genuinely different MATLAB convention from
       the visually-similar radial-profile aggregation in
       ``grid.py::fast_gridness`` (whose MATLAB caller never specifies a fill
       value, so it defaults to 0) -- confirmed by directly reading both
       ``.m`` sources rather than assumed from the similar shape. A bin can
       also come back NaN despite having contributing pixels, if ANY of them
       was itself NaN (typically from ``power_spectrum_2d``'s own
       ``max_rad`` circle-mask, whose boundary can disagree by up to ~0.5
       bins with this function's rounded-radius rebinning): MATLAB's plain
       ``mean()`` (not ``nanmean()``) propagates a single NaN to the whole
       bin, and this is reproduced by summing with ``np.bincount`` directly
       over the (possibly-NaN) pixel values, not by ignoring NaNs.
    """
    ps2d = np.asarray(ps2d)
    if ps2d.ndim != 2 or ps2d.shape[0] != ps2d.shape[1]:
        raise ValueError(
            f"to_polar_fft expects a square 2-D power spectrum, got shape {ps2d.shape!r}"
        )
    fft_size = ps2d.shape[0]
    cent = fft_size // 2  # 0-based centre -- matches power_spectrum_2d's fftshift convention

    max_rad = int(max_rad)
    if max_rad:
        if max_rad > cent:
            # MATLAB indexing ps2D(y+centCord, x+centCord) out of range would
            # hard-error here; Python's negative-index wraparound would
            # instead silently return wrong data -- raise explicitly (new
            # code's own robustness concern, not a MATLAB-fidelity issue).
            raise ValueError(f"max_rad ({max_rad}) exceeds the array's half-size ({cent})")
        x = np.arange(-max_rad, max_rad + 1)
        y = np.arange(-max_rad, 1)
    else:
        x = np.arange(fft_size) - cent
        y = np.arange(-cent, 1)

    rows = y + cent
    cols = x + cent
    ps2d_sel = ps2d[np.ix_(rows, cols)]  # top half of ps2D (+ max_rad crop)

    xx, yy = np.meshgrid(x, y)
    # y increases downward in image/FFT coordinates, but angle should be
    # measured the standard anticlockwise-from-+x way -- negate the row
    # offset before atan2 (matching the cart2pol convention already
    # established in grid_phase.py::cross_sac_props / field_shuffle.py).
    th = np.arctan2(-yy.astype(float), xx.astype(float))
    r = np.hypot(xx.astype(float), yy.astype(float))
    # abs() undoes IEEE-754 signed-zero: on the y==0 row, -yy is -0.0, and
    # atan2(-0.0, negative x) returns -pi (not +pi) -- the MATLAB comment's
    # "change -180 to 180 [no other -ve vals]" refers to exactly this; every
    # other value is already >= 0 here since y (hence -yy) never exceeds 0.
    th = np.abs(np.degrees(th))

    r_subs = _matlab_round(r / r_bin).astype(np.int64)  # 0-based radial bin
    n_th_bins = int(round(180.0 / th_bin))
    th_subs = np.mod(_matlab_round(th / th_bin).astype(np.int64), n_th_bins)  # 0-based angular bin
    n_r_bins = int(r_subs.max()) + 1

    flat_idx = (th_subs * n_r_bins + r_subs).ravel()
    values = ps2d_sel.ravel()
    n_bins = n_th_bins * n_r_bins
    sums = np.bincount(flat_idx, weights=values, minlength=n_bins)
    counts = np.bincount(flat_idx, minlength=n_bins)
    ps_pol = np.full(n_bins, np.nan)
    hit = counts > 0
    ps_pol[hit] = sums[hit] / counts[hit]  # NaN propagates through hit bins too (see docstring)
    ps_pol = ps_pol.reshape(n_th_bins, n_r_bins)
    ps_pol = ps_pol[:, : max_rad + 1]  # unconditional crop -- see the warning above

    th_out = np.arange(n_th_bins) * th_bin
    # frq labels the max_rad+1 retained columns; exact for r_bin=1 (the only
    # value the MATLAB author recommends/the caller ever uses -- see the
    # r_bin docstring above and FINDINGS.md for the r_bin != 1 caveat).
    frq_out = np.arange(max_rad + 1) * r_bin

    if sigma:
        if np.isscalar(sigma):
            r_sig = th_sig = float(sigma)
        else:
            sigma = tuple(sigma)
            r_sig, th_sig = sigma if len(sigma) == 2 else (sigma[0], sigma[0])
        ps_pol = _smooth_polar_fft(ps_pol, r_sig, th_sig)

    return ps_pol, th_out, frq_out


# --------------------------------------------------------------------------
# main_spatial_fft2.m -- Recording-level orchestrator
# --------------------------------------------------------------------------
def grid_fft_components(recording, electrode, cluster, params=None, n_shuffles=None, rng=None):
    """Significant 2-D FFT spatial-frequency components of a cell's ratemap.

    A ``Recording``-level redesign of ``main_spatial_fft2.m`` (like
    :func:`~barrylab.analysis.decoding.decode_crunch_data` and
    :func:`~barrylab.analysis.field_shuffle.field_shuffle_significance` --
    validate the component functions (:func:`power_spectrum_2d`,
    :func:`to_polar_fft`) individually; this orchestrator is intentionally
    **not** bit-exact, and inherits genuine shuffle randomness on top).

    Bins an *unsmoothed* ratemap (``smooth_ratemap(..., kern_size=1)`` is
    exactly a boxcar-of-width-1, i.e. a no-op smooth that still marks
    unvisited bins NaN -- matching ``make_smooth_ratemap(binXY,binSpk,1)`` in
    the MATLAB), gets its 2-D FFT and polar FFT, then compares the real polar
    FFT against ``n_shuffles`` **temporal** shuffles (spikes' position-sample
    indices circularly rotated by a random amount) to threshold out
    non-significant power before finding angular peaks.

    .. note::
       **Temporal-shuffle bound, a deliberate design choice, not a MATLAB
       transcription.** ``main_spatial_fft2.m``'s own shuffle
       (``sf_shuf_spk_ind``) draws ``randShft = ceil(rand()*(nPosPt-minShuf))
       + minShuf`` -- one-sided: it only guarantees the shift is *at least*
       ``minShuf`` samples, but says nothing about staying far from
       ``nPosPt`` itself. Since the shift then wraps via
       ``mod(spkInd+randShft, nPosPt)``, a draw close to ``nPosPt`` produces
       an almost-zero *effective* rotation after wraparound -- silently
       defeating the whole point of ``minShuf``, for a genuine (if narrow)
       slice of the draw range. Rather than reproduce this, ``rng`` here
       draws from a **two-sided** bound, ``[min_rot_samp, n_pos -
       min_rot_samp]``, exactly matching
       :func:`~barrylab.analysis.field_shuffle.field_shuffle_significance`'s
       already-established temporal-shuffle convention in this package (per
       the porting brief's explicit instruction to reuse that exact
       pattern), so every draw is genuinely bounded away from a
       near-zero-effective-rotation on both sides. Also, like that function,
       generalises the MATLAB's hardcoded 50 Hz (``f.minShuf*50``) to
       ``recording.pos_sample_rate``.

    Parameters
    ----------
    recording : Recording
    electrode, cluster : int
    params : FFTParams, optional
        Analysis parameters (defaults to ``FFTParams()``).
    n_shuffles : int, optional
        Override ``params.n_it`` (the real MATLAB default, 1000, is
        expensive -- pass a small number, e.g. 5-10, for quick/test runs).
    rng : int, numpy.random.RandomState, or None
        Seed or random state for the temporal shuffle.

    Returns
    -------
    dict with:
        ``ps2d`` -- the real cell's 2-D FFT power spectrum.
        ``ps_polar`` -- the real cell's (smoothed) polar FFT.
        ``ps_polar_thresholded`` -- ``ps_polar`` minus the shuffle's 95th
            percentile, with non-positive values set to NaN.
        ``theta_deg`` -- angles (deg) labelling ``ps_polar``'s rows.
        ``angular_profile`` -- ``ps_polar_thresholded`` collapsed (nanmean)
            over radius -- one value per angle.
        ``shuffle_max_power`` -- ``(n_shuffles,)`` each shuffle's max polar
            power (MATLAB ``max()`` ignores NaN by default -- reproduced with
            ``np.nanmax``, since the polar FFT legitimately contains NaN).
        ``crit95`` -- the shuffle distribution's 95th percentile
            (``np.nanpercentile``, matching MATLAB ``prctile``'s own
            NaN-removal behaviour).
        ``significant`` -- bool, whether any angular bin exceeded threshold
            (MATLAB's ``~any(angMean>0)`` early-return, redesigned as a flag
            on an always-populated dict rather than ``deal([])`` on 4
            separate outputs -- a deliberate, natural consequence of this
            port's dict-return redesign, not a MATLAB behaviour change).
        ``peak_angle_deg`` -- angle(s) (deg) of significant Fourier
            components (empty if not ``significant``).
        ``peak_wavelength_cm`` -- wavelength (cm) of each significant
            component's dominant radial frequency (shorter than the grid
            scale for a perfectly regular grid: ``wavelength = scale *
            cos(30deg)``).
    """
    if params is None:
        params = FFTParams()
    n_shuffles = params.n_it if n_shuffles is None else n_shuffles
    rng = _as_random_state(rng)
    max_rad = params.max_rad

    elec = recording.electrode(electrode)
    spk_idx = elec.spike_pos_index(cluster)
    n_pos = recording.n_pos

    bin_occ = occupancy_map(recording, params.bin_size_cm)
    bin_spk = spike_map(recording, spk_idx, params.bin_size_cm)
    rm, _, _ = smooth_ratemap(
        bin_occ, bin_spk, 1, "boxcar"
    )  # kern_size=1 -> unsmoothed, NaN-marked

    th_sig_bins = params.th_sig / params.th_bin  # smoothing sigma: degrees -> angular bins
    r_sig_bins = params.r_sig / params.r_bin  # smoothing sigma: FFT components -> radial bins

    ps2d = power_spectrum_2d(rm, params.fft_size, params.smth_2d_sig, max_rad)
    ps_pol, th, _ = to_polar_fft(
        ps2d, max_rad, params.th_bin, params.r_bin, (r_sig_bins, th_sig_bins)
    )

    min_rot_samp = int(round(params.min_shuf_s * recording.pos_sample_rate))
    rot_low, rot_high = min_rot_samp, n_pos - min_rot_samp

    shuf_max_power = np.full(n_shuffles, np.nan)
    for h in range(n_shuffles):
        rot_amount = rng.randint(rot_low, rot_high + 1)
        rot_spk_idx = (spk_idx + rot_amount) % n_pos
        shuf_bin_spk = spike_map(recording, rot_spk_idx, params.bin_size_cm)
        shuf_rm, _, _ = smooth_ratemap(bin_occ, shuf_bin_spk, 1, "boxcar")
        shuf_ps2d = power_spectrum_2d(shuf_rm, params.fft_size, params.smth_2d_sig, max_rad)
        shuf_ps_pol, _, _ = to_polar_fft(
            shuf_ps2d, max_rad, params.th_bin, params.r_bin, (r_sig_bins, th_sig_bins)
        )
        shuf_max_power[h] = np.nanmax(shuf_ps_pol)  # MATLAB max() ignores NaN by default

    crit95 = np.nanpercentile(shuf_max_power, 95)  # MATLAB prctile() also ignores NaN

    ps_pol_thresh = ps_pol - crit95
    ps_pol_thresh = np.where(ps_pol_thresh <= 0, np.nan, ps_pol_thresh)
    # A row can legitimately be all-NaN (no radial component in that angular
    # bin survived thresholding) -> nanmean warns "Mean of empty slice" and
    # correctly returns NaN for it; likewise the NaN vs. 0 comparison just
    # below is expected to hit NaN rows. Both are benign/expected here
    # (matching the same all-nan-slice suppression pattern in lfp.py).
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        ang_mean = np.nanmean(ps_pol_thresh, axis=1)

    result = {
        "ps2d": ps2d,
        "ps_polar": ps_pol,
        "ps_polar_thresholded": ps_pol_thresh,
        "theta_deg": th,
        "angular_profile": ang_mean,
        "shuffle_max_power": shuf_max_power,
        "crit95": crit95,
        "significant": False,
        "peak_angle_deg": np.array([]),
        "peak_wavelength_cm": np.array([]),
    }
    with np.errstate(invalid="ignore"):
        any_significant = np.any(ang_mean > 0)
    if not any_significant:
        return result

    snd_peak_thresh = np.nanmax(ang_mean) * params.snd_peak_coef  # MATLAB max() ignores NaN
    # find_peaks can't handle NaN neighbours sensibly (comparisons with NaN
    # are always False, which would suppress detection of a genuine peak
    # adjacent to a NaN gap) -- substitute -inf, which can never register as
    # a peak and behaves numerically sanely as a neighbour either side.
    ang_finite = np.where(np.isnan(ang_mean), -np.inf, ang_mean)
    peak_idx, _ = find_peaks(
        ang_finite, height=snd_peak_thresh, distance=max(1.0, params.min_peak_dist / params.th_bin)
    )
    if peak_idx.size == 0:
        return result

    # For each significant angular peak, the dominant radial frequency
    # component (MATLAB max()-ignores-NaN again -- nanargmax). A returned
    # angular peak is guaranteed to have >=1 finite (>0) entry in its row of
    # ps_pol_thresh (that's what made angMean>0 there), so nanargmax cannot
    # hit an all-NaN row.
    peak_rows = ps_pol_thresh[peak_idx, :]
    peak_freq_comp = np.nanargmax(peak_rows, axis=1)  # 0-based == MATLAB's (peakFComp-1)
    with np.errstate(divide="ignore"):
        wavelength_cm = params.fft_size / ((1.0 / params.bin_size_cm) * peak_freq_comp)

    result["significant"] = True
    result["peak_angle_deg"] = th[peak_idx]
    result["peak_wavelength_cm"] = wavelength_cm
    return result
