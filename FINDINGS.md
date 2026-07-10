# Findings from porting `Universal_Matlab`

Notes gathered while porting the DACQ loader + rate map + Skaggs information and
validating against the shipped `gridCell1838_t2c3_14-12-10.mat` fixture. Nothing
here blocks the port; a couple of items are decisions for the PI.

## 1. The example fixture is older than the code (2 mismatches, both benign)

The reference `.mat` validates the loader chain perfectly, but two steps do not
bit-match because the fixture predates the current MATLAB:

**(a) Head-position combination — DECISION NEEDED.**
The current `postprocess_pos_data.m` combines the two LEDs with a *tracking-fraction
weighted* mean (a change AJ added 2010-07-14 to be robust to a poorly tracked
light). The fixture's `data.pos.xy` was made with a plain **midpoint** of the two
LEDs — our port matches that to 1e-13. For this trial the back LED is tracked only
53% of the time, so the two conventions differ by up to ~8 px (~1 bin).

- `barrylab` implements both via `load_dacq(..., head_pos_weighting=...)`.
- Default is `"tracking"` (faithful to the current repo). `"midpoint"` reproduces
  the fixture / legacy behaviour.
- **Question for the PI:** which is canonical going forward? The weighted mean is
  more robust to a lost light but biases the estimated head position toward the
  better-tracked LED rather than the geometric head centre.

**(b) Rate-map smoothing — fixture used a "precursor".**
The example README says `smthRm` was made with "make_smooth_ratemap.m (or its
precursor)". Empirically the fixture matches a **Gaussian (σ≈2 bins, r=0.98)**, not
the current **5-bin boxcar** (`easy_make_ratemap.m`, r=0.94). Our port faithfully
reproduces the *current* boxcar `make_smooth_ratemap.m`; the binning underneath is
proven bit-exact (brute-force cross-check + exact totals). So the residual is the
fixture's smoothing, not our code.

**Consequence / recommendation.** The shipped `.mat` is enough to validate the
whole loader chain (which it does, to machine precision) but not the smoothing or
the default head position. Regenerate fresh fixtures from the *current* MATLAB with
`validation/generate_fixtures.m` and commit them, so every ported function has a
current golden-master test.

## 2. Minor issues spotted in the MATLAB (candidates to fix when porting fully)

- **`bin_data.m` line ~182: `if nargoout==2`** — `nargoout` is a typo for
  `nargout`. In MATLAB this errors if that branch is ever reached (the `pxd`
  path with two outputs), so that branch is effectively dead/broken. Not ported.
- **`make_smooth_ratemap.m` line ~77: `binSpk = lower(binSpk)`** — `lower()` on a
  numeric array is a no-op; harmless but clearly a leftover. Dropped in the port.
- **`read_DACQ.m` tetrode-count loop** — `ntetrodes_used` is incremented for all 16
  channels regardless of `collectMask`, so the variable does not mean what its name
  says. It is unused downstream; the port derives the tetrode list from
  `collectMask_*` directly.
- **`histnd.m` out-of-range warning** — the ">10% bad points" test sums bin
  *indices* (`sum(bad_low_points)`) rather than counting points, so the warning
  threshold is not what it appears. Cosmetic; the binning itself is correct.
- **`key_value` prefix matching** — the default (non-`exact`) match is by string
  prefix, so a key like `EEG_ch_1` also matches `EEG_ch_10..16`. The MATLAB callers
  pass `'exact'` where this matters; the port preserves that.
- **`multiGridness.m` line ~163: `innerRadi`** — the expanding-gridness inner radius
  is `mean(sqrt(sum(closestPeaksCoordCentral.^2)))/2`. The `sum` has no dimension
  argument, so it sums *down columns* (over peaks) rather than per-peak
  (`sum(...,2)`), which does not match the stated intent ("half of mean distance to
  the six central peaks"). Likely a latent bug. The port replicates the MATLAB
  behaviour so numbers match; flag for the PI whether to fix it (only affects the
  expanding gridness, which maximises over radii anyway).

## 3. Head-direction + place-field stats port, v0.2/v0.3 — dependency check

Prompted by a direct request to audit cross-function dependencies before reusing
ported code (not just port each function in isolation):

- **`led2v1_dir_offset.m` calls the *wrong* circular mean — genuinely worth
  knowing.** `cb_circMean.m`'s own header says it was written "to replace
  unreliable circ_mean", yet `led2v1_dir_offset.m` still calls the toolbox
  `circ_mean`, not `cb_circMean`. Checked via `grep -rl cb_circMean` across the
  whole MATLAB repo: `cb_circMean` has **no callers at all** outside its own
  file — every consumer, including this one, was written (or left) against the
  old function. This looks like an inherited inconsistency rather than a
  deliberate choice. Practically: the two functions compute the *same* mean
  angle (both reduce to `atan2(sum sin, sum cos)`), but `circ_mean` returns it
  in `(-π, π]` (MATLAB's `angle()` convention) while `cb_circMean` returns
  `[0, 2π)`. `led2v1_dir_offset` needs the signed range (a small negative bias
  must stay negative, not wrap to ~355°). **Fix:** `circular_mean()` now takes
  an explicit `wrap={"0_2pi", "-pi_pi"}` parameter (default `"0_2pi"`, matching
  `cb_circMean.m` exactly for bit-for-bit continuity with existing callers);
  `led2v1_dir_offset` explicitly requests `wrap="-pi_pi"`. This turns a
  silent-corruption risk into a documented, tested choice at each call site —
  see the docstrings in `analysis/direction.py`.
- **`pcFieldStats.m`'s `nansum(rmS)==0` check works, but only by luck.**
  `rmS` is a 2-D matrix; `nansum` without a dimension argument sums down
  columns, returning a row vector, not a scalar — so this reads like a bug
  (intent was clearly `rmS(:)`, i.e. the whole-matrix sum). MATLAB's `if` on a
  vector requires *every* element to be truthy, which is exactly equivalent to
  "grand total is zero" **only because firing rates are non-negative** (if
  every column-sum is 0 and all values are ≥0, every value must be 0, and vice
  versa). So the result is correct, just fragile — ported as the clear,
  equivalent `np.nansum(rm) == 0` (flattened) rather than replicating the
  column-wise quirk.
- **`get_true_PPM_boundingBox.m` under-crops by one row and column — a real
  off-by-one, faithfully reproduced.** `ijWidth = boundingBox([4 3]) - 1` then
  `iRange = ijCorner(1):ijCorner(1)+ijWidth(1)-1` has one `-1` too many: it
  should enumerate `r1:r2` (all rows of the bounding box) but actually gives
  `r1:r2-1`, dropping the last row (and, symmetrically, the last column). This
  doesn't affect the returned PPM value (computed from the filled *area*, not
  from these ranges), only the crop ranges returned for callers who use them to
  window the ratemap. Ported faithfully — see `true_ppm_bounding_box()`'s
  docstring/comment in `analysis/spatial.py` — since goal 1 is matching
  results, not silently fixing bugs; worth a PI decision on whether to correct
  it once parity is established.

## 4. Grid-cell analysis port (SAC + gridness), v0.2

Ported `xPearson.m` -> `spatial_autocorrelogram`, `multiGridness.m` -> `gridness`
(standard + expanding), and `sacProps.m` -> `grid_properties`. The MATLAB
image-processing calls (`imregionalmax`, `bwlabel`, `regionprops`, `imrotate`) are
reproduced with `scipy.ndimage` (no scikit-image dependency). `imregionalmax` is
approximated as a 3×3 maximum-filter equality restricted to positive bins; on real
(sharp-peaked) SACs this finds the same peaks, but could differ from MATLAB on
pathological plateaus.

**Validation (structural, not yet bit-exact).** No golden gridness fixture ships
with the toolbox, so the port is validated by behaviour: the SAC is exact by
construction (central bin = 1, point-symmetric), a synthetic hexagonal grid scores
high gridness, and — the strong check — the port reproduces the example README's
ground truth: it scores the documented grid cells (tetrode 2 cells 1 and 3) at
gridness 1.25 and 0.85 and the non-grid cells (2c4, 2c5) at −0.20 and 0.19. Grid
scale comes out at 40 cm. To lock these bit-exact, `validation/generate_fixtures.m`
now also dumps `sac`, `gridness_std/exp`, `grid_scale/orientation/mean_r` from the
current MATLAB; add a regression test asserting `np.allclose` once a fixture exists.

## 5. Place x direction ML model (pxd), v0.3

Ported `pxd.m` (the alternating-maximum-likelihood factorial fit) and
`wrapper_pxd.m` (raw-data binning + smoothing wrapper) as `analysis/pxd.py`.
Also added a fresh general N-D `histnd` (`analysis/_binning.py`) since this is
the first analysis needing true 3-D (direction x posY x posX) binning — cross-
checked against the existing validated 2-D binning inside `maps.py` and found
**exactly identical** (max diff 0.0), so the new implementation is trustworthy.

**A real bug in `wrapper_pxd.m`, faithfully reproduced.** The spatial bin range
is built as `[max(y), max(y)]` instead of `[max(y), max(x)]` — a copy-paste
error (the corresponding *min*-value row correctly uses `min(y)` then
`min(x)`). For a non-square tracked window this mis-scales the X extent used
for binning. Per goal 1, replicated exactly (see `place_by_direction()`'s
docstring) rather than silently fixed.

**Smoothing derivation worth recording.** `wrapper_pxd.m` smooths its two
output marginals with ad hoc inline code (edge-corrected boxcar for the
position map, plain circular boxcar for direction) rather than calling the
already-existing `make_smooth_ratemap.m`. Worked out that both are exactly
expressible as calls to our already-ported `smooth_ratemap()`: pass the
binary visited-mask as `bin_pos` and `p * visited_mask` as `bin_spk` and
`smooth_ratemap`'s own edge-correction reduces to precisely the same
edge-corrected mean the wrapper computes by hand (the self-normalising ratio
`smooth(visited)/smooth(visited) = 1` cancels out the "rate" division,
leaving just the edge-corrected smooth of `p`). Confirmed this algebraically
before implementing — reuses validated code instead of duplicating the
smoothing logic a second time.

**Validation.** Structural + cross-consistency, no golden fixture available
for this one: the recovered place map correlates r=0.94 with the
already-validated plain rate map for the same cell (r1838 t2c3); the
recovered direction map correlates r=0.85 with the already-validated plain
head-direction tuning curve. Also a synthetic test with known ground-truth
place x direction factors, recovered with r>0.8/0.7 respectively by the ML
fit, plus an exact check of the model's own normalisation invariant (total
predicted spikes matches actual, per marginal).

## 6. Bayesian decoding, v0.4

Ported `decode_calcPost.m` -> `bayesian_posterior` (naive ML Poisson population
decoder), `decode_getSpkCnt.m` -> `spike_counts_per_window` (windowed, optionally
overlapping/Gaussian-smoothed spike counting), `decode_processPost.m` ->
`decode_maximum_likelihood` (ML decode + accuracy metrics), and
`decode_calcBayesPost.m` -> `bayesian_posterior_with_prior` (full Bayesian
decoder with a spatial occupancy prior and a sequential movement/"closeness"
prior, Zhang et al. 1998 eq. 36). `decode_crunchData.m` is *not* a port but a
fresh `crunch_decode()` orchestrator redesigned around `Recording` (train/test
split, build rate maps, decode) — validate the four functions above
individually, this convenience wrapper is intentionally not bit-exact.

**Two bugs found and fixed in my own new code (not MATLAB-fidelity issues —
`crunch_decode`/`bayesian_posterior_with_prior` are original code, so these are
just implementation bugs, fixed directly rather than "faithfully reproduced"):**
1. `crunch_decode`'s first draft computed occupancy over the *whole* session
   while spikes were correctly restricted to the training split, diluting the
   rate-map denominator with test-period dwell time. Fixed by reusing
   `spike_map` (which already bins an arbitrary position-index subset) for a
   train-restricted occupancy count too, instead of the whole-session
   `occupancy_map`.
2. `bayesian_posterior_with_prior`'s per-step peak-finding used plain
   `np.argmax`, not `np.nanargmax`. NaN in a numpy array is comparison-broken
   (`nan > x` is `False`), so with real (partially-unvisited) rate maps a
   sizeable NaN-contaminated posterior would make `argmax` return the *first*
   NaN cell as the "decoded" peak — corrupting every subsequent step of the
   sequential movement-prior filter. MATLAB's `max()` ignores NaN by default,
   so this was a genuine divergence from the source, not a faithful port of
   anything. Confirmed the fix: before, per-timestep posterior sums drifted
   away from 1 in the presence of NaN rate-map bins; after, they're exact.

**A subtler finding worth documenting so it isn't later mistaken for a bug.**
`decode_processPost.m`'s docstring claims ties in the ML decode are broken "at
random" (`tmpInd=tmpInd(randperm(length(tmpInd),1))`), but MATLAB's `max()`
already collapses to a single (first) index before that line runs, so the
"random" selection is drawing from a length-1 list — a no-op. Ties always
resolve to the first maximum in MATLAB's column-major linear order. Reproduced
faithfully via `order='F'` (not genuine random tie-breaking) per goal 1.

**Validation.** `bayesian_posterior` is bit-exact against a hand-computed
Poisson PMF (`scipy.stats.poisson.pmf`, matches to 1e-6) and normalises exactly
to 1. `spike_counts_per_window` matches hand-counted windows (including
overlap) exactly, and Gaussian smoothing preserves total spike mass. The
decisive check for `decode_maximum_likelihood`/`bayesian_posterior`: decoding
**synthetic data with clean, non-periodic place fields and known ground
truth** recovers position far better than chance (mean error 2.25 bins vs a
~9.4-bin chance level on a 20x20 grid). Real r1838 data (via `crunch_decode`)
decodes only modestly better than a naive chance baseline when using the
documented grid cells (t2c1/t2c3) — this is very likely genuine
neuroscience, not a bug: grid cells have multiple periodic firing fields, so
their population activity is intrinsically ambiguous about *absolute*
position (a well-documented real limitation), unlike the clean single-field
synthetic case which decodes cleanly. The synthetic test is the real
correctness proof; the real-data numbers are reported for context, not as a
pass/fail bar.

## 7. Phase precession: circular-linear regression

Ported `circ_regress_cb.m` -> `circular_linear_regression()`
(`analysis/phase_precession.py`): fits a line relating a linear variable
(in-field position) to a circular one (theta phase) by maximising the mean
resultant length of the phase residual over slope, reports a
circular-circular correlation (Kempter et al. 2012, via the already-ported
`circular_mean`), and estimates significance by re-fitting on `k` shuffles of
the phase relative to position. Confirmed this correctly calls the lab's own
`cb_circMean` (not the buggy toolbox `circ_mean` — see the earlier finding
about `led2v1_dir_offset.m`; this one has no such landmine).

**Validation oracle resolved a real puzzle, not a bug.** The MATLAB docstring
embeds Fisher (1993) dataset B.20, claiming it should give `R=-0.316,
p=0.057` when run through the function. Running the full port on that data
gave `R=-0.4985`, a real mismatch worth investigating rather than dismissing.
Re-reading the comment closely: it says the book's R value comes from "a
circ-circ correlation **after taking 2*atan(t)**" — i.e. Fisher's textbook
answer correlates a *fixed* transform of the linear variable (`2*atan(x)`)
against the phase, which is a fundamentally different computation from what
`circ_regress_cb` actually does (which fits an *optimal* slope for
`mod(slope*x, 2*pi)`, not the fixed `2*atan(x)` transform) — the two are not
expected to agree except by coincidence. Confirmed this interpretation
directly: `kempter_corr(2*atan(x), t) = -0.3160`, matching the book value to
4 decimal places. So the underlying circular-circular correlation formula
(and the `circular_mean` it depends on) is now bit-exact validated against an
independent textbook oracle; the full pipeline is separately validated by
recovering a known, strong synthetic phase-precession slope accurately
(within 5%) with a highly significant permutation p-value, while
random/uncorrelated data reliably gives a far less significant result by
comparison.

**Not yet ported:** the `lineTraj_*` group (replay/line-trajectory fitting,
~500 lines across 4 files, with a shipped `pMat_data.mat` oracle per
`PORTING_PLAN.md`) — deferred to cover more breadth across v0.5/v0.6 given
time constraints; still on the plan.

## 8. Elliptical grid correction, v0.5

Ported `gridEllipse_fit.m` -> `fit_grid_ellipse()` (fits an ellipse to a grid
cell's six central SAC peaks via least-squares conic fitting — a MathWorks
File Exchange algorithm embedded as a subfunction in the original, ported
faithfully as `_fit_ellipse()`) and `gridEllipse_correct.m` ->
`regularise_elliptic_sac()` (rotates/rescales/rotates-back an elliptical SAC
toward circular, so gridness can be recomputed on the corrected map;
`regularise_eliptic_grid.m` is a byte-identical MATLAB duplicate, ported
once). Reused the existing `grid.py::_sac_peaks` peak-finder rather than
duplicating it a third time (new code, not a refactor of shipped code).

**A synthetic-test convention puzzle, resolved.** A first synthetic test
constructed a known ellipse with a standard (CCW) rotation matrix and got the
correct semi-axes back but the *wrong-signed* orientation. Investigated
rather than adjusted the test: `gridEllipse_fit.m`'s own plotting subcode
uses `R = [cos_phi sin_phi; -sin_phi cos_phi]` — the **transpose** of a
standard CCW rotation matrix — so the algorithm's own internal convention for
`phi` is CW-flavoured relative to a naive assumption. Reconstructing the
synthetic test with MATLAB's own `R` convention recovers `a`, `b`, **and**
`phi` bit-exactly (`atol=1e-6`) — confirming the port is correct and the
mismatch was purely in the test's assumed convention, not the code.

**A units mismatch inherited from the original, documented rather than
silently "fixed".** `fit_grid_ellipse()`'s `orientation` is in radians, but
`regularise_elliptic_sac()`'s `orient_deg` is in degrees — this split exists
in the MATLAB too (`gridEllipse_fit.m` returns radians; `gridEllipse_correct.m`
expects degrees; `brandonRegEcc.m` explicitly `rad2deg`-converts when chaining
them). Kept identical to preserve the exact contract, with an explicit
docstring note and conversion example so it isn't a silent trap.

**A real off-by-one fixed (new code, not a MATLAB-fidelity issue).**
`gridEllipse_correct.m`'s final centre-crop (`sizeDif=round(...)`) can land
one pixel short after an odd rotated-size difference, because MATLAB's
`round()` rounds half-away-from-zero while a naive numpy port would round
half-to-even (see finding on `multiGridness.m`'s `innerRadi` for the same
class of MATLAB/numpy rounding divergence). Since the MATLAB code's own
comment states its explicit intent ("keep only the central portion... so
that it matches the original size"), implemented the crop to directly slice
to the exact target shape rather than replicate the rounding ambiguity —
this fulfils the code's own documented purpose rather than changing what any
value means.

**Validation.** `_fit_ellipse` is bit-exact against a known synthetic
ellipse (using the algorithm's own rotation convention, `atol=1e-6`) and
correctly recovers a circle (`a==b`) for circularly-arranged points. On real
r1838 grid cells, `fit_grid_ellipse` returns sensible eccentricity
(0.4-0.55), scale, and orientation; `regularise_elliptic_sac` runs cleanly
(correct output shape, low NaN fraction) and gridness stays finite and
sane after correction. The MATLAB source itself flags this whole method as
"not fully tested... requires further testing before publishing results" —
correction does not reliably *improve* gridness on both test cells (it did
for t2c1, was roughly flat for t2c3), consistent with that original caveat
rather than indicating a port defect.

## 9. LFP loading infrastructure, v0.6a — the theta/LFP unlock

Added ``.eeg``/``.egf`` reading to the DACQ loader (`load_dacq(...,
load_lfp=True)`) and a `Recording.lfp` field, unblocking the whole
theta/LFP/phase-precession analysis family. Ported `eeg_convert_to_volts.m`
faithfully. `.eeg` (250 Hz, signed 8-bit) reuses the existing
`data_start`/`data_end` framing infrastructure directly (`read_axona_file`),
since it turns out to use exactly the same marker convention as `.pos` and
the tetrode files.

**`.egf` (4.8 kHz, signed 16-bit) also uses that same marker framing —
simpler and more robust than the MATLAB reader's approach.** `getegf.m`
doesn't search for `data_start` at all; it hardcodes "read exactly 8 header
lines, extract sample rate, read 2 more lines, extract sample count, then
skip forward 10 bytes" — which only works if the header always has exactly
10 lines before the marker. Confirmed against the real example `.egf`
header that this coincides with the marker-search approach for
standard-format files, so reused the same robust `read_axona_file` +
`little_endian_view` (a new helper — `.egf` is the one DACQ file that is
little-endian; everything else is big-endian) infrastructure instead of
replicating the fragile fixed-line-count parsing.

**Deliberately did not port `getegf.m`'s decimation step.** MATLAB's
`decimate()` uses a specific Chebyshev IIR filter design; reproducing it
bit-exact in another language/library isn't something `scipy.signal.decimate`
(a different default filter) can guarantee, and goal 1 is "same results," not
"a plausible-looking approximation." `_read_egf`/`load_dacq` return the raw,
undecimated 4.8 kHz signal — an honest, unambiguous choice — rather than
silently offer a resampled signal that looks right but may not numerically
match MATLAB's output. A caveated, explicitly-approximate decimation helper
can be added later if a bit-exact comparison becomes necessary.

**Validation.** Structural + physiological plausibility (no MATLAB fixture
for LFP data yet): loaded sample counts imply a duration that matches the
session's position-derived duration to within 1 second for every channel on
r1838 (3 EEG channels). A crude FFT of the first 60s of the primary channel
peaks at 9.58 Hz within the 4-12 Hz band — squarely in the classic
hippocampal theta range for a foraging rat, a meaningful physiological
sanity check that the byte layout and endianness are being read correctly
(a parsing error would far more likely produce noise or an implausible peak
than a textbook theta frequency).

## 10. LFP power spectrum, v0.6b

Ported `eeg_power_spectra.m` -> `power_spectrum()` (`analysis/lfp.py`): an
FFT-based power spectrum with Gaussian smoothing, theta-band peak detection,
and signal-to-noise, optionally restricted to specific position samples (e.g.
only fast-running periods). Derived the `pos2use`-to-EEG-sample index mapping
fresh in 0-based terms rather than transcribing MATLAB's 1-based offset
formula, to avoid an off-by-one.

**`eeg_instant_freq_power.m` (bandpass + Hilbert instantaneous phase/frequency)
and `detect_ripples.m` were deferred at the time this section was written;
both are now ported — see §20 for the full write-up.**

**Validation.** Synthetic sinusoid embedded in noise: recovered peak
frequency within 0.5 Hz and a clear, high signal-to-noise ratio. On real LFP
(r1838 ch1, the `.egf` channel loaded in v0.6a): peak at 9.14 Hz (whole
session) / 9.38 Hz (a position-restricted subset), both cleanly in the
6-12 Hz theta band, with signal-to-noise of 14-16x — a textbook-clean
hippocampal theta signature for a foraging rat, and a strong end-to-end
sanity check spanning the new LFP loader (v0.6a) through this analysis.

## 11. Watershed field segmentation, v0.7 (partial)

The PI chose to add `scikit-image` (a new optional `[fieldshuffle]` extra —
scipy has no watershed transform) rather than vendor one or substitute a
simpler peak+threshold scheme. Ported the field-**segmentation** stage:
`watershed_fields()` (`watershed(-wsRm)` from `fs_masterAnalysis.m`, followed
by `fs_fixWatershed.m`'s ridge-bin reassignment) and `field_bin_locations()`
(`fs_binMembLoc.m`: each field's peak and every other bin's polar offset from
it). Lives in `analysis/field_shuffle.py`, not eagerly imported by the
top-level package (matching the existing lightweight-by-default pattern for
`[neuropixels]`/`[nwb]` — a DACQ-only user shouldn't need scikit-image).

**Confirmed a real API gap in scikit-image's `watershed()` before writing
any port code, rather than assume.** MATLAB's `watershed()` leaves ridge
lines unlabelled (0) by design — `fs_fixWatershed.m`'s whole purpose is
reassigning those. `skimage.segmentation.watershed()` defaults to
`watershed_line=False`, which assigns *every* pixel to a basin with no
ridge left to reassign — silently making the whole two-stage MATLAB design
pointless if used as-is. Verified with a tiny synthetic test that
`watershed_line=True` restores the 0-labelled-ridge behaviour MATLAB
depends on, and used `connectivity=2` to match MATLAB's default 8-connectivity
for 2-D images (skimage's default `connectivity=1` is 4-connected).

**A real bug in my own new code, found and fixed via validation, not
MATLAB-fidelity.** A first draft additionally forced every unvisited-position
bin to label 0 *after* running watershed, which doesn't match MATLAB (which
only zeroes the *rate map* before watershedding and lets the transform
itself decide those bins' fate) — this corrupted the ridge-only semantics
`fs_fixWatershed` depends on and crashed `field_bin_locations` on an
all-NaN field slice. Fixed by removing the post-hoc overwrite; validated by
confirming every bin then belongs to exactly one field (`sum(field sizes) ==
total bins`) and that `_fix_watershed`'s reassignments are always computed
against the *original* (pre-reassignment) labels — verified with a
purpose-built test where a naive in-place implementation would pick the
wrong neighbour.

**Not yet ported: `fs_shufFields.m` (the actual randomised field-reallocation
shuffle) and the `fs_masterAnalysis.m`/`identify_grid_cell.m` orchestration.**
These need substantially more careful treatment than the segmentation stage
here — `fs_shufFields.m` is a genuinely intricate iterative algorithm
(random field-centre reallocation, per-field random orientation, polar-order
bin placement, round-robin interleaving across fields), and the
orchestration additionally depends on `fs_fastGridness.m` (a *different*,
Langston-2011-specific expanding-gridness variant — not the same as the
already-shipped `gridness()`'s expanding measure, which follows
`multiGridness.m`'s simpler heuristic) and `brandonRegEcc.m` (not yet
ported). Rather than deliver a rushed, under-validated version of the most
algorithmically complex code in the toolbox, this is left for a dedicated
pass. Still tracked in `PORTING_PLAN.md`.

**Validation.** `_fix_watershed` unit-tested directly, including a
purpose-built case distinguishing "uses frozen original labels" from a naive
in-place implementation. `watershed_fields` + `field_bin_locations` validated
on synthetic two-peak data (correct peak recovery, correct field-size
accounting) and on the real r1838 grid cell (16 fields found via unmarked
watershed — plausible over-segmentation for a marker-free watershed, a known
general property of the method itself rather than a port defect, though not
cross-checked against MATLAB directly; every one of the map's bins accounted
for in exactly one field).

## 12. The local MATLAB clone was stale — a real upstream bug fix found and applied

While adding a note to `Universal_Matlab`'s README linking to this port (see
its repo), discovered the local clone used for the whole port so far was 7
commits behind the actual GitHub `master`. Diffed every changed file against
what's already been ported (`decode_calcPost.m` was the only overlap; the
other new/changed files — `decode_cosSim.m`, `decode_calcPost_twoPhoton.m`,
`makeRM_kilosort.m`, `circ_circ_corr.m` extraction, a `default_vars.m`
directory move — aren't relevant to anything ported yet).

**`decode_calcPost.m` was fixed upstream in 2020** (commit `d12abbc`,
"changed the min value that is add to prevent values being 0"): the
zero-avoidance epsilon `smallV` changed from `eps^8` (~4.6e-124 — so
vanishingly tiny it doesn't actually prevent the zero-probability problem
its own comment describes, since adding it to a rate map bin leaves the
bin numerically indistinguishable from exactly 0) to plain `eps`
(~2.22e-16 — still tiny, but numerically meaningful for its stated purpose).
Since this is the *author's own* later correction to the gold-standard
repo, not a discrepancy between my port and some external reference, this
is a case for actually updating the port (unlike goal-1's "preserve MATLAB
bugs faithfully" rule, which is about not silently diverging from the
authoritative implementation — here the authoritative implementation moved,
and the old value was never the intended one). `bayesian_posterior`'s
`_SMALL_V_CALC_POST` updated from `eps**8` to `eps` to match.

**A genuine inconsistency in the current MATLAB itself, left as-is.**
`decode_calcBayesPost.m` has the *identical* originating comment and
originally-identical `smallV=eps.^8` line as `decode_calcPost.m` — but the
2020 fix was applied only to `decode_calcPost.m` and never propagated to its
near-twin. Checked directly against the current upstream file to confirm
this isn't an artifact of my own stale clone. Since `decode_calcBayesPost.m`
genuinely still has the old value in the actively-maintained repo,
`bayesian_posterior_with_prior`'s `_SMALL_V_BAYES_POST` is kept at `eps**8`
to match its actual current source — the two Python constants are now
deliberately different, one per MATLAB function, rather than a single shared
value, because the two MATLAB functions have themselves diverged.

**Practical note for future sessions:** the local MATLAB clone can drift
behind the real upstream between sessions (this lab's other members/PI
actively commit to it). Worth an occasional `git fetch && git log
HEAD..origin/master` check on `Universal_Matlab` before trusting it as the
current ground truth, especially before porting a function whose file has
seen no local activity in a while.

## 13. Line-trajectory replay scoring (`lineTraj_*`), v0.4 (deferred piece)

Ported the previously-deferred `lineTraj_*` group (`analysis/line_trajectory.py`):
`best_line` (`lineTraj_decode.m`), `reget_line_band` (`lineTraj_reget_lnBnd.m`),
`line_residual` (`lineTraj_get_res.m`). Finds the best constant-speed straight
line through a decoded posterior matrix (space x time) — used to score
candidate replay trajectories.

**Biggest risk was `filter2`/`conv2`/the inline `il_filter2` kernel
anchoring**, flagged in advance in `PORTING_PLAN.md`. Verified the two
translations numerically against a brute-force reimplementation *before*
writing the port (not just assumed from the MATLAB docs' stated equivalence):
`il_filter2(a,b) = conv2(b, rot90(a,2), 'valid')` translates to
`scipy.signal.convolve2d(b, np.rot90(a,2), mode='valid')`; and
`filter2(kern, X)` with a symmetric all-ones column kernel in `'same'` mode
(used for the line's band) is exactly a rolling sum along axis 0, i.e.
`scipy.ndimage.convolve1d(X, np.ones(k), axis=0, mode='constant', cval=0.0)`.
Both matched a brute-force loop implementation exactly.

**No MATLAB runtime available on this machine**, so unlike every other ported
function so far there's no way to generate independent reference numbers from
the shipped `Decoding/pMat_data.mat` fixture to compare against. Used the
oracle `lineTraj_testScript.m` itself relies on instead: rebuilding the best
line's band from its own reported gradient/intercept
(`reget_line_band`) and rescoring it against the pMat (`line_residual`) must
reproduce the original residual exactly. Ran this across all 230 events in
`pMat_data.mat` (118 yield a fit; the other 112 have no valid data bins, per
the fixture's own `READ ME.txt` extraction recipe) — all 118 round-trip to
within 1e-15 (float64 machine precision). A bug in any of the three functions,
or in the shared offset/index arithmetic between them, would very likely have
broken this, so it's a strong (if not literally MATLAB-comparative) check.
**This is a validation gap worth closing if anyone has MATLAB access:** run
`lineTraj_testScript.m` on `pMat_data.mat` and compare the printed `bstRes`
values event-by-event against this port's `best_line(...)[0]`.

**A genuine MATLAB quirk found and faithfully replicated (not fixed):**
`lineTraj_decode.m`'s own comment says "if there is no data... bestRes will
be 0", and it has an explicit `if bstRes == 0` branch meant to null out
`bstSpd`/`bstYInterc`/`bstTstLnBnd` to `nan` for that case. But for a
literally all-zero `pMat`, `bstRes` is computed as `0/sum(pMat) = 0/0`, which
under IEEE 754 is `NaN`, not `0` — so the `== 0` check silently never fires
for the very case the comment describes. `bstRes` still ends up `NaN` (any
caller checking `isnan(bstRes)` behaves correctly), but `bstSpd`/`bstYInterc`
are left as whatever the first (arbitrarily tied) candidate gradient happened
to be, rather than also being `nan`. Reproduced exactly in `best_line`
(`line_trajectory.py`), pinned by
`test_best_line_all_zero_pmat_leaves_spd_uninitialised_nan_quirk`. Mostly
dormant in practice since real callers (e.g. the shipped test script) already
filter out all-zero/no-data pMats via `bins2use` before calling
`lineTraj_decode` at all.

`min(bins2use):max(bins2use)` in `pMat_data.mat READ ME.txt`'s own extraction
recipe (also used by `lineTraj_testScript.m`) can include *interior* NaN/zero
columns between the first and last valid one — not just the contiguous outer
trim it might look like. `line_trajectory.best_line`/`line_residual` already
zero out NaNs before use, matching how the MATLAB script's own next line
(`pMat(isnan(pMat))=0`) handles it, so this is a non-issue, just worth noting
for anyone assuming `bins2use` selects a fully-valid contiguous block.

## 14. `fs_fastGridness.m` -> `fast_gridness`, v0.5 — five confirmed divergences from `multiGridness.m`

Ported `fs_fastGridness.m` into `grid.py::fast_gridness`. It only has one real
caller (`fs_masterAnalysis.m`, itself not yet ported — the FieldShuffle
significance-testing orchestration), but the porting plan groups it with v0.5
as a standalone, independently-validatable analysis function, and its own
header comment describes it as "a cut down version of autocorr props" — i.e.
a second, faster implementation of the *same measurement* already ported as
`gridness`/`grid_properties` (`multiGridness.m`/`sacProps.m`).

**It is not numerically interchangeable with `gridness`/`grid_properties`,
even on the identical SAC.** Reading both MATLAB sources side by side turned
up five confirmed algorithmic divergences (not assumed — each checked by
re-reading both `.m` files):

1. Peak detection: `fs_fastGridness.m` calls `imregionalmax(sac)` directly
   (`sac` already has NaN substituted to -1). `multiGridness.m` first masks
   non-positive bins to -1 (`tmp(sac<=0)=-1`) *before* calling
   `imregionalmax`. So `fs_fastGridness.m` can in principle pick up spurious
   "peaks" among negative-valued or NaN-substituted SAC regions that
   `multiGridness.m` deliberately excludes.
2. Candidate peak list: `fs_fastGridness.m` builds it straight from
   `find(peaksSAC)` -- one candidate per *pixel* of any local-max plateau.
   `multiGridness.m` goes through `bwlabel`+`regionprops` first -- one
   candidate per *connected plateau* (its centroid). A multi-pixel plateau
   is therefore one peak to `multiGridness.m` and several to
   `fs_fastGridness.m`.
3. Every radius/mask comparison in `fs_fastGridness.m` uses a
   distance-from-centre grid **rounded** to the nearest integer bin
   (`dstFrmCentre=round(dstFrmCentre)`, needed for its precomputed
   `allMsks` disc-mask family). `multiGridness.m`'s equivalent
   (`centreRadi`) is left continuous throughout. Boundary pixels of a
   radius-thresholded mask can therefore differ by one ring between the two.
4. Grid scale: `fs_fastGridness.m`'s `gScl` is the **mean** distance to the
   closest ~6 peaks (`mean(dstFrmCent(2:end))`). `grid_properties`'s
   `scale` is the **median** (`sacProps.m`'s `median(...)`).
5. Expanding gridness picks a different central/inner radius entirely.
   `multiGridness.m` (Killian 2012 style, already in `gridness`) uses half
   the mean peak distance (`round(mean(peakDist)/2)`).
   `fs_fastGridness.m` (Langston et al. 2011 style) instead computes a
   *radial mean-correlation profile* (mean SAC value at each integer radius)
   and takes the first regional minimum or first negative crossing of that
   profile as the inner radius, falling back to a hard-coded 5 bins if
   neither exists.

Also confirmed while porting (not a divergence from `multiGridness.m`, but
worth noting): the half-height field-extent sub-step (`sf_findPeakExtent`)
*is* byte-for-byte identical between the two files, so `fast_gridness` reuses
the already-validated `_peak_fields`/`_rotational_gridness` helpers built for
`gridness` for that piece and for the final rotated-correlation step —
genuine shared code, not an assumption of equivalence. One subtlety this
reuse required getting right: `_peak_fields` resolves *overlapping*
half-height fields by keeping the numerically larger field-label id, and
MATLAB assigns those ids as the peaks' *original* (pre-distance-sort)
`find()`-order positions, not a fresh 1..7 relabelling by distance — using
sequential ids instead would silently change which field wins an overlap.

**Validation.** No MATLAB runtime available (same constraint as §13), so no
bit-exact comparison against real MATLAB output. Used: (a) the existing
synthetic hexagonal-lattice fixture (from `test_grid.py`) — `fast_gridness`
scores it high on both standard and expanding gridness, same as `gridness`;
(b) a single-Gaussian-blob non-grid fixture — scores low/NaN on both, same
qualitative call as `gridness`; (c) a real validated grid cell
(`gridCell1838_t2c3`) — `fast_gridness` and `gridness`/`grid_properties`
agree to within ~0.05 gridness and ~0.4 bins of scale, both correctly call it
a strong grid cell. Given the five confirmed divergences above, exact
agreement was never expected; the ballpark cross-check is there to catch a
port bug that would make the two wildly disagree, not to assert equivalence.

## 15. `xSacProps.m` -> `cross_sac_props`, v0.5 — a likely real MATLAB fragility found

Ported `xSacProps.m` into a new `grid_phase.py::cross_sac_props`, plus a
public `grid.py::cross_correlogram` (the already-private `_pearson_all_lags`
already supported two distinct maps — the plan's "confirm the ported SAC
supports cross-correlation mode" note just needed a thin public wrapper).

**Validated the core geometry (cart2pol/wrapTo2Pi/distance-then-angle
sorting/1-based-vs-0-based indexing) is an exact translation** using a
hand-built xSAC with peaks placed at exactly-known pixel positions on a
strictly-varying (no ties, no flat regions) background — every returned
field (`phase_mag`, `phase_ang`, `peak_coord`, `peak_corr`, and all three
`axis_ang_dist` rows) matched the known ground truth exactly
(`test_cross_sac_props_recovers_known_peak_geometry_exactly`).

**Then found something worth flagging per the standing "watch for bugs in
the original MATLAB" instruction — likely a genuine fragility in
`xSacProps.m` itself, not just this port.** Running the same function on a
*realistic* xSAC (two synthetic hex-grid rate maps, deliberately offset from
each other, cross-correlated) gave a clearly wrong answer: `peak_corr` came
back as exactly `0.0` and `phase_mag` didn't match the true injected offset.
Root cause, confirmed by direct inspection of the intermediate arrays (not
assumed): `xSacProps.m` zero-clips every negative correlation bin
(`xsac(xsac<0)=0`) *before* calling `imregionalmax`. For a real correlogram
(58-62% of bins were negative/zero-clipped in these test cases — real
correlation surfaces spend most of their area away from genuine structure)
this creates one or more huge flat regions at exactly `0`, and
`imregionalmax` treats a flat region as trivially "locally maximal"
everywhere within it — so `find()` returns thousands of spuriously-tied,
meaningless zero-valued "peak" candidates alongside the handful of genuine
ones. Since every field this function returns is chosen by proximity
(closest peak to the xSAC's geometric centre, then closest six peaks to
*that* peak), and a zero-plateau boundary pixel is almost always nearer than
the next genuine hexagonal peak, the true structure gets swamped. Manually
separating the candidate peaks by value (`==0` vs `>0`) confirmed the
*genuine* peaks were exactly where a correctly-working algorithm should put
them (matching the known grid scale) — it's specifically the
proximity-search step that gets misled, not the peak values themselves.

This only manifested because a *realistic* (not toy) synthetic case was
tried in addition to the clean ground-truth one — a good example of why
"validated on a clean synthetic case" and "validated on realistic data" are
different tiers (see the validation-discipline checklist). Faithfully
preserved (not fixed) per goal 1: `cross_sac_props` zero-clips exactly as
`xSacProps.m` does. Documented prominently in the function's own docstring
(check `peak_corr` before trusting the other fields) and pinned by
`test_cross_sac_props_swamped_by_zero_clipped_background_on_realistic_input`,
so this is a known, understood limitation rather than a silent trap for
whoever next calls this function on real data.

**Requires >= 7 detected peaks, exactly as unforgivingly as the MATLAB.**
`xSacProps.m` indexes `dist(2:7)` unconditionally, which errors in MATLAB if
fewer than 7 peaks exist. Python's slicing doesn't error the same way (`x[1:7]`
on a shorter array silently truncates), which would have made `cross_sac_props`
silently return a malformed `axis_ang_dist` with fewer than 3 rows instead of
failing — added an explicit `ValueError` for this precondition so the failure
mode matches MATLAB's in spirit (a clear error) even though the literal
mechanism differs.

## 16. `brandonRegEcc.m` -> `brandon_regularise`, v0.5 — a self-acknowledged rough method

Ported `brandonRegEcc.m` into `grid_ellipse.py::brandon_regularise`, alongside
the already-ported least-squares ellipse fit/correction
(`fit_grid_ellipse`/`regularise_elliptic_sac`). The MATLAB source's own
docstring is candid that this is a *worse* method kept only as a fallback:
"a bit clunky," "doesn't actually work that well — certainly not as well as"
the least-squares fit, and it flags its own first axis-length formula as
possibly needing a sqrt of a negative number ("I've just ignored the
negative sign").

**Confirmed by reading the source directly (not assumed from the similar
name) that peak detection here follows `fs_fastGridness.m`/`xSacProps.m`'s
convention** — raw `imregionalmax` with no positive-value restriction — NOT
`multiGridness.m`'s convention (used by `fit_grid_ellipse`). Since this
"unrestricted `imregionalmax` + column-major `find()`" pattern is now needed
by three ported functions (`fast_gridness`, `cross_sac_props`,
`brandon_regularise`), consolidated it into a shared
`grid.py::_raw_regional_max_peaks` helper — a refactor of code written this
session (not a change to older, separately-validated code), verified
behaviour-preserving by re-running the existing `fast_gridness`/
`cross_sac_props` test suites unchanged before and after.

**The MATLAB's acknowledged possibly-negative-sqrt is real and was handled
carefully, including an asymmetry between the two output methods.** MATLAB's
`sqrt` promotes a negative real to a complex result; a naive `numpy.sqrt`
would instead give `nan` with a warning — used `numpy.lib.scimath.sqrt` to
match. Method 1's call wraps its axis lengths in `abs()` (`abs([m1MajR,
m1MinR])`, which for a complex number in MATLAB returns the magnitude) but
**method 2's call does not** (`[m2MajR, m2MinR]`, bare) — the identical
failure mode is treated differently between the two methods, which reads
like an oversight rather than a deliberate choice. Faithfully preserved: for
method 2, a complex result's real part is used before comparison/regularisation
(since `regularise_elliptic_sac` sorts its `ab_scale` argument, undefined for
complex numbers in Python, unlike MATLAB's real-part-based complex
ordering) — this only changes behaviour in an already-degenerate corner the
original author's own comment says isn't well understood, not a
divergence from any well-defined MATLAB behaviour.

**Validation — structural, matching how the *better* least-squares method
was itself validated (no MATLAB runtime, and this method is explicitly
weaker).** No crashes on real or synthetic input; correctly returns the
input SAC unchanged (both methods) when fewer than 4 peaks are found (using
the same strictly-monotonic-ramp construction established in `test_grid.py`
to avoid the flat-background tie confound — see §14). A genuine,
algebraically-derivable check: for a perfectly regular (non-elliptical)
hexagonal lattice, the six central peaks are equidistant from the centre, so
`closeR == distR`; substituting into method 2's minor-axis formula gives
`sin²θ/((1/R²)−(cos²θ/R²)) = R²` exactly, independent of θ — i.e. method 2's
major/minor estimates must come out equal, triggering
`regularise_elliptic_sac`'s own `a==b` shortcut and returning the SAC
completely unchanged. Confirmed exactly on a synthetic regular hex grid. On
a deliberately elongated synthetic hex grid, gridness improves after
correction for at least one of the two methods (method 2 sometimes leaves
the SAC unchanged even here, and method 1 shows the improvement) — consistent
with, not contradicting, the source's own "doesn't always work" caveat.

## 17. `fs_shufFields.m` -> `shuffle_fields`, v0.7 — the field-reallocation shuffle, plus a real MATLAB fragility

Ported `fs_shufFields.m` into `field_shuffle.py::shuffle_fields`: randomly
reallocates every watershed field's centre to a new bin, then works through
each field's own member bins (in order of closeness to its own peak, with a
random per-field rotational offset broken tie) placing each one at the
closest still-available bin to the field's *new* centre (also with an
independent random per-field rotational offset), processing all fields
round-robin (one bin from each in turn) rather than one field fully before
the next, so no single field can systematically monopolise the best
positions.

**Randomness means bit-exact comparison is categorically impossible here**,
unlike everywhere else in this codebase — MATLAB's unseeded `rand`/`randperm`
run a different RNG algorithm than anything numpy can reproduce, regardless
of translation care. Given that, `shuffle_fields` operates throughout in
0-based **row-major** (C-order) flat indices — matching
`field_bin_locations`'s own established convention — rather than replicating
MATLAB's 1-based column-major linear indices, since doing so would only add
translation risk for zero validation benefit (the specific index-numbering
scheme is an implementation detail; what's faithfully replicated is the
algorithm's *logic*). Accepts a `numpy.random.RandomState` or seed for
reproducibility, which MATLAB's own unseeded calls don't offer at all.

**Two sort-key orders needed care and were double-checked against the
source rather than assumed.** `sortrows(cell2mat(fldBinLoc(nn,4:5)),[2,1])`
concatenates *(theta, radius)* as columns 1–2, so column-spec `[2,1]` sorts
by radius (column 2) primary, theta (column 1) secondary — matching the
code's own comment ("based first on dist to peak then orientation").
`np.lexsort` takes its *last* key as primary, so this is
`np.lexsort((theta, radius))`, not the more visually-obvious-looking
`(radius, theta)` — got this backwards on a first draft and caught it by
re-deriving the column order from `cell2mat`'s concatenation order rather
than trusting an initial intuition, before writing any code.

**A genuine MATLAB fragility found: silent output corruption when
`staticBin` over-constrains the problem, faithfully preserved as a raised
error instead.** If a field's ranked candidate-destination list is fully
exhausted (every bin already claimed) before all its member bins are
reallocated — which the greedy, first-available round-robin approach can
hit even when *globally* there are exactly enough free bins overall, if
another field's earlier greedy picks consumed more than its "fair share" —
MATLAB's `find(vector,1)` returns `[]`, and `shufInds([])=x` is a silent
no-op in MATLAB. That means MATLAB would quietly return an **incomplete or
malformed permutation** (some bin's value never placed anywhere, its
target slot left at the preallocated `0`) rather than fail loudly. Confirmed
this is easy to trigger in a genuinely over-constrained synthetic case (two
watershed fields tiling an entire small ratemap with no unlabelled
background slack, plus a `staticBin` scattered evenly enough to leave zero
give). Given the alternative is either a bare, uninformative `IndexError` or
silently returning a corrupted permutation a caller might not notice,
`shuffle_fields` raises an explicit `ValueError` instead — a case where
"faithfully replicate" is served better by a loud failure than by
reproducing MATLAB's specific silent-corruption mechanism (this is the same
judgement call as `cross_sac_props`'s `<7`-peaks `ValueError`, see §15). In
realistic usage (`staticBin` = genuinely unvisited bins, typically a modest
and usually spatially-compact minority of a well-explored arena, not an
adversarially-scattered majority) this should not arise.

**Validation — structural, the only kind possible for a genuinely randomised
algorithm.** Confirmed: the output is always a valid permutation of every
bin (no bin lost or duplicated); static bins always stay exactly in place;
identical results given the same seed (reproducibility); field centres are
relocated to a different bin than their original position; and — the
algorithm's actual stated purpose — bins belonging to a field remain
clustered *around the field's own new centre* after reallocation rather than
scattering across the whole map (checked quantitatively: mean post-shuffle
distance from a field's own bins to its own new centre is well under half
the map's diagonal, on a synthetic two-field case). Also runs end-to-end
without error on the real, watershed-segmented r1838 grid cell.

## 18. `fs_masterAnalysis.m`/`identify_grid_cell.m` -> `field_shuffle_significance`/`is_grid_cell`, v0.7 done — two more genuine MATLAB bugs, plus a real capacity fragility discovered on real data

Ported the last of the FieldShuffle orchestration: `field_shuffle_significance`
(`fs_masterAnalysis.m`) and `is_grid_cell` (`identify_grid_cell.m`), both
redesigned around `Recording` like `crunch_decode` (validate the component
functions individually — `fast_gridness`, `brandon_regularise`,
`shuffle_fields` — this end-to-end orchestrator is not bit-exact, and
inherits `shuffle_fields`' genuine per-call randomness on top). Also added
the small `_brandon_gridness` helper (port of the `sf_getBrandonG`
subfunction embedded in `fs_masterAnalysis.m`: best-of-two-methods
Brandon-corrected gridness, using `np.fmax` rather than `max()`/
`np.maximum` to match MATLAB `max()`'s NaN-ignoring semantics — the same
recurring MATLAB-vs-numpy NaN gotcha hit several times elsewhere in this
project).

**Two further genuine MATLAB bugs found by reading the source directly**
(per the standing "watch for bugs in the original" instruction), neither
reproduced — both make more sense to route around than to faithfully
replicate, for reasons specific to each:

1. **`fs_masterAnalysis.m`'s temporal-shuffle wraparound assigns the wrong
   variable.** `rotSpkPos=mod(spkPos+rotAmnt,trlLgthPosSamp);
   rotSpkPos(rotSpkPos==0)=nLpForShuf;` — MATLAB is 1-based, so a
   mod-wrapped result of exactly `0` is invalid and should remap to the
   *last* valid index, `trlLgthPosSamp` (the trial length). Instead it
   substitutes `nLpForShuf` — the number of shuffle *iterations* (typically
   100), an unrelated quantity — silently placing a small number of rotated
   spikes near the very start of the trial instead of wrapping correctly.
   **Not reproduced, and not really "fixable" as a faithful port either:**
   in 0-based Python indexing this whole workaround is moot, since `0` is
   already a perfectly valid index — `(spk_idx + rotation) % n_pos` is used
   directly with no special case, which sidesteps the bug's relevance
   entirely rather than needing to decide whether to replicate it.
2. **`identify_grid_cell.m` references an undefined variable.** It calls
   `fs_masterAnalysis(data.pos, spkPos, ...)`, but the function's own
   parameter is named `pos`, not `data.pos` — `data` is never defined
   anywhere in the function. This would raise "undefined variable" in
   MATLAB on every single invocation, meaning `identify_grid_cell.m` was
   very likely never actually run as a standalone function call (only
   possibly ever exercised by pasting its body into a script/workspace
   context where a variable happened to be named `data`). There is no
   well-defined numerical behaviour here to preserve — only a guaranteed
   crash — so `is_grid_cell` simply calls `field_shuffle_significance` with
   the actual `recording` argument, as any correct implementation must.

**A significant capacity fragility found and fixed, discovered only because
real (not just synthetic) data was tried in addition to synthetic
fixtures.** Running the new `is_grid_cell` across all four documented
example cells (t2c1, t2c3 = grid; t2c4, t2c5 = non-grid), the two grid
cells worked immediately but t2c4 (non-grid) deterministically raised
`shuffle_fields`'s "ran out of available destination bins" error on
**every one of 30 different seeds tried** — not a rare, unlucky draw.
Root-caused by direct inspection of the field-segmentation numbers (not
assumed): unmarked watershed always tiles the *entire* ratemap into
catchment basins, including never-visited regions, and t2c4 happens to have
a tiny, marginal field whose single highest-rate bin is itself an unvisited
(static) bin. Both MATLAB's `fs_shufFields.m` and a literal translation of
it are vulnerable to this: every field's peak gets relocated to a brand-new
slot from the movable-bin pool regardless of whether that peak coincides
with a bin flagged static (nothing in `fs_binMembLoc.m`/
`field_bin_locations` checks for this), so a static-coincident peak
consumes a slot from supply while contributing nothing to declared demand —
deterministically over-constraining the allocation by exactly one slot.
Confirmed by direct arithmetic on the real data: demand exceeded supply by
exactly 1 for t2c4's segmentation, matching the "always fails" empirical
result exactly. **Fixed in `shuffle_fields`** (not merely documented, since
this is squarely the kind of correctness issue worth resolving rather than
faithfully reproducing MATLAB's silent-corruption fallback): a field whose
peak is itself a static bin is now treated as entirely static (all its bins
stay in place) rather than relocated — sensible on its own terms too, since
such a "field" is watershed noise from a never-visited area, not a genuine
firing field. Confirmed this resolves the concrete case (30/30 seeds
succeed on t2c4 afterward) and added a second, purely defensive
`ValueError` for a distinct (and, after this fix, essentially
unreachable-via-self-consistent-input) failure mode: more genuine fields
than movable bins to give them new centres. Full detail, including the
worked arithmetic, in `field_shuffle.py::shuffle_fields`'s docstring.

This is a good illustration of why real, non-grid, "should fail the test"
cells are worth deliberately including in validation, not just the
documented positive (grid-cell) examples — the fragility here specifically
affects sparse/patchy-firing cells, exactly the population `identify_grid_cell`
exists to correctly classify as non-grid, and would have gone undetected
testing only the two grid cells.

**Validation.** Structural, as for the rest of the FieldShuffle chain (no
MATLAB runtime, plus genuine randomness in the field-shuffle step means
bit-exact comparison is impossible regardless): `field_shuffle_significance`
runs end-to-end on the real example trial and its standard/expanding
gridness for the "no correction" case matches `fast_gridness` called
directly on the same ratemap exactly (a strong internal-consistency check,
since both paths independently build the same ratemap and call the same
gridness function); `_brandon_gridness` matches an independently-computed
best-of-both-methods calculation exactly; `is_grid_cell` runs cleanly and
returns correctly-shaped output across all four documented example cells
(both grid and non-grid) with multiple shuffle counts and seeds, including
the t2c4 case that specifically exercised the capacity-fragility fix above.

## 19. Things confirmed correct and bit-exact

`.pos`/tetrode/`.cut` binary readers (big-endian, 216-byte spikes, 20-byte pos
records); LED swap + speed filtering; NaN interpolation; 20-tap boxcar position
smoothing (scipy `origin=0` matches MATLAB `imfilter`); `timestamp/timebase`;
1-based→0-based `pos_sample`; cluster selection; and `histnd` spatial binning.
Line-trajectory best-fit-line round-trip self-consistency (118/118 real replay
events from `pMat_data.mat`, see §13). Cross-SAC phase geometry
(`cross_sac_props`) against a hand-built xSAC with exactly-known peak
positions (see §15). Brandon-method ellipse regularisation's regular-hexagon
identity (`m2MinR == m2MajR` algebraically — see §16).

## 20. `eeg_instant_freq_power.m` / `detect_ripples.m`, v0.6b complete — a confirmed doc/code default mismatch, and a "bug" that traces out harmless

Ported the last two theta/LFP functions: `eeg_instant_freq_power.m` ->
`instant_freq_power()` (Blackman-windowed FIR bandpass + `filtfilt` +
Hilbert transform -> instantaneous phase/amplitude/frequency, matching
`scipy.signal.firwin(..., window="blackman", fs=eeg_freq)`'s default
Nyquist-normalisation and passband scaling against MATLAB's
`fir1(eegFreq, Wn/nyquist, blackman(eegFreq+1))`) and `detect_ripples.m` ->
`detect_ripples()` (Csicsvari-style ripple detection built on top of it).
Both live in `analysis/lfp.py` alongside the already-ported `power_spectrum`.

**The MATLAB's `interp2`-based instantaneous-frequency re-centring was
traced through by hand and replaced with an equivalent direct computation,
not transcribed literally.** `diff(phase)` gives a length-(N-1) array of
"between-sample" frequency estimates; the source re-centres this onto all N
original samples via `interp2([freqEEG, ones(...)], ones(N+1,1),
[0.5:1:N+0.5]')` — a two-argument `interp2` call that (traced step by step:
implicit grid X=1:2/Y=1:(N-1), query xi always column 1, yi at half-integer
positions 0.5..N+0.5) is exactly "linearly interpolate the diff array at its
own midpoints, with the two out-of-domain endpoint queries (0.5 and N+0.5)
returning NaN". Reproduced directly: `freq[0]=freq[-1]=nan`,
`freq[i]=mean(d[i-1],d[i])` for interior `i`. Confirmed numerically
equivalent (not just theoretically) via the synthetic 8Hz-carrier test below.

**A confirmed, real doc/code mismatch: `detect_ripples.m`'s default
`ripThresh` is 5, not the 4 its own docstring claims.** The function header
comment says "The detection threshold is defined as mean + std * ripThresh.
If not supplied default is 4", explicitly citing "Csisvari uses 4" — but the
actual code sets `dfRipThresh = 5` with its own separate inline comment
("Freyja uses 3 here - Csisvari uses 4"), i.e. the chosen default is
*stricter* than either cited literature value and contradicts the
function's own documented behaviour. This is a real inconsistency in the
original MATLAB (not a porting ambiguity), confirmed by reading both lines
side by side. Per the standing policy of matching actual behaviour over
stale documentation, `detect_ripples()`'s default is 5, matching the code
that actually runs — prominently flagged in the Python docstring so a user
who has read the MATLAB header isn't misled. Pinned by
`test_detect_ripples_default_threshold_matches_documented_code_default`.

**`deal(zeros(length(ripStr)))` does build an n×n matrix instead of the
evidently-intended n×1 column — traced through by hand and confirmed
numerically harmless, contrary to the initial `PORTING_PLAN.md` note calling
it a bug to fix.** `zeros(n)` with a scalar `n` is MATLAB shorthand for an
n-by-n matrix, so `valid`/`pkPower`/`meanFreq` all come back n² in size, not
the clearly-intended n×1. But the following loop only ever writes
`valid(k)`/`pkPower(k)`/`meanFreq(k)` for `k=1..n` — and for an n×n matrix
stored column-major, linear indices `1..n` fall entirely within column 1,
i.e. exactly the same memory positions a correctly-shaped n×1 vector would
occupy. Every later read (`valid(n)` inside the loop, then
`logical(valid)` and `ripStr(valid)`/`pkPower(valid)`/`meanFreq(valid)` via
logical indexing, which MATLAB resolves via `find()` over the full
flattened array) only ever touches those same first-`n` linear positions,
since columns 2..n of the oversized matrix are never written and so stay
`false`/`0`. Traced this out by hand before writing any test, precisely
because the porting plan had already flagged it as a bug needing a fix —
confirming it *first* avoided baking an unnecessary behavioural change into
the port. The only real cost in the original is a wasted, potentially large
`n²` memory allocation, not wrong values; `detect_ripples()` simply uses
correctly-shaped 1-D arrays throughout (there's no MATLAB output to
preserve here, since the outputs are provably bit-identical).

**The trailing-ripple drop-off is real and is faithfully preserved.**
`ripStr = ripStr(1:length(ripEnd))` keeps only the first `length(ripEnd)`
starts; if the trace ends while still above the low threshold (one more
rise than fall), that final, unterminated ripple is silently dropped rather
than reported with the trace's end as its boundary. This is well-defined
(not a crash), so it's kept as-is per goal 1, documented prominently in
`detect_ripples()`'s docstring.

**The empty-candidate-list crash is real, realistic, and *is* worked
around.** `ripStr(1)>ripEnd(1)` indexes element 1 of what can be a
genuinely empty array whenever nothing anywhere in the trace crosses the
low threshold — not a contrived edge case, but the expected outcome for any
short or genuinely ripple-free recording (confirmed directly:
`test_detect_ripples_no_activity_returns_empty` hits exactly this path on
pure background noise). MATLAB errors here; "no ripples found" is a
well-defined, useful result, so `detect_ripples()` returns empty arrays
instead, per the standing crash-avoidance policy.

**Also confirmed, cosmetic only:** `ripMinDur=60/1000`'s own inline comment
("i.e. 40ms = 0.004") is wrong twice over — 60ms is not 40ms, and 40ms is
0.04s, not 0.004s. The *value* (0.06s = 60ms) is what's actually used and
what's ported; the comment is simply stale/miscopied, no code consequence.

**Validation.** `instant_freq_power`: a synthetic 8Hz sinusoid (theta band,
60s @ 250Hz) recovers `mean(freq)=8.00Hz` to <0.001Hz away from every-sample
noise floor, and the pos-rate-downsampled version agrees to the same
precision; the first/last samples are confirmed `nan` as expected from the
interp2-equivalent centring; an explicit test confirms the deliberately
non-integer downsampling ratio raises `ValueError`, and another confirms a
trace with a few trailing samples beyond an exact block multiple is handled
by dropping them rather than crashing. On real r1838 LFP (channel 1, `.egf`
@ 4800Hz, first 60s): `nanmean(freq)=9.11Hz`, squarely in the theta band and
closely matching `power_spectrum`'s independently-computed FFT peaks for
the same channel (9.14-9.58Hz, per §10/§9) — strong cross-consistency
between the two very differently-implemented theta measures (FFT peak vs.
Hilbert instantaneous frequency).

`detect_ripples`: three synthetic 150-250Hz oscillatory bursts (200Hz, 150ms
each) injected into low-amplitude background noise (1200Hz trace, well
above the 500Hz-plus Nyquist requirement the MATLAB source itself
documents) are all three recovered, with start times within ~100ms of the
true injected onset (the low-threshold expansion pads the window either
side, as intended by the algorithm's own design) and `mean_freq` within
15Hz of the true 200Hz carrier. A pure-noise trace with no injected bursts
returns zero ripples (exercising the empty-candidate-list fix). On real
r1838 egf data (channel 1, first 120s): 23 candidate ripple-band events
with durations 65-190ms (comfortably above the 60ms minimum) and mean
frequencies clustering 194-215Hz — squarely inside, not at the edges of,
the requested 150-250Hz band, a meaningful sanity check that the filter
design and edge-indexing are correct (a bug here would much more likely
produce frequencies piled up at the band edges or wildly outside it than a
clean mid-band distribution).

**A design decision worth flagging:** `instant_freq_power` accepts
`pos_sample_rate=None` to skip the position-rate downsampling step
entirely. `detect_ripples` always calls it this way, since ripple detection
never needs the downsampled output and the egf-rate signal's length isn't
guaranteed to divide evenly by 50 (unlike the always-integer-ratio case for
standard 250Hz `.eeg`) — this sidesteps a spurious `ValueError` for a
computation that was never going to use that output anyway, rather than
computing and discarding it.
