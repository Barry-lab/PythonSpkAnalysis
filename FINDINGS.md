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
and `detect_ripples.m` deferred, not ported.** Both are more involved
(FIR bandpass design + `filtfilt` + Hilbert transform; ripple detection has
several documented bugs per `PORTING_PLAN.md`) and time did not allow the
same careful read-source/port/validate cycle used elsewhere. Still tracked.

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

## 14. Things confirmed correct and bit-exact

`.pos`/tetrode/`.cut` binary readers (big-endian, 216-byte spikes, 20-byte pos
records); LED swap + speed filtering; NaN interpolation; 20-tap boxcar position
smoothing (scipy `origin=0` matches MATLAB `imfilter`); `timestamp/timebase`;
1-based→0-based `pos_sample`; cluster selection; and `histnd` spatial binning.
Line-trajectory best-fit-line round-trip self-consistency (118/118 real replay
events from `pMat_data.mat`, see §13).
