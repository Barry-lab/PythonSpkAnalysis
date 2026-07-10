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

- `barrelab` implements both via `load_dacq(..., head_pos_weighting=...)`.
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

## 3. Grid-cell analysis port (SAC + gridness), v0.2

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

## 4. Things confirmed correct and bit-exact

`.pos`/tetrode/`.cut` binary readers (big-endian, 216-byte spikes, 20-byte pos
records); LED swap + speed filtering; NaN interpolation; 20-tap boxcar position
smoothing (scipy `origin=0` matches MATLAB `imfilter`); `timestamp/timebase`;
1-based→0-based `pos_sample`; cluster selection; and `histnd` spatial binning.
