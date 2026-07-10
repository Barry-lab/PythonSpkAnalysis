# barrylab — Master Porting Plan

Porting the Barry lab MATLAB toolbox (`Universal_Matlab`) to the Python package **barrylab**.
Architecture: `io/` loaders fill a `core/Recording` → `analysis/` pure functions consume a Recording or
NumPy arrays and never touch disk.

Current package layout (for module-home mapping below):
`src/barrylab/io/axona.py`, `core/recording.py`, `analysis/{maps,spatial,grid,direction}.py`.
New analysis modules proposed here (`place_fields`, `pxd`, `decoding`, `line_trajectory`,
`grid_ellipse`, `spatial_fft`, `lfp`, `field_shuffle`, `_stats`) are additive.

Already ported (not re-planned): `load_dacq`, `Recording`/`Electrode`, `occupancy_map`, `spike_map`,
`rate_map`, `smooth_ratemap` (boxcar/gaussian), `spatial_information`, `adaptive_binning`,
`spatial_autocorrelogram` (xPearson), `gridness`, `grid_properties`, `grid_score`.

---

## 1. Overview

What remains is roughly five families of work: (a) **remaining spatial statistics** — place-field
metrics, spatial coherence/correlation, PPM correction; (b) **head-direction & place×direction (pxd)**
analysis; (c) **population decoding & replay** (Bayesian decoder + line-trajectory fit); (d) **grid
extras** — elliptical correction, 2D-FFT grid components, cross-correlogram phase, expanding gridness;
and (e) **LFP/theta** — instantaneous phase/power, power spectra, ripple detection, phase precession.

The dominating dependency is that **families (a)–(d) need no new infrastructure** — they are pure array
math over positions, directions, ratemaps and autocorrelograms already available on the Recording — so
they can be ported and golden-fixture-validated immediately. **Family (e) is gated on one new-infra
item: EEG/EGF loading into the DACQ loader plus an `lfp` field on Recording.** That single unlock
(`read_eeg_file`, low effort) enables every theta/oscillation/ripple analysis. One further blocker,
**watershed segmentation** (no scipy/skimage equivalent available), gates only the FieldShuffle
sub-project and can be deferred to the end.

Recommended order therefore front-loads the self-contained, deterministically-validatable analysis
(shared utilities → place-field stats → head-direction/pxd → decoding → grid extras), then adds the
EEG/LFP infra and its theta analyses, and finally the watershed-gated FieldShuffle and the optional
extra loaders/viz layer.

---

## 2. Recommended porting order

Phases are coherent, independently shippable batches. Infra-free, fixture-validatable work first;
infra-dependent work strictly after its prerequisite.

**v0.2 — Shared utilities + place-field & stability stats.**
Foundational utilities (`nan_corr`, `histnd`, direction/speed binning) plus the headline spatial
metrics (`pcFieldStats`, `spatialCoherance`, `spatialCorrelation`, `get_true_PPM_boundingBox`). No
infra; all deterministic and fixture-able; `spatialCorrelation` (remapping/stability) is a headline
analysis. Ships the small shared bits everything else reuses.

**v0.3 — Head-direction & place×direction. STATUS: DONE.**
`circular_mean`, `directional_kldiv` (+ KL via scipy), `directional_fft`, `led2v1_dir_offset`, and the
`pxd`/`histnd`/`wrapper_pxd` unit are all shipped in `analysis/direction.py` and `analysis/pxd.py`.

**v0.4 — Population decoding, replay & circular–linear regression. STATUS: decoding core + phase precession DONE; lineTraj_* deferred.**
Bayesian decoder trio (`decode_calcPost`, `decode_getSpkCnt`, `decode_processPost`),
`decode_calcBayesPost`, the `lineTraj_*` group (strong shipped oracle `pMat_data.mat`), the
Recording-level `crunch_decode` orchestrator, and `circ_regress_cb` (phase-precession circular–linear
regression — infra-free, built-in Fisher B.20 oracle). No infra; grouped because it is the biggest
concentration of column-major/1-based index risk and shares validation discipline.

**v0.5 — Grid extras: elliptical correction, 2D-FFT, cross-SAC, expanding gridness.**
`gridEllipse_fit`, `regularise_elliptic_sac` (single copy), `xSacProps`, the 2D-FFT subpackage
(`fftVars`, `get_2d_fft`, `get_polar_fft`, `main_spatial_fft2`), and `fs_fastGridness`
(expanding/Langston gridness). Niche siblings (`brandonRegEcc`, `sacProps90deg`) only if elliptical /
square-grid analysis is wanted. No infra, **but** confirm the ported SAC supports cross-correlation
mode (needed by `xSacProps`).

**v0.6 — NEW INFRA: EEG/EGF loading → then LFP/theta analysis.**
*v0.6a (infra):* add `.eeg` (int8, 250 Hz) and `.egf` (int16, 4.8 kHz) reading to the DACQ loader, an
`lfp` field on Recording, and `eeg_convert_to_volts`. *v0.6b (analysis, after 6a):*
`eeg_instant_freq_power`, `eeg_power_spectra`, `detect_ripples` (ripples require `.egf`). This is the
only phase blocked on new infra; everything above is independent of it.

**v0.7 — Watershed-gated FieldShuffle (resolve watershed first).**
`fs_fixWatershed`, `fs_binMembLoc`, `fs_shufFields`, `fs_masterAnalysis`, `identify_grid_cell`. Blocked
on a watershed segmentation decision (vendor one / add skimage / replace with peak+threshold). Last
because it is a self-contained research sub-project with the heaviest single blocker.

**Optional / parallel — extra loaders & viz.**
`read_INP_file`, `read_stm_file`, `read_KwikRPi`, `write_cut_file`, and a `viz/` layer
(`image_rm`, `tintColorMap`). Independent of the above; port on demand when the data types / export /
plotting are actually in scope.

---

## 3. Per-module tables

Consumes = **R** (Recording) or **A** (arrays). Oracle: **GF** = golden fixture via
`generate_fixtures.m` on shipped Example_data feasible; **SYN** = synthetic/analytic known-result;
**STRUCT** = structural only; **MAT** = shipped `.mat` oracle.

### Analysis_Dir — head direction
| function | barrylab target | consumes | needs-infra | oracle | priority |
|---|---|---|---|---|---|
| `cb_circMean` | `analysis/direction.py::circular_mean` | A | no | GF+SYN | core |
| `kldiv` | use `scipy.stats.entropy(base=2)` (don't port verbatim; keep sym/js only if needed) | A | no | SYN | skip/superseded |
| `dir_kldiv` | `direction.py::directional_kldiv` | A | no | GF+SYN | useful |
| `circ_FTT` | `direction.py::directional_fft` | A | no | GF+SYN | useful |
| `led2v1_dir_offset` | `direction.py::led2v1_dir_offset` | A | no | GF+SYN | niche |
| *(upstream gap)* `dir_rate_map` builder | `direction.py::dir_rate_map` | R | no | GF | core (enables the two above) |

Watch-points: `circular_mean` returns `[0,2π)` vs CircStats `(-π,π]` — `led2v1` needs the re-wrap.
`circ_FTT` 1-based→0-based (`power(7)`→`power[6]`; use `relFreq` directly, not `relFreq+1`). KL is
base-2 (bits). `imfilter 'symmetric'`→ndimage `reflect`, `'circular'`→`wrap`.

### Analysis_Freq — theta / LFP
| function | barrylab target | consumes | needs-infra | oracle | priority |
|---|---|---|---|---|---|
| `getegf` (+ `.egf` reader) | `io/axona.py::_read_egf` | path | **YES (lfp)** | GF | core (infra) |
| `eeg_convert_to_volts` | `io/axona.py::_eeg_to_volts` | A+meta | **YES (lfp)** | SYN+GF | core (infra) |
| `eeg_instant_freq_power` | `analysis/lfp.py::instant_freq_power` | A | **YES (lfp)** | GF+SYN | core (after infra) |
| `eeg_power_spectra` | `lfp.py::power_spectrum` | A | **YES (lfp)** | GF+SYN | core/useful (after infra) |
| `detect_ripples` | `lfp.py::detect_ripples` | A | **YES (egf)** | SYN+GF (fix bugs) | useful (after infra) |
| `circ_regress_cb` | `analysis/phase_precession.py::circ_lin_regress` | A | **no** | SYN (Fisher B.20 R≈−0.316) | core (portable now) |

Bugs to fix on port: `getegf` decimate factor vs comment (code wins: /4, /2 default IIR);
`detect_ripples` `deal(zeros(n))` makes n×n matrix, empty-`ripStr` crash, trailing-ripple drop, default
`ripThresh=5`; `eeg_power_spectra` float-equality argmax → `np.argmax`; reshape-to-50 Hz assumes exact
multiples (parametrise on `pos_sample_rate`).

### Analysis_Spatial_Field — place-field / spatial stats
| function | barrylab target | consumes | needs-infra | oracle | priority |
|---|---|---|---|---|---|
| `pcFieldStats` | `analysis/place_fields.py::field_stats` | A | no | GF+SYN | useful |
| `spatialCoherance` | `place_fields.py::spatial_coherence` | A | no | GF (matches Tint) | useful |
| `spatialCorrelation` | `place_fields.py::spatial_correlation` | A | no | GF (same-size case exact) | core |
| `histnd` | `analysis/_binning.py::histnd` | A | no | GF | core (shared) |
| `pxd` | `analysis/pxd.py::pxd` | A | no | GF (float-tol) | useful |
| `wrapper_pxd` | `pxd.py::pxd_ratemap` (strip plotting) | R | no | GF | useful |
| `skaggsInfo`/`skaggsAdaptiveBinning` | — | — | — | — | **DONE** (spatial_information / adaptive_binning) |

Bugs: `wrapper_pxd` `dataRng` uses column-2 max twice for X (should be column-1) — mis-scales X;
`pcFieldStats` no-firing check `nansum(rmS)` on a matrix (intent was `rmS(:)`); `pxd` denominator floor
substitutes `tol` inside the log term only — preserve. `spatialCorrelation` per-cell branch returns a
vector (docstring says mean — averaging line is commented out); preserve, expose `reduce=` option.
`bwlabel(_,4)`→`scipy.ndimage.label` with 4-conn structure.

### Analysis_Grid — remaining
| function | barrylab target | consumes | needs-infra | oracle | priority |
|---|---|---|---|---|---|
| `gridEllipse_fit` | `analysis/grid_ellipse.py::fit_grid_ellipse` | A | no | GF | useful |
| `gridEllipse_correct` | `grid_ellipse.py::regularise_elliptic_sac` | A | no | GF | useful |
| `regularise_eliptic_grid` | *(byte-identical dup)* alias to above | A | no | — | **skip (dup)** |
| `brandonRegEcc` | `grid_ellipse.py::brandon_regularise` | A | no | GF (low value) | niche |
| `sacProps90deg` | `grid.py::square_grid_props` (or `symmetry=4` flag) | A | no | GF+STRUCT | niche |
| `xSacProps` | `analysis/grid_phase.py::cross_sac_props` | A | no (needs SAC cross mode) | GF | useful |
| `fftVars` | `analysis/spatial_fft.py::FFTParams` (dataclass) | — | no | — | useful |
| `get_2d_fft` | `spatial_fft.py::power_spectrum_2d` | A | no | GF | useful (subpkg core) |
| `get_polar_fft` | `spatial_fft.py::to_polar_fft` | A | no | GF | useful (subpkg core) |
| `main_spatial_fft2` | `spatial_fft.py::grid_fft_components` | R | no | GF (deterministic parts) + seeded shuffle | useful |
| `main_spatial_fft` (v1) | — | — | — | — | **skip (superseded by v2)** |
| `fs_fastGridness` | `grid.py::expanding_gridness` | A | no | GF | useful (if expanding gridness not already in grid_score) |
| `fs_fixWatershed` | `analysis/field_shuffle.py::fix_watershed` | A | **watershed** | STRUCT | useful (gated) |
| `fs_binMembLoc` | `field_shuffle.py::bin_member_loc` | A | no | GF | useful (gated w/ shuffle) |
| `fs_shufFields` | `field_shuffle.py::shuffle_fields` | A | **watershed** | seeded | useful (the novel algo) |
| `fs_masterAnalysis` | `field_shuffle.py::field_shuffle_gridness` | R | **watershed** | seeded + real-cell GF | useful (gated) |
| `identify_grid_cell` | `field_shuffle.py::is_grid_cell` | R | **watershed** | STRUCT | niche |
| `fs_findFields`, `main`, `wraper_scale_range`, `genSynData_*`, `interpRes_*`, `testRealData_*` | — | — | — | — | **skip** (drivers/scripts/superseded — see §5) |

Cross-cutting MATLAB-isms (whole module): `imregionalmax`→`ndimage.maximum_filter(x,3)==x` + plateau
collapse via `label`+`center_of_mass`; `imrotate` is CCW-positive vs scipy `rotate` CW-positive
(**negate the angle**, verify on fixture); centre `ceil(size/2)` (1-based) → `size//2` (0-based);
`cart2pol` with y inverted; `bwperim`→erosion XOR. **skimage unavailable** → all covered by
`scipy.ndimage` except watershed.

### Decoding
| function | barrylab target | consumes | needs-infra | oracle | priority |
|---|---|---|---|---|---|
| `decode_calcPost` | `analysis/decoding.py::bayesian_posterior` | A | no | GF (bit-exact) | core |
| `decode_getSpkCnt` | `decoding.py::spike_counts_per_window` | A | no | GF | core |
| `decode_processPost` | `decoding.py::decode_ml` | A (+bin index) | no | GF (avoid ties) | core |
| `decode_calcBayesPost` | `decoding.py::bayesian_posterior_with_prior` | A | no | GF (unique maxima) | useful |
| `decode_crunchData` | `decoding.py::crunch_decode` (rewrite around Recording) | R | no | end-to-end approx | useful |
| `lineTraj_get_res` | `line_trajectory.py::line_residual` | A | no | MAT | useful |
| `lineTraj_define_line` | `line_trajectory.py::_define_line` | A | no | STRUCT (via decode) | useful |
| `lineTraj_reget_lnBnd` | `line_trajectory.py::_reget_line_band` | A | no | MAT roundtrip | useful |
| `lineTraj_decode` | `line_trajectory.py::best_line` | A | no | MAT (`pMat_data.mat`) | useful |
| `lineTraj_testScript` | → convert to pytest roundtrip | — | — | — | **skip (script)** |
| `*.m~` backups | — | — | — | — | **skip (editor backups)** |

Biggest risk (whole module): column-major flatten + 1-based `ind2sub`/`sub2ind` in every ML-decode and
line-mask fn — standardise on `order='F'` + `np.unravel_index`. Also `filter2`/`conv2` kernel anchoring,
MATLAB round-half-away vs numpy round-half-even, and a latent overlap inconsistency between
`decode_getSpkCnt` (supports overlap) and `decode_processPost` (recomputes width from `winLim`,
reshapes pos contiguously → mis-assigns true position for overlapping windows).

### Read_Bin_Data + Analysis_General — loaders & misc
| function | barrylab target | consumes | needs-infra | oracle | priority |
|---|---|---|---|---|---|
| `read_eeg_file` | `io/axona.py::_read_eeg` | path | **YES (lfp field)** | GF | core (unblocks theta) |
| *EGF reader* (no `.m` here) | `io/axona.py::_read_egf` | path | **YES (lfp field)** | GF/STRUCT | useful (write fresh) |
| `read_INP_file` | `io/axona.py::read_inp` | path | optional (events) | GF if `.inp` present | useful |
| `read_stm_file` | `io/axona.py::read_stm` | path | optional (events) | GF if `.stm` present | niche |
| `read_KwikRPi` | `io/kwik.py::load_kwik` (separate loader) | path | **YES (lfp)** + `postprocess_pos_data` | STRUCT (needs Kwik data) | useful/niche |
| `write_cut_file` | `io/axona.py::write_cut` | A | no | roundtrip vs Tint | useful |
| `bin_data` (dir/speed/pxd branches) | `analysis/_binning.py::bin_by_direction/bin_by_speed` (+pxd path) | A | no | GF | core (dir/speed); useful (pxd) |
| `bin_data` (position branch) | — | — | — | — | **DONE** (occupancy_map/spike_map) |
| `cb_bin_data` | — deprecation shim | — | — | — | **skip** |
| `histnd` | `analysis/_binning.py::histnd` | A | no | GF | core (shared; may already be implicit in map ports — check) |
| `pxd` | `analysis/pxd.py::pxd` | A | no | GF (float-tol) | useful |
| `easy_make_polarplot` | `direction.py::directional_ratemap` (sugar) | R | no | GF | niche/skip |
| `image_rm` | `viz/plots.py::image_rm` | A | no | — | niche (viz) |
| `tintColorMap` | `viz/colormaps.py::TINT_CMAP` (constant) | — | no | — | niche (viz) |
| `get_true_PPM_boundingBox` | `spatial.py::true_ppm_bounding_box` | A | no | GF | useful |
| `nan_corr` | `analysis/_stats.py::nan_corr` | A | no | SYN+GF | core (shared util) |

Bugs to fix: `write_cut_file` hard-coded `Exact_cut_for: 1500728a` literal (ignore `fileStem`);
`bin_data` pxd `/50` should use `sample_rate` + `nargoout` typo; `read_KwikRPi` EEG-channel `while true`
loop logic + `strcmp('eeg',varargin)` (string-vs-cell) option detection.

---

## 4. New infrastructure needed

1. **EEG/EGF loading + `Recording.lfp` field (the one big unlock).** Add to the DACQ loader:
   - `.eeg` reader — signed **int8**, ~250 Hz (`read_eeg_file`; low effort, port first).
   - `.egf` reader — **int16** little-endian, ~4.8 kHz (`getegf` model; **no `.m` source in
     Read_Bin_Data** — write from format spec using `read_eeg_file` as template; needed for ripples,
     whose 150–250 Hz band the 250 Hz `.eeg` cannot represent).
   - `eeg_convert_to_volts` — counts→volts using `ADC_fullscale_mv` + per-channel `gain_ch_N` (note the
     `gain_ch_{eegCh-1}` 0-indexed key; keep the −1). Apply at load so `lfp` stores volts.
   - `Recording.lfp` field holding at least `{samples (1-D, clean), sample_rate, channel, gain,
     adc_fullscale_mv}` per channel. Store samples 1-D so analysis fns skip the recurring
     `[value timestamp]` two-column auto-transpose. Parametrise all 50 Hz-hardcoded pos coupling
     (`eegFreq/50`, 50 Hz downsample) on `Recording.pos_sample_rate`.
   - **Unblocks:** `instant_freq_power`, `power_spectrum`, `detect_ripples`, and the LFP path of the
     Kwik loader. `circ_regress_cb` does **not** need it (portable now).

2. **Watershed segmentation (gates FieldShuffle only).** scipy has **no watershed** and **scikit-image
   is not available**. Decide one of: (a) vendor a small marker-controlled watershed seeded from
   `scipy.ndimage` local maxima; (b) add scikit-image as a dependency (weigh against the project's
   no-skimage stance); (c) replace field segmentation with a peak+threshold+Voronoi scheme. Must be
   resolved before `fs_shufFields`/`fs_masterAnalysis`. Nothing outside FieldShuffle needs it.

3. **SAC cross-correlation mode (small).** `xSacProps` consumes a *cross*-correlogram of two ratemaps.
   Confirm the ported `spatial_autocorrelogram`/`xPearson` supports two-input cross mode; if only
   auto mode was ported, add cross mode first.

4. **Dependencies to add (all already implied by scipy/numpy stack):** `scipy.special`
   (`gamma`/`gammaln` for decoding and pxd), `scipy.optimize.minimize_scalar` (circ regression),
   `h5py` (only if the Kwik loader is in scope). No scikit-image (unavailable). No other new deps.

5. **Small shared modules to create early:** `analysis/_stats.py::nan_corr`,
   `analysis/_binning.py::{histnd, bin_by_direction, bin_by_speed}`, and a `dir_rate_map` builder —
   these back several downstream functions and should land in v0.2.

---

## 5. Drop / skip list

- **`kldiv.m`** — replace with `scipy.stats.entropy(p1,p2,base=2)`; keep `'sym'`/`'js'` only if a caller
  needs them (only `dir_kldiv` uses the default mode). Don't port verbatim.
- **`regularise_eliptic_grid.m`** — byte-for-byte identical to `gridEllipse_correct.m`; port once, alias.
- **`main_spatial_fft.m`** — superseded by `main_spatial_fft2` (v2 thresholds in polar space; v1 uses
  obsolete `bin_pos_data`/`goPar`). Author's own header says use v2.
- **`cb_bin_data.m`** — deprecation shim; call `bin_data` directly.
- **`easy_make_polarplot.m`** — pure glue over `bin_data` + `smooth_ratemap`; port only as optional
  sugar, otherwise skip.
- **FieldShuffle drivers/scripts** — `fs_findFields` (exploratory driver, obsolete `autoCorrProps`),
  `main`, `wraper_scale_range`, `genSynData_*` (synthetic-data generators; `genSynData_spkTrain` is
  niche-useful only as a test-data generator), `interpRes_*` (plotting), `testRealData_masterLoop`
  (hardcoded-path population driver), **`testRealData_readDACQ`** (cut-down DACQ loader — entirely
  superseded by ported `load_dacq`). `.mat`/`.zip`/RatPaths are data, not code (but usable as fixtures).
- **`lineTraj_testScript.m`** — convert to a pytest roundtrip against `pMat_data.mat`; don't port.
- **`*.m~`** (`decoded_getSpkCnt.m~`, `decoded_getSpkCntPerWindow.m~`) — editor backups; ignore.
- **`image_rm.m` / `tintColorMap.m`** — plotting; defer to a `viz/` layer, not analysis.

---

## 6. Validation strategy (goal 1: functional parity)

Primary method: **golden-master fixtures via `validation/generate_fixtures.m`** (already present) on the
shipped `Example_data`, compared in pytest against the Python port. Complement with **synthetic/analytic
oracles** where a closed-form answer exists, and **structural checks** where output is stochastic or
data is unavailable.

**Bit-/tight-tolerance validatable (deterministic array→array):** `circular_mean`, `directional_fft`,
`directional_kldiv`, `pcFieldStats`, `spatialCoherance`, `spatialCorrelation` (validate the
correlation/masking logic on **same-size** inputs where interp2 is identity; resize path gets looser
tolerance because `interp2` node placement is hard to match exactly), `histnd`, `nan_corr`,
`get_true_PPM_boundingBox` (subject to matching fill-holes/bbox edge conventions), `get_2d_fft`,
`get_polar_fft`, `gridEllipse_fit`, the whole decoding trio (`decode_calcPost` bit-exact;
`decode_processPost`/`calcBayesPost` on inputs with **unique** posterior maxima to avoid tie-break
divergence), and the `lineTraj_*` group against shipped **`pMat_data.mat`** (prime target — compare
`bstRes`/`bstSpd`/`bstYInterc` per event).

**Float-convergence tolerance:** `pxd` (iterative ML — pin `max_iter`/`accuracy`/`tol` and compare with
tolerance).

**Built-in analytic oracles:** `circ_regress_cb` (Fisher dataset B.20, R≈−0.316, p≈0.057);
`directional_fft` (`cos(6θ)` → power peak at bin 6); `spatialCoherance` (matches Tint exactly).

**Seeded-RNG / structural only:** the shuffle-based analyses — `main_spatial_fft2` (validate the
deterministic 2D/polar transform, then the shuffle threshold with a fixed seed or compare pre-shuffle
spectra only), `fs_shufFields`/`fs_masterAnalysis` (seed RNG; validate the deterministic real-cell
gridness against a fixture, shuffles structurally), `detect_ripples` (thresholds are data-relative —
regression fixture on `.egf` plus a synthetic inject-ripples-into-noise test), `decode_crunchData`
(end-to-end approximate; validate its components individually instead).

**Infra fixtures:** `read_eeg_file`/`eeg_convert_to_volts` — bytewise/known-result unit tests on shipped
`.eeg` + hand-picked counts. Kwik/INP/STM need their own sample files (likely absent) → structural only.

**Discipline note:** for every function with a known MATLAB bug (§3), generate the fixture from the
**corrected** MATLAB call and document the divergence, or keep both a bug-replicating and a corrected
test. Watch the recurring failure sources: 1-based↔0-based indexing, column-major (`order='F'`) flatten,
`imfilter` boundary modes (`symmetric`→`reflect`, `circular`→`wrap`), `imrotate` sign, `filter2`/`conv2`
anchoring, and MATLAB round-half-away vs numpy round-half-even.

---

## 7. Toward goal 2 (rationalise / optimise — AFTER parity)

Redundancies and refactor opportunities spotted across modules, to revisit **only after** functional
parity is proven, and only in ways that preserve identical results:

- **Peak-finding on SACs is duplicated** across `gridEllipse_fit`, `brandonRegEcc`, `xSacProps`,
  `sacProps90deg`, `fs_fastGridness` and the already-ported `grid_properties`/`gridness` — every one
  re-implements "`imregionalmax` → label/centroid → recentre → sort by distance → take N closest". Factor
  one `find_sac_peaks(sac, n_closest, exclude_centre)` helper with a single plateau-collapse policy.
- **SAC/ratemap preprocessing** (`real`, `nan→−1`, `≤0→−1/0`, centre = `ceil(size/2)`) recurs verbatim
  in most grid functions — extract a shared `prep_sac()`.
- **Poisson likelihood core** is written twice in decoding (`decode_calcPost` with the `/obs!` term,
  `decode_calcBayesPost` without) — unify behind one function with a `use_factorial` flag (the specs
  already suggest this).
- **`pxd` appears in two places** (Analysis_Spatial_Field and Read_Bin_Data) — port once to
  `analysis/pxd.py`; `histnd` likewise is duplicated (Read_Bin_Data + the pxd wrapper) and may already
  be implicit inside the ported occupancy/spike maps — consolidate to one `_binning.histnd`.
- **Edge-corrected smoothing** (smooth signal ÷ smoothed validity mask) is re-implemented in
  `get_2d_fft`, `get_polar_fft`, `decode_getSpkCnt`, `eeg_power_spectra`, `wrapper_pxd` — factor a
  single `edge_normalised_filter()`; consider routing `wrapper_pxd`'s positional smoothing through the
  already-ported `smooth_ratemap` (author flagged its own smoothing as incorrect).
- **`nan_corr`-style NaN-safe Pearson** likely already exists inside `spatial_autocorrelogram`/`gridness`
  (xPearson) — dedupe against the new shared `_stats.nan_corr`.
- **Expanding/Langston gridness** in `fs_fastGridness` overlaps the ported `gridness`/`grid_score` — if
  not already included, add it there as an option rather than a parallel implementation; likewise
  `sacProps90deg` is better as a `symmetry=4` flag on the existing gridness than a separate file.
- **Two-column `[value timestamp]` LFP auto-transpose** recurs in three Freq functions — eliminate at the
  source by having the loader store `lfp.samples` as a clean 1-D array.
- **Hardcoded 50 Hz pos coupling** (`/50`, `eegFreq/50`, `50*60`) appears across pxd, FFT, LFP and
  FieldShuffle — parametrise uniformly on `Recording.pos_sample_rate`.
- **Drop all `parfor`/`gcp`/`goPar`** parallelism during the literal port (numerically irrelevant);
  reintroduce via `joblib`/`multiprocessing` later only where profiling shows it matters.
