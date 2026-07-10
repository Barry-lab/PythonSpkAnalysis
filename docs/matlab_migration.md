# Migrating from `Universal_Matlab`

If you know the lab's MATLAB toolbox, you already know `barrylab`. The analyses
are the same; only the packaging changed. `barrylab` is a modern, validated
Python successor that reads the same historic Axona/DACQ data and is built to
grow into Neuropixels without rewriting your analysis code. Where a function was
ported, its numerics were checked against MATLAB output rather than merely for
internal consistency — the loader chain and spatial binning are validated
**bit-exact** (see `FINDINGS.md`).

The big conceptual shift: instead of passing loose structs and arrays between
functions, you load a trial once into a single `Recording` object and hand that
object to every analysis. A rate map does not care whether the spikes came from a
2010 tetrode or a 2026 Neuropixel — it only sees a `Recording`.

## Function map

| `Universal_Matlab` | `barrylab` | Notes |
|---|---|---|
| `read_DACQ.m` | `load_dacq(path)` | public entry point |
| `read_pos_file.m` | `io.axona._read_pos` | internal (via `load_dacq`) |
| `read_tetrode_file.m` | `io.axona._read_tetrode` | internal |
| `read_cut_file.m` | `io.axona._read_cut` | internal |
| `postprocess_pos_data.m` | `io._position.postprocess_pos_data` | internal |
| `led_swap_filter` / `led_speed_filter` | `io._position` | internal |
| `key_value.m`, `read_binary_data.m`, `find_word.m` | `io._axona_files` | internal |
| `bin_data.m`, `histnd.m` | `occupancy_map`, `spike_map` (`histnd` = `_histnd_2d`) | public maps |
| `make_smooth_ratemap.m` | `smooth_ratemap` (and `rate_map` wrapper) | public |
| `easy_make_ratemap.m` | `rate_map(recording, electrode, cluster)` | public |
| `skaggsInfo.m` | `spatial_information(rates, dwell_times)` | public |
| `skaggsAdaptiveBinning.m` | `adaptive_binning(dwell_times, spikes, alpha=200)` | public |
| `default_vars.m` | `barrylab.params` | subset ported so far |
| `xPearson.m` (`Analysis_General`) | `spatial_autocorrelogram(rate_map)` | public — the SAC |
| `multiGridness.m` | `gridness(sac)` → `(standard, expanding)` | public |
| `sacProps.m` | `grid_properties(sac)` | public — scale/orientation/peaks |
| `fs_fastGridness.m` | `analysis.grid.fast_gridness(sac)` → `(standard, expanding, scale)` | public — NOT bit-comparable to `gridness`/`grid_properties`, see its docstring |
| `xPearson.m` (cross mode) | `analysis.grid.cross_correlogram(rm_a, rm_b)` | public — cross-SAC between two rate maps |
| `xSacProps.m` | `analysis.grid_phase.cross_sac_props(xsac)` | public — see docstring: a likely real MATLAB fragility on realistic data (check `peak_corr`), see `FINDINGS.md` |
| `Analysis_Grid` (2-D FFT, field shuffle orchestration) | — | **not yet ported** (roadmap v0.5/v0.7) |
| `cb_circMean.m` | `circular_mean(angles, weights, wrap=...)` | public — circular mean; see `wrap` note below |
| `dir_kldiv.m` / `kldiv.m` | `directionality(tuning_curve)` | public — directionality index |
| `circ_FTT.m` | `directional_fft(tuning_curve)` | public — multimodal/hexadirectional |
| `led2v1_dir_offset.m` | `led2v1_dir_offset(led2_dir, led1_dir)` | public — 2-LED/1-LED heading bias |
| (HD tuning curve) | `head_direction_map` / `head_direction_stats` | public — bins by direction |
| `pxd.m` | `pxd(spikes, times)` | public — place x direction ML fit |
| `wrapper_pxd.m` | `place_by_direction(recording, electrode, cluster)` | public — bins + fits + smooths |
| `histnd.m` (general N-D case) | `_binning.histnd(data, ranges, bin_size)` | internal — used by `pxd` |
| `pcFieldStats.m` | `field_stats(rate_map, min_field_bins, min_rate)` | public — place-field size/rate |
| `spatialCoherance.m` | `spatial_coherence(rate_map)` | public — smoothness/signal quality |
| `spatialCorrelation.m` | `spatial_correlation(rm1, rm2, stack_maps=True)` | public — remapping/stability |
| `nan_corr.m` | `nan_corr(var1, var2, corr_type="pearson")` | public — NaN/Inf-safe correlation |
| `get_true_PPM_boundingBox.m` | `true_ppm_bounding_box(occupancy, bin_px, env_m)` | public — infer true PPM |
| `decode_calcPost.m` | `bayesian_posterior(obs_spk, rate_maps, t_bin)` | public — naive ML Poisson decoder |
| `decode_getSpkCnt.m` | `spike_counts_per_window(...)` | public — windowed spike counts |
| `decode_processPost.m` | `decode_maximum_likelihood(post, recording, bin_px, win_lim)` | public — ML decode + error |
| `decode_calcBayesPost.m` | `bayesian_posterior_with_prior(...)` | public — full Bayes + movement prior |
| `decode_crunchData.m` | `crunch_decode(recording, cells, ...)` | public — Recording-level orchestrator (redesign, not bit-exact) |
| `circ_regress_cb.m` | `circular_linear_regression(x, t, k, pp_dir)` | public — phase precession fit |
| `lineTraj_decode.m` | `line_trajectory.best_line(pmat, t_bin_sz, s_bin_sz, y_range)` | public — best-fit replay trajectory |
| `lineTraj_reget_lnBnd.m` | `line_trajectory.reget_line_band(...)` | public — rebuild a line-band mask for another pMat |
| `lineTraj_get_res.m` | `line_trajectory.line_residual(pmat, line_band)` | public — score a line-band mask against a pMat |
| `gridEllipse_fit.m` | `fit_grid_ellipse(sac)` | public — ellipse fit to 6 central peaks |
| `gridEllipse_correct.m` / `regularise_eliptic_grid.m` | `regularise_elliptic_sac(sac, ab_scale, orient_deg)` | public — de-ellipse a SAC |
| `brandonRegEcc.m` | `grid_ellipse.brandon_regularise(sac)` → `(rsac1, rsac2)` | public — a self-admittedly weaker alternative ellipse correction, see its docstring |
| `read_eeg_file.m` | `io.axona._read_eeg` | internal — 250Hz .eeg reader |
| `getegf.m` (raw read only, no decimation — see `FINDINGS.md`) | `io.axona._read_egf` | internal — 4.8kHz .egf reader |
| `eeg_convert_to_volts.m` | `io.axona.eeg_convert_to_volts(samples, channel, settings, kind)` | public — counts to volts |
| (LFP loading) | `load_dacq(..., load_lfp=True)` -> `Recording.lfp` | public — unlocks theta/LFP analysis |
| `eeg_power_spectra.m` | `power_spectrum(eeg, pos2use, ...)` | public — FFT power spectrum, theta peak/s2n |
| `eeg_instant_freq_power.m`, `detect_ripples.m` | — | **not yet ported** (`PORTING_PLAN.md` v0.6b, deferred) |
| `fs_masterAnalysis.m` watershed step + `fs_fixWatershed.m` | `field_shuffle.watershed_fields(rate_map, occupancy)` | public — needs `barrylab[fieldshuffle]` |
| `fs_binMembLoc.m` | `field_shuffle.field_bin_locations(labels, rate_map)` | public — needs `barrylab[fieldshuffle]` |
| `fs_shufFields.m`, `fs_masterAnalysis.m`/`identify_grid_cell.m` orchestration | — | **not yet ported** (`PORTING_PLAN.md` v0.7, deferred) |
| 2-D FFT grid analysis, cross-SAC | — | **not yet ported** (`PORTING_PLAN.md` v0.5) |

Everything marked *internal* runs inside `load_dacq`; you do not call it directly,
but it is exposed for validation and testing.

**A dependency gotcha worth knowing about `circular_mean`.** `circular_mean()`
defaults to `wrap="0_2pi"`, matching `cb_circMean.m` exactly. But if you're
porting or comparing against a MATLAB function that instead calls the toolbox
`circ_mean` (as `led2v1_dir_offset.m` does — see `FINDINGS.md` §3), that
function's angles come back in `(-π, π]`, not `[0, 2π)`. The mean angle is
identical either way; only the range differs. Pass `wrap="-pi_pi"` when you
need the signed convention. See the full porting plan and dependency notes in
[`PORTING_PLAN.md`](../PORTING_PLAN.md).

## How you used to / how you do it now

Load a trial, make a rate map for a cell, and read out its spatial information.

**MATLAB**

```matlab
data = read_DACQ('83812141');
rm   = easy_make_ratemap(data, 2, 3);         % tetrode 2, cluster 3
info = skaggsInfo(rm, bin_data(data, 'dwell'));
```

**barrylab**

```python
import barrylab as bl

rec  = bl.load_dacq("/data/r1838/83812141")
rm   = bl.rate_map(rec, electrode=2, cluster=3)     # 5-bin boxcar by default
bits_per_spike, bits_per_sec = bl.spatial_information(rm, bl.occupancy_map(rec))
```

`load_dacq` accepts any file of the trial or the shared basename; the extension is
stripped for you. `rate_map` defaults to a 5-bin boxcar (the lab standard) but
takes `smoothing="gaussian"`, `smoothing=None`, and `smooth_size=...` if you want
to override it. The `Recording` also carries `rec.xy`, `rec.direction`,
`rec.speed`, `rec.pos_times`, and `rec.electrodes`; reach a tetrode with
`rec.electrode(2)` and list its clusters with `.cluster_ids()`.

## Two behavioural differences vs the old example fixture

The shipped `gridCell1838_t2c3_14-12-10.mat` fixture predates the current MATLAB,
so two steps deliberately differ from it. Both are **PI decisions**, documented in
full in `FINDINGS.md`; the underlying binning is proven bit-exact, so these are
choices about convention, not correctness.

1. **Head-position weighting.** The current `postprocess_pos_data.m` combines the
   two LEDs with a *tracking-fraction weighted* mean (`head_pos_weighting=
   "tracking"`, the `barrylab` default), whereas the fixture used a plain
   **midpoint**. Pass `load_dacq(..., head_pos_weighting="midpoint")` to reproduce
   the legacy fixture exactly. For this trial the back LED is tracked only ~53% of
   the time, so the conventions can differ by up to ~1 bin.

2. **Rate-map smoothing.** `rate_map` defaults to the current 5-bin boxcar
   (`easy_make_ratemap.m`), which correlates with the fixture at r≈0.94. The
   fixture itself was made with a *precursor* smoother that empirically matches a
   Gaussian (σ≈2 bins, r≈0.98). Use `smoothing="gaussian"` if you need to match
   that older map; the residual is the fixture's smoothing, not the port.

The recommendation in `FINDINGS.md` is to regenerate fresh golden-master fixtures
from the current MATLAB (`validation/generate_fixtures.m`) so every ported
function has a current numerical test.

## Minor MATLAB issues we found

While porting we spotted a handful of small, harmless bugs in the original
toolbox — a `nargoout`/`nargout` typo in `bin_data.m`, a no-op `lower()` on a
numeric array in `make_smooth_ratemap.m`, a mis-named tetrode-count loop in
`read_DACQ.m`, and a cosmetic out-of-range warning threshold in `histnd.m`. None
affect the numerics you rely on, and each is handled or dropped in the port.
See section 2 of `FINDINGS.md` for the details.
