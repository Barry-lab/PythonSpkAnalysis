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
| `xPearson.m` (`Analysis_General`) | `spatial_autocorrelogram(rate_map)` / `cross_correlogram(a, b)` | public — `xPearson` is the general Pearson-correlation-at-all-lags primitive (one arg = autocorrelation = the SAC; two args = cross-correlation). These two wrappers are its autocorr and cross-corr uses. |
| `multiGridness.m` | `gridness(sac)` → `(standard, expanding)` | public |
| `sacProps.m` | `grid_properties(sac)` | public — scale/orientation/peaks |
| `sacProps90deg.m` | `analysis.grid.square_grid_props(sac)` | public — 4-fold (square) grid variant of `grid_properties` |
| (convenience) | `grid_score(recording, electrode, cluster)` → `(gridness, scale_cm)` | public — one call: rate map → SAC → standard gridness + scale |
| `fs_fastGridness.m` | `analysis.grid.fast_gridness(sac)` → `(standard, expanding, scale)` | public — NOT bit-comparable to `gridness`/`grid_properties`, see its docstring |
| `xPearson.m` (cross mode) | `analysis.grid.cross_correlogram(rm_a, rm_b)` | public — cross-SAC between two rate maps |
| `xSacProps.m` | `analysis.grid_phase.cross_sac_props(xsac)` | public — see docstring: a likely real MATLAB fragility on realistic data (check `peak_corr`), see `FINDINGS.md` |
| `fftVars.m` | `analysis.spatial_fft.FFTParams` | public — parameter dataclass; `max_rad` is a derived property |
| `get_2d_fft.m` | `analysis.spatial_fft.power_spectrum_2d(rate_map, fft_size, smth_sig, max_rad)` | public — mean-normalised, zero-padded 2-D FFT power spectrum |
| `get_polar_fft.m` | `analysis.spatial_fft.to_polar_fft(ps2d, max_rad, th_bin, r_bin, sigma)` | public — polar (angle/radius) reorganisation of a 2-D FFT power spectrum |
| `main_spatial_fft2.m` | `analysis.spatial_fft.grid_fft_components(recording, electrode, cluster, ...)` | public — Recording-level redesign like `decode_crunch_data`, not bit-exact; `main_spatial_fft.m` (v1) superseded, not ported |
| `cb_circMean.m` | `circular_mean(angles, weights, wrap=...)` | public — circular mean; see `wrap` note below |
| `dir_kldiv.m` / `kldiv.m` | `directionality(tuning_curve)` | public — directionality index |
| `circ_FTT.m` | `directional_fft(tuning_curve)` | public — multimodal/hexadirectional |
| `led2v1_dir_offset.m` | `led2v1_dir_offset(led2_dir, led1_dir)` | public — 2-LED/1-LED heading bias |
| (Rayleigh vector) | `resultant_vector(tuning_curve)` → `(r, preferred_deg)` | public — mean resultant length + preferred direction |
| `bin_data('dir')` | `directional_occupancy(recording)` | public — dwell time per direction bin |
| (HD tuning curve) | `head_direction_map` / `head_direction_stats` | public — bins by direction |
| `pxd.m` | `pxd(spikes, times)` | public — place x direction ML fit |
| `wrapper_pxd.m` | `place_by_direction(recording, electrode, cluster)` | public — bins + fits + smooths |
| `histnd.m` (general N-D case) | `_binning.histnd(data, ranges, bin_size)` | internal — used by `pxd` |
| `pcFieldStats.m` | `field_stats(rate_map, min_field_bins, min_rate)` | public — place-field size/rate |
| `spatialCoherance.m` | `spatial_coherence(rate_map)` | public — smoothness/signal quality |
| `spatialCorrelation.m` | `spatial_correlation(rm1, rm2, stack_maps=True)` | public — remapping/stability |
| `nan_corr.m` | `nan_corr(var1, var2, corr_type="pearson")` | public — NaN/Inf-safe correlation |
| `get_true_PPM_boundingBox.m` | `true_ppm_bounding_box(occupancy, bin_px, env_m)` | public — infer true PPM |
| `decode_calcPost.m` | `decode_calc_post(obs_spk, rate_maps, t_bin)` | public — naive ML Poisson decoder |
| `decode_getSpkCnt.m` | `decode_get_spike_counts(...)` | public — windowed spike counts |
| `decode_processPost.m` | `decode_process_post(post, recording, bin_px, win_lim)` | public — ML decode + error |
| `decode_calcBayesPost.m` | `decode_calc_bayes_post(...)` | public — full Bayes + movement prior |
| `decode_crunchData.m` | `decode_crunch_data(recording, cells, ...)` | public — Recording-level orchestrator (redesign, not bit-exact) |
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
| `eeg_instant_freq_power.m` | `instant_freq_power(eeg, filter_range, eeg_freq, pos_sample_rate)` | public — bandpass + Hilbert instantaneous phase/amp/freq |
| `detect_ripples.m` | `detect_ripples(eeg, eeg_freq, rip_freq, rip_thresh)` | public — ripple/SWR detection; **default `rip_thresh` is 5, matching the MATLAB code (its own docstring wrongly claims 4 — see `FINDINGS.md`)** |
| `fs_masterAnalysis.m` watershed step + `fs_fixWatershed.m` | `field_shuffle.watershed_fields(rate_map, occupancy)` | public — needs `barrylab[fieldshuffle]` |
| `fs_binMembLoc.m` | `field_shuffle.field_bin_locations(labels, rate_map)` | public — needs `barrylab[fieldshuffle]` |
| `fs_shufFields.m` | `field_shuffle.shuffle_fields(shape, field_bin_loc, static_bin, rng=...)` | public — involves genuine randomness, see docstring/`FINDINGS.md` |
| `fs_masterAnalysis.m` | `field_shuffle.field_shuffle_significance(recording, electrode, cluster, ...)` | public — Recording-level redesign like `decode_crunch_data`, not bit-exact |
| `identify_grid_cell.m` | `field_shuffle.is_grid_cell(recording, electrode, cluster, ...)` | public — see docstring/`FINDINGS.md` for 2 genuine MATLAB bugs found, not reproduced |
| `read_INP_file.m` | `io.axona.read_inp(path)` | public — digital-timestamp (button/DIO) events; validated structurally, no example `.inp` ships with `Example_data` |
| `read_stm_file.m` | `io.axona.read_stm(path)` | public — stimulator-timestamp events; validated structurally, no example `.stm` ships with `Example_data` |
| `write_cut_file.m` | `io.axona.write_cut(spike_assignment, dir_path, file_stem, tetrode_no)` | public — writes a Tint-compatible `.cut`; fixes 2 genuine MATLAB bugs, see `FINDINGS.md` |

Everything marked *internal* runs inside `load_dacq`; you do not call it directly,
but it is exposed for validation and testing.

**A dependency gotcha worth knowing about `circular_mean`.** `circular_mean()`
defaults to `wrap="0_2pi"`, matching `cb_circMean.m` exactly. But if you're
porting or comparing against a MATLAB function that instead calls the toolbox
`circ_mean` (as `led2v1_dir_offset.m` does — see `FINDINGS.md` §3), that
function's angles come back in `(-π, π]`, not `[0, 2π)`. The mean angle is
identical either way; only the range differs. Pass `wrap="-pi_pi"` when you
need the signed convention. See the per-function validation notes and
dependency gotchas in [`FINDINGS.md`](../FINDINGS.md).

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
