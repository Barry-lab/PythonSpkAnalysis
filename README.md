# barrylab

Spatial-coding analysis for tetrode (Axona/DACQ) and, soon, Neuropixels
recordings — rate maps, grid cells, head direction, Bayesian decoding, replay
scoring, phase precession, and LFP/theta analysis, all built around one
common `Recording` object so the same analysis code runs unchanged regardless
of what acquisition system the spikes came from.

This is a Python port of the Barry lab's MATLAB toolbox,
[`Universal_Matlab`](https://github.com/Barry-lab/Universal_Matlab), done
with an AI-assisted, source-read-and-validated process rather than by hand —
see [Provenance & validation](#provenance--validation) below before relying
on it for a paper. `Universal_Matlab` remains the reference implementation
until a lab member has reviewed this port properly.

## Install

```bash
pip install -e ".[dev,plot]"          # or: uv pip install -e ".[dev,plot]"
```

Requires Python ≥ 3.9. Everything below works out of the box for DACQ data;
`[fieldshuffle]` (scikit-image, for the watershed-based field-shuffle test)
and `[neuropixels]`/`[nwb]` (for the not-yet-built Neuropixels/NWB front
ends) are separate extras so a DACQ-only install stays lightweight.

## Quickstart

```python
import barrylab as bl

rec  = bl.load_dacq("/data/r1838_10-12-14/83812141")     # any file of the trial
rmap = bl.rate_map(rec, electrode=2, cluster=3)            # 2-D firing rate map (Hz)
occ  = bl.occupancy_map(rec)                               # matching dwell-time map
bits_per_spike, bits_per_s = bl.spatial_information(rmap, occ)
```

`load_dacq` reads the whole trial (position, all tetrodes, cluster cuts,
optionally LFP) into a single `Recording` — every analysis function below
just takes that object (or a rate map/array derived from it) and a couple of
parameters. See [`docs/getting_started.md`](docs/getting_started.md) for a
walked-through tour (install → load a trial → first rate map → a figure) and
[`examples/quickstart.py`](examples/quickstart.py) for the runnable script
version.

## What's here

Grouped the same way the original MATLAB toolbox is (`Analysis_Grid`,
`Analysis_Dir`, `Decoding`, `Analysis_Freq`, ...) so it's easy to find the
Python equivalent of a function you already know. Every snippet below is
runnable once you have a `rec` from `load_dacq` (and `import numpy as np` for
the couple that need it) — see [`docs/matlab_migration.md`](docs/matlab_migration.md)
for the complete function-by-function map to the MATLAB originals.

### Rate maps & spatial information

The foundation everything else builds on: bin occupancy and spikes over
space, smooth, and quantify how spatially informative a cell's firing is.

```python
rmap = bl.rate_map(rec, 2, 3, smoothing="gaussian", smooth_size=2)  # sigma in bins
bits_per_spike, bits_per_s = bl.spatial_information(rmap, bl.occupancy_map(rec))
```

### Grid cells

A spatial autocorrelogram (SAC), its gridness (how hexagonally-regular the
firing pattern is), grid scale/orientation, elliptical-distortion correction,
and cross-cell phase comparison for two cells from the same grid module.

```python
sac = bl.spatial_autocorrelogram(rmap)
standard_g, expanding_g = bl.gridness(sac)          # two classic gridness measures
props = bl.grid_properties(sac)                     # scale, orientation, peak geometry

# correct for elliptical distortion before recomputing gridness
ellipse = bl.fit_grid_ellipse(sac)                  # note: orientation is in RADIANS
corrected_sac = bl.regularise_elliptic_sac(
    sac, ellipse["ab_scale"], orient_deg=np.degrees(ellipse["orientation"])
)

# relative spatial phase between two cells of the same grid module
xsac = bl.cross_correlogram(rmap_a, rmap_b)
phase = bl.cross_sac_props(xsac)                    # phase_mag, phase_ang, ...
```

`bl.fast_gridness` is a second, faster gridness implementation (standard +
expanding + scale in one call) used by the field-shuffle significance test
below — it is a genuinely different algorithm from `gridness`, not just an
optimisation, so don't expect the two to agree bit-for-bit (see its
docstring).

### Head direction

Directional tuning curves, resultant vector, directionality index, and
hexadirectional (grid-cell) modulation of firing by heading.

```python
hd_map = bl.head_direction_map(rec, 2, 3)
stats = bl.head_direction_stats(rec, 2, 3)          # tuning curve + resultant length/direction/directionality in one dict
directionality = bl.directionality(hd_map)          # KL-divergence-based index, on its own
hex_mod = bl.directional_fft(hd_map)                # 6-fold (and other) modulation
```

### Place-field statistics

Field size/rate/count, spatial coherence, and map-to-map correlation (for
remapping/stability analyses).

```python
fields = bl.field_stats(rmap, min_field_bins=9, min_rate=1.0)
coherence = bl.spatial_coherence(rmap)
stability = bl.spatial_correlation(rmap_session1, rmap_session2)
```

### Place × direction

A joint place-by-direction firing model, fit by alternating maximum
likelihood — useful when a cell's spatial and directional tuning aren't
independent.

```python
place_map, direction_map = bl.place_by_direction(rec, electrode=2, cluster=3)
```

### Bayesian population decoding

Decode position from a population's spike counts — a naive maximum-likelihood
Poisson decoder, a full Bayesian version with a movement prior, and a
Recording-level convenience wrapper that builds rate maps on a training split
and decodes a held-out test split.

```python
result = bl.decode_crunch_data(rec, cells=[(2, 3), (2, 1), (3, 4)], t_bin=0.25)
print(result["err"].mean(), "cm mean decode error")
```

Lower-level pieces (`bl.decode_calc_post`, `bl.decode_calc_bayes_post`,
`bl.decode_get_spike_counts`, `bl.decode_process_post`) are all public
too if you want to build your own pipeline.

### Replay / trajectory scoring

Given a decoded posterior over time (e.g. during a sharp-wave ripple), find
the best constant-speed straight-line trajectory through it — a common way
to score candidate replay events.

```python
best_res, best_speed, best_y_intercept, line_mask = bl.best_line(
    posterior_matrix, t_bin_sz=0.005, s_bin_sz=2, y_range=14
)
```

### Phase precession

Circular-linear regression relating position-within-field to theta phase,
with permutation-based significance.

```python
fit = bl.circular_linear_regression(position_in_field, theta_phase, k=1000)
print(fit["slope_opt"], fit["corr"], fit["psim"])   # rad/unit-x, circular-circular R, permutation p-value
```

### LFP / theta

Load raw LFP alongside position/spikes (`load_dacq(..., load_lfp=True)`),
then get an FFT-based power spectrum with automatic theta-band peak
detection.

```python
rec = bl.load_dacq(path, load_lfp=True)
eeg = rec.lfp[1]["samples"]                          # channel 1's raw trace
spectrum = bl.power_spectrum(eeg, eeg_freq=rec.lfp[1]["sample_rate"])
print(spectrum["max_freq"], "Hz theta peak, s2n =", spectrum["s2n"])
```

### Grid-cell significance testing (field shuffle)

Is a cell's gridness higher than chance? Segment its rate map into firing
fields (watershed), then compare its real gridness against temporal-shuffle
and field-shuffle null distributions. Needs the `[fieldshuffle]` extra
(`pip install "barrylab[fieldshuffle]"`), so it lives under
`barrylab.analysis.field_shuffle` rather than the top-level namespace:

```python
from barrylab.analysis.field_shuffle import is_grid_cell

result = is_grid_cell(rec, electrode=2, cluster=3, n_shuffles=100)
print(result["is_grid"])   # (4, 2) bool: 4 gridness variants x [temporal, field] shuffle
```

## The `Recording` object

```
RAW FILES              FRONT-END LOADERS          COMMON OBJECT        ANALYSIS
--------------------   ------------------------   ------------------   ------------------
DACQ / Axona       ->  load_dacq()            \
Neuropixels (sorted +  load_neuropixels()      >  Recording        ->  everything above
  a position track)    [experimental]         /   (spikes, xy,
                                                   speed, dir, lfp)
```

A loader's only job is to fill a
[`Recording`](src/barrylab/core/recording.py) — a plain container you can
poke at directly (`rec.xy`, `rec.speed`, `rec.direction`, `rec.pos_times`,
`rec.electrode(n).spike_pos_index(cluster)`, ...). Every analysis function
consumes that object (or arrays derived from it) and nothing else, so the same
analysis code runs on DACQ and Neuropixels data alike. See
[`docs/architecture.md`](docs/architecture.md) for the full design, and
[`docs/adding_a_loader.md`](docs/adding_a_loader.md) if you're writing a new
front end.

### Neuropixels (experimental)

> ⚠️ **The Neuropixels front end is experimental and far less tested than the
> DACQ loader.** It has no golden-master validation and is exercised only
> against synthetic/mock data (open Neuropixels spatial datasets are huge and
> need GPU spike-sorting this package doesn't do). Treat results with care,
> cross-check against a method you trust, and **please improve/fix it** — it's
> a starting point, not a finished pipeline. Bug reports and PRs welcome.

`barrylab` doesn't spike-sort or read raw Neuropixels files — that's what
Kilosort + [SpikeInterface](https://spikeinterface.readthedocs.io) are for.
`load_neuropixels` takes the **already-sorted** output plus a **separate
position track** (a Neuropixels probe records neural data only — position
comes from your VR rig / camera / DeepLabCut, on its own clock, which you must
have already aligned to the neural clock's `t=0`):

```python
import barrylab as bl

# `sorting` = a Kilosort/phy output folder (read via SpikeInterface), or any
# SpikeInterface Sorting object. `positions` = your (n_pos, 2) behavioural
# track, sampled at `pos_sample_rate`.
rec = bl.load_neuropixels(
    sorting="/data/mouse42/kilosort_output",
    positions=xy,                 # (n_pos, 2) x, y from your tracker
    pos_sample_rate=50.0,         # Hz of the position track
    pixels_per_metre=100.0,       # spatial scale (use 100 if xy is already in cm)
)

# from here it's identical to DACQ: every analysis just takes the Recording
rmap = bl.rate_map(rec, electrode=1, cluster=7)   # cluster = a sorted unit id
gridness, _ = bl.gridness(bl.spatial_autocorrelogram(rmap))
```

Needs the `[neuropixels]` extra for the SpikeInterface reader
(`pip install "barrylab[neuropixels]"`); if you already have a `Sorting`
object in hand, no extra install is needed. Good public datasets to try it on
include the [One-shot entorhinal maps](https://data.mendeley.com/datasets/rgtk6jygjc/1)
Neuropixels VR + open-field recordings.

## Documentation

- [`docs/getting_started.md`](docs/getting_started.md) — install, load a trial, make your first rate map and figure.
- [`docs/architecture.md`](docs/architecture.md) — the three-layer design (`io` → `core` → `analysis`) and the `Recording` object.
- [`docs/matlab_migration.md`](docs/matlab_migration.md) — the complete `Universal_Matlab` function → `barrylab` function map, for anyone migrating existing analysis scripts.
- [`docs/adding_a_loader.md`](docs/adding_a_loader.md) — how to write a new front-end loader (e.g. for a new acquisition system).
- [`CONTRIBUTING.md`](CONTRIBUTING.md) — dev setup, running tests, coding conventions.

## Provenance & validation

Every ported function was checked numerically before being considered
done — against a bit-exact MATLAB fixture where one exists, otherwise a
synthetic case with known ground truth, an independent textbook/analytic
oracle, or cross-consistency against an already-validated function (real
data alone was never treated as sufficient). Any MATLAB bugs or quirks found
along the way are documented rather than silently fixed, unless replicating
one would leave the Python function simply broken. Full detail — the exact
validation used for every function, and every MATLAB bug found — is in
[`FINDINGS.md`](FINDINGS.md); remaining/deferred pieces are tracked in
[`PORTING_PLAN.md`](PORTING_PLAN.md).

**No lab member has independently reviewed this port line-by-line yet.** If
you use it for real analysis, sanity-check results that matter against the
MATLAB toolbox, and open an issue if anything looks off.

## License

MIT — see [LICENSE](LICENSE).
