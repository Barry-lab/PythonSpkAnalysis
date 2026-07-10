# barrylab

Spatial-coding analysis for the Barry lab, in Python. A modern, validated
successor to the [`Universal_Matlab`](https://github.com/Barry-lab/Universal_Matlab)
toolbox that reads **historic Axona/DACQ tetrode data today** and is built to
grow into **Neuropixels (SpikeGLX & Open Ephys)** without rewriting a single
line of analysis code. If you are new to the lab and coming from MATLAB, this is
where to start.

> **This is an AI-assisted port of `Universal_Matlab`, not yet independently
> checked by a human. It needs checking.** Every ported function is validated
> numerically against MATLAB output where a reference exists (bit-exact
> fixture, synthetic ground truth, or a textbook/analytic oracle — see
> [Validation](#validation) and [`FINDINGS.md`](FINDINGS.md)), and any MATLAB
> bugs/quirks found along the way are documented there rather than silently
> "fixed" — but no lab member has yet reviewed the port line-by-line.
> `Universal_Matlab` remains the reference implementation and should be
> treated as ground truth until someone has been through this repo properly.
> If you use this package for real analysis, please sanity-check your results
> against the MATLAB toolbox and open an issue if anything looks off. This
> section will be updated as the port progresses and gets reviewed.
>
> **Status: in progress.** DACQ loading, spatial binning, rate maps, Skaggs
> info, grid/gridness, head direction, place-field stats, Bayesian decoding,
> phase precession, grid-ellipse correction, LFP loading + power spectrum, and
> watershed field segmentation are implemented and validated (see
> [Validation](#validation)). Neuropixels loading and several remaining
> analyses are on the [roadmap](#roadmap); see
> [`PORTING_PLAN.md`](PORTING_PLAN.md) for the full current status.

## One glance

```python
import barrylab as bl

rec = bl.load_dacq("/data/r1838_10-12-14/83812141")   # any file of the trial
rmap = bl.rate_map(rec, electrode=2, cluster=3)        # 2-D firing rate map (Hz)
occ = bl.occupancy_map(rec)                            # matching dwell-time map
bits_per_spike, bits_per_s = bl.spatial_information(rmap, occ)
```

See [`examples/quickstart.py`](examples/quickstart.py) for a runnable version
that also plots the trajectory and rate map.

## The idea: one backend, many frontends

```
RAW FILES            FRONT-END LOADERS          COMMON OBJECT        ANALYSIS
------------------   ------------------------   ------------------   ------------------
DACQ / Axona     ->  load_dacq()            \
Neuropixels/SpikeGLX  (SpikeInterface)       >  Recording        ->  rate maps, Skaggs
Neuropixels/OpenEphys (SpikeInterface)      /   (spikes, xy,          info, gridness,
                                                 speed, dir, lfp)     HD tuning, decoding
```

A loader's only job is to fill a [`Recording`](src/barrylab/core/recording.py).
Every analysis consumes that object and nothing else, so a rate map does not
care whether the spikes came from a 2010 tetrode or a 2026 Neuropixel.

## Install

```bash
# from a clone, with uv (recommended)
uv pip install -e ".[dev,plot]"
# or plain pip
pip install -e ".[dev,plot]"
```

Requires Python >= 3.9. Neuropixels and NWB support are optional extras
(`.[neuropixels]`, `.[nwb]`) so a DACQ-only user stays lightweight.

## Layout

```
src/barrylab/
  core/recording.py     the common Recording / Electrode objects
  io/
    _axona_files.py      byte-level Axona header + binary reader
    _position.py         DACQ position post-processing (swap/speed/interp/smooth)
    axona.py             load_dacq(): the DACQ front-end loader
  analysis/
    maps.py              binning + rate maps (histnd, smooth_ratemap)
    spatial.py           Skaggs spatial information + adaptive binning
  params.py              default parameters (from default_vars.m)
tests/                   pytest suite (synthetic + golden-master)
validation/              numerical validation against MATLAB, + fixture generator
```

## Validation

We follow a **golden-master** philosophy: each ported function is checked
numerically against the lab's MATLAB output, not merely for internal
consistency. The example trial ships with a reference `.mat` so we can validate
stage by stage.

```bash
python3 validation/validate_against_matlab.py     # 15/15 checks pass
pytest                                            # 106 tests
```

Two steps intentionally do **not** bit-match this particular fixture because it
predates the current code — both are documented PI decisions in
[`FINDINGS.md`](FINDINGS.md). Regenerate fresh fixtures from the current MATLAB
with [`validation/generate_fixtures.m`](validation/generate_fixtures.m).

## Documentation

- [`docs/getting_started.md`](docs/getting_started.md) — install, load a trial, make your first rate map.
- [`docs/architecture.md`](docs/architecture.md) — the three-layer design and the `Recording` object.
- [`docs/matlab_migration.md`](docs/matlab_migration.md) — where each `Universal_Matlab` function lands.
- [`docs/adding_a_loader.md`](docs/adding_a_loader.md) — how to write a new front-end loader.
- [`CONTRIBUTING.md`](CONTRIBUTING.md) — dev setup, tests, and coding conventions.
- [`FINDINGS.md`](FINDINGS.md) — validation notes, cross-function dependency gotchas, and known MATLAB quirks.
- [`PORTING_PLAN.md`](PORTING_PLAN.md) — the master plan for porting everything still remaining, module by module.

## Roadmap

Full detail and per-function status in [`PORTING_PLAN.md`](PORTING_PLAN.md).

- **v0.1** — DACQ loader, `Recording`, rate maps, Skaggs info, MATLAB validation harness. ✓
- **v0.2** — spatial autocorrelogram + gridness ✓, head-direction tuning ✓, place-field/coherence/correlation stats ✓.
- **v0.4** — Bayesian decoding (naive ML + full Bayes w/ movement prior) ✓, phase precession ✓.
- **v0.5** — grid ellipse correction ✓; 2-D FFT grid analysis + cross-SAC still to come.
- **v0.6** — LFP/EEG loading + power spectrum ✓; instantaneous freq/power and ripple detection still to come.
- **v0.7 (in progress)** — watershed field segmentation ✓ (`barrylab[fieldshuffle]`); the field-reallocation shuffle itself and `identify_grid_cell` orchestration remain.
- **v0.3 / v1.0** — Neuropixels front end (SpikeGLX + Open Ephys via SpikeInterface + Kilosort4), NWB export, plugin loader API, human review of the whole port.

## License

MIT — see [LICENSE](LICENSE).
