# barrelab

Spatial-coding analysis for the Barry lab, in Python. A modern, validated
successor to the [`Universal_Matlab`](https://github.com/Barry-lab/Universal_Matlab)
toolbox that reads **historic Axona/DACQ tetrode data today** and is built to
grow into **Neuropixels (SpikeGLX & Open Ephys)** without rewriting a single
line of analysis code. If you are new to the lab and coming from MATLAB, this is
where to start.

> **Status: early / proof-of-concept (v0.1).** The DACQ loader, spatial binning,
> rate maps and Skaggs spatial information are implemented and validated against
> the lab's MATLAB output (see [Validation](#validation)). Neuropixels loading,
> gridness, head-direction and LFP analyses are on the [roadmap](#roadmap).

## One glance

```python
import barrelab as bl

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

A loader's only job is to fill a [`Recording`](src/barrelab/core/recording.py).
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
src/barrelab/
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
pytest                                            # 17 tests
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
- [`FINDINGS.md`](FINDINGS.md) — validation notes and the two intentional MATLAB differences.

## Roadmap

- **v0.1 (now)** — DACQ loader, `Recording`, rate maps, Skaggs info, MATLAB validation harness.
- **v0.2** — gridness/SAC, head-direction tuning, border/speed scores, plotting.
- **v0.3** — Neuropixels front end (SpikeGLX + Open Ephys via SpikeInterface + Kilosort4).
- **v0.4 → v1.0** — LFP/theta/phase precession, Bayesian decoding, NWB export, plugin loader API.

## License

MIT — see [LICENSE](LICENSE).
