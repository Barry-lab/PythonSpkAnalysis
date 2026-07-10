# barrelab

Spatial-coding analysis for the Barry lab, in Python. A modern, validated
successor to the [`Universal_Matlab`](https://github.com/Barry-lab/Universal_Matlab)
toolbox that reads **historic Axona/DACQ tetrode data today** and is built to
grow into **Neuropixels (SpikeGLX & Open Ephys)** without rewriting a single
line of analysis code.

> **Status: early / proof-of-concept (v0.1).** The DACQ loader, spatial binning,
> rate maps and Skaggs spatial information are implemented and validated against
> the lab's MATLAB output (see [Validation](#validation)). Neuropixels loading,
> gridness, head-direction and LFP analyses are on the [roadmap](#roadmap).

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

Neuropixels and NWB support are optional extras (`.[neuropixels]`, `.[nwb]`) so a
DACQ-only user stays lightweight.

## Quickstart

```python
import barrelab as bl

rec = bl.load_dacq("/data/r1838_10-12-14/83812141")   # any file of the trial
print(rec)                       # Recording(name='83812141', n_pos=60050, electrodes=[1..8])

rmap = bl.rate_map(rec, electrode=2, cluster=3)        # 2-D firing rate map (Hz)
bits_per_spike, bits_per_s = bl.spatial_information(rmap, bl.occupancy_map(rec))
```

See [`examples/quickstart.py`](examples/quickstart.py).

## Layout

```
src/barrelab/
  core/recording.py     the common Recording / Electrode objects
  io/
    _axona_files.py      byte-level Axona header + binary reader
    _position.py         DACQ position post-processing (swap/speed/interp/smooth)
    axona.py             load_dacq(): the DACQ front-end loader
  analysis/
    maps.py              binning + rate maps (histnd, make_smooth_ratemap)
    spatial.py           Skaggs spatial information + adaptive binning
  params.py              default parameters (from default_vars.m)
tests/                   pytest suite (synthetic + golden-master)
validation/              numerical validation against MATLAB, + fixture generator
```

## Validation

Each ported function is checked against the lab's MATLAB output. The example
trial ships with a reference `.mat` (`gridCell1838_t2c3`) holding the fully
loaded `data` struct and a golden rate map, so we can validate stage by stage:

```bash
python3 validation/validate_against_matlab.py     # 15/15 checks pass
```

Bit-exact (machine precision) against MATLAB: the `.pos`/tetrode/`.cut` readers,
the position pipeline (in `midpoint` mode), the spike→position mapping, cluster
selection, and the spatial binning (also cross-checked brute-force).

Two steps intentionally do **not** bit-match this particular fixture because the
fixture is older than the current code — see [`FINDINGS.md`](FINDINGS.md). To
lock these down, regenerate fresh fixtures from the current MATLAB with
[`validation/generate_fixtures.m`](validation/generate_fixtures.m) and commit them.

## Roadmap

- **v0.1 (now)** — DACQ loader, `Recording`, rate maps, Skaggs info, MATLAB validation harness.
- **v0.2** — gridness/SAC, head-direction tuning, border/speed scores, plotting; docs site.
- **v0.3** — Neuropixels front end (SpikeGLX + Open Ephys via SpikeInterface + Kilosort4).
- **v0.4 → v1.0** — LFP/theta/phase precession, Bayesian decoding, NWB export, plugin loader API.

## License

MIT — see [LICENSE](LICENSE).
